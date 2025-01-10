#!/usr/bin/env python
# ruff: noqa: I002

"""
Parse a dependencies specification, applying preprocessing selectors like conda-build.

Multiple dependency files can be merged together.
General example of dependencies.yaml:

    build:
    - {{ compiler('cxx') }}   #[not bootstrap]
    - cxx-compiler            #[bootstrap]
    - cmake                   #[prebuilt_cctbx]
    - setuptools
    host:
    - cctbx-base >=2024       #[prebuilt_cctbx]
    - pip
    - python
    - libcxx                  #[bootstrap and osx]
    run:
    - ffbidx                  #[linux]
    - numpy >=1.21.5,<2       #[bootstrap]
    - numpy                   #[not bootstrap]
    - tabulate
    - tqdm
    test:
    - dials-data >=2.4.72
    - pytest

Only the last "#[expression]" is stripped and evaluated, which allows
passing further selection expressions through for onward processing.

Expression values supported:
- Platform e.g. "osx", "win" or "linux"
- "bootstrap": Whether we are targeting bootstrap.py. This allows us
  to have advanced conda-build syntax (e.g. jinja templates) while
  ignoring these when generating the bootstrap dependency list.
- "prebuilt_cctbx": Whether we are building cctbx. This is true for
  bootstrap with CMake and for conda-build builds.
- Compound expressions with "and" e.g. "osx and bootstrap". This will
  work with "not"-expressions, but nothing more complex.
"""

import logging
import os
import re
import sys
from collections import namedtuple

VALID_SECTIONS = ["build", "host", "run", "test"]  # type: list[SectionName]
Dependency = namedtuple("Dependency", ["name", "version", "raw_line"])

try:
    from typing import Any, Literal, TypeAlias  # noqa: F401

    SectionName = Literal["build", "host", "run", "test"]  # type: TypeAlias
    Dependencies = dict[SectionName, list[Dependency]]  # type: TypeAlias

except ImportError:
    pass


re_selector = re.compile(r"# *\[([^#]+)]$")
re_pin = re.compile(r"""{{ *pin_compatible *\( *['"]([^'"]+)['"]""")


def _native_platform():
    # type: () -> Literal["osx", "win", "linux"]
    """Gets the native platform name for selection purposes"""
    if sys.platform == "darwin":
        return "osx"
    elif os.name == "nt":
        return "win"
    elif sys.platform.startswith("linux"):
        return "linux"


def _split_dependency_line(line):
    """Split a single line into (name, version, raw_line) parts"""
    # type: (str) -> Dependency

    # Lines that are templated get ignored here
    if "{" in line:
        return Dependency(None, None, line)
    pending = line
    # Strip off the comment/selector
    if "#" in line:
        pending = pending[: pending.index("#")].strip()
    # If we have a version spec and no space, this is an error
    if " " not in pending and (set(pending) & set("><=!")):
        raise RuntimeError(
            "Error: Versioned requirement '%s' has no space" % (pending,)
        )
    vers = None
    if " " in pending:
        pending, vers = pending.split(" ", 1)
    return Dependency(pending, vers, line)


def _merge_dependency_lists(source, merge_into):
    # type: (list[Dependency], list[Dependency]) -> None
    """
    Merge two lists of dependencies into one unified list.

    This will replace unversioned dependencies with versioned
    dependencies, merge dependencies with identical versions, and
    leave in place depenencies with versions specified.

    Lines from the source list that don't have a dependency name
    will be added as long as they don't have a duplicate line in the
    target list.
    """
    indices = {x[0]: i for i, x in enumerate(merge_into)}
    for pkg, ver, line in source:
        if pkg is None:
            # Lines that don't define a package always get added, as long
            # as we don't have an identical line already.
            if not any(x.raw_line == line for x in merge_into):
                merge_into.append(Dependency(pkg, ver, line))
        elif pkg in indices:
            # This already exists in the target. Should we replace it?
            other_ver = merge_into[indices[pkg]][1]
            if not other_ver and ver:
                logging.debug(
                    "Merging '{}' over {}".format(line, merge_into[indices[pkg]])
                )
                merge_into[indices[pkg]] = Dependency(pkg, ver, line)
            elif other_ver and ver and ver != other_ver:
                raise RuntimeError(
                    "Cannot merge conflicting requirements for %s: '%s' and '%s' - only know how to merge if these are the same, or one is unbound."
                    % (pkg, ver, other_ver)
                )
        else:
            merge_into.append(Dependency(pkg, ver, line))
            indices[pkg] = len(merge_into) - 1


class DependencySelectorParser(object):
    """
    Parse simple conda-build selectors syntax, with optional variables.

    Supported:
    - Variables linux, osx, win, in addition to anything passed into __init__
    - Variable inversion e.g. "not osx"
    - Basic "And" combinations e.g. "bootstrap and not osx"
    """

    def __init__(self, **kwargs):
        self._vars = dict(kwargs)
        if kwargs.get("platform", None) is None:
            kwargs["platform"] = _native_platform()
        self._vars.update(
            {
                "osx": kwargs["platform"] == "osx",
                "linux": kwargs["platform"] == "linux",
                "win": kwargs["platform"] == "win",
            }
        )

    def _parse_expression(self, fragment, pos=0):
        # type: (str, int) -> bool
        """Recursively parse an expression or fragment of an expression."""
        if fragment in self._vars:
            return self._vars[fragment]
        if " and " in fragment:
            left, right = fragment.split(" and ", 1)
            return self._parse_expression(left, pos) and self._parse_expression(
                right, pos + fragment.index(" and ")
            )
        if fragment.startswith("not "):
            return not self._parse_expression(fragment[4:].strip(), pos + 4)
        raise ValueError("Could not parse selector fragment '" + fragment + "'")

    def preprocess(self, data):
        # type: (str) -> str
        """Apply preprocessing selectors to raw file data"""
        output_lines = []
        for line in data.splitlines():
            match = re_selector.search(line)

            if match:
                if self._parse_expression(match.group(1)):
                    output_lines.append(line)
            elif re_pin.search(line):
                # Ignore pin_compatible dependencies
                continue
            else:
                output_lines.append(line)
        return "\n".join(output_lines)

    def parse_file(self, filename):
        # type: (str) -> Dependencies
        """
        Parse a dependency file into a structured dictionary.

        The dictionary has structure:
        {
            "section": [
                ("dependency_name", "dependency_version", "raw_line"),
                ...
            ]
        }
        """
        with open(filename, "rt") as f:
            data = self.preprocess(f.read())
        output = {}  # type: Dependencies
        current_section = None  # type: SectionName | None
        for n, line in enumerate(data.splitlines()):
            if "#" in line:
                line = line[: line.index("#")]
            line = line.strip()
            if line.endswith(":"):
                new_section = line[:-1].strip()
                assert new_section in VALID_SECTIONS
                current_section = new_section
                output[current_section] = []
            elif line.startswith("-"):
                if not current_section:
                    raise RuntimeError(
                        "Error parsing "
                        + filename
                        + ":"
                        + str(n + 1)
                        + "; No current section on line '"
                        + line
                        + "'"
                    )
                assert current_section in VALID_SECTIONS
                req = _split_dependency_line(line[1:].strip())
                the_list = output.setdefault(current_section, [])
                the_list.append(req)
            else:
                if line:
                    raise RuntimeError(
                        "Error parsing "
                        + filename
                        + ":"
                        + str(n + 1)
                        + "; Uncategorised line '"
                        + line
                        + "'"
                    )
        return output

    def parse_files(self, filenames):
        # type: (list[str | os.PathLike]) -> Dependencies
        """Parse and merge multiple dependency files."""
        reqs = {}  # type: Dependencies
        for source in filenames:
            source_reqs = self.parse_file(str(source))
            # Now, merge this into the previous results
            for section, items in source_reqs.items():
                _merge_dependency_lists(items, reqs.setdefault(section, []))
        return reqs


def preprocess_for_bootstrap(paths, prebuilt_cctbx, platform=None, sections=None):
    # type: (list[str | os.PathLike], bool, str | None, list[SectionName]|None) -> list[str]
    """
    Do dependency file preprocessing, intended for bootstrap.py.

    Args:
        paths: List of dependency list files to merge
        prebuilt_cctbx: Whether this is processing for a prebuilt CCTBX
            distribution, or not.
        platform:
            The platform to process the dependencies for. Default: Current.
        sections:
            Which sections to process (build, host, run, test). Default: All.

    Returns:
        A list of dependency strings, suitable for passing to conda/mamba install.
    """
    parser = DependencySelectorParser(
        prebuilt_cctbx=prebuilt_cctbx,
        bootstrap=True,
        platform=platform or _native_platform(),
    )
    reqs = parser.parse_files(paths)
    merged_req = []
    for section, items in reqs.items():
        if section in sections or not sections:
            _merge_dependency_lists(items, merged_req)

    output_lines = []
    for pkg, ver, _ in sorted(merged_req, key=lambda x: x[0]):
        if pkg == "python":
            # Bootstrap handles this dependency implicitly
            continue
        output_lines.append("conda-forge::" + pkg + (ver or ""))
    return output_lines


def test_parser():
    parser = DependencySelectorParser(bootstrap=True, prebuilt_cctbx=False)
    assert parser._parse_expression("osx")
    assert parser._parse_expression("bootstrap")
    assert parser._parse_expression("osx and bootstrap")
    assert not parser._parse_expression("linux and bootstrap")
    assert not parser._parse_expression("prebuilt_cctbx and osx and not bootstrap")


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "--conda-build",
        help="Generate structured conda-build-style output",
        action="store_true",
    )
    parser.add_argument(
        "-p",
        "--platform",
        choices=["osx", "linux", "win"],
        help="Choose the target for handling bootstrap dependency lists. Default: {}".format(
            _native_platform()
        ),
    )
    parser.add_argument(
        "--prebuilt-cctbx",
        help="Mark as using prebuilt cctbx. Implied by conda-build.",
        action="store_true",
    )
    parser.add_argument(
        "--build",
        help="Include build section in output",
        dest="sections",
        action="append_const",
        const="build",
    )
    parser.add_argument(
        "--host",
        help="Include build section in output",
        dest="sections",
        action="append_const",
        const="host",
    )
    parser.add_argument(
        "--test",
        help="Include build section in output",
        dest="sections",
        action="append_const",
        const="test",
    )
    parser.add_argument(
        "--run",
        help="Include build section in output",
        dest="sections",
        action="append_const",
        const="run",
    )
    parser.add_argument(
        "-v", "--verbose", help="Show debugging output", action="store_true"
    )
    parser.add_argument("sources", nargs="+", help="Dependency files to merge")
    args = parser.parse_args()
    if not args.sections:
        args.sections = VALID_SECTIONS

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s"
    )
    if not args.conda_build:
        print(
            "\n".join(
                preprocess_for_bootstrap(
                    args.sources,
                    prebuilt_cctbx=args.prebuilt_cctbx,
                    platform=args.platform,
                    sections=args.sections,
                )
            )
        )
    else:
        if args.platform:
            sys.exit("Error: Can only specify platform with plain-list mode.")
        deps = DependencySelectorParser(bootstrap=False, prebuilt_cctbx=True)
        reqs = deps.parse_files(args.sources)
        output = []
        for section in VALID_SECTIONS:
            if section not in reqs or not reqs[section] or section not in args.sections:
                continue
            output.append(section + ":")
            output.extend(
                "    - " + x.raw_line
                for x in sorted(
                    reqs[section],
                    key=lambda x: (0 if x.raw_line.startswith("{{") else 1, x.raw_line),
                )
            )

        print("\n".join(output))
