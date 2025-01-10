#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 4 -*-

# Running bootstrap requires a minimum Python version of 2.7.

# To download this file:
# wget https://raw.githubusercontent.com/dials/dials/main/installer/bootstrap.py
# or
# curl https://raw.githubusercontent.com/dials/dials/main/installer/bootstrap.py > bootstrap.py

from __future__ import absolute_import, division, print_function

import argparse
import multiprocessing.pool
import os
import platform
import re
import shutil
import socket as pysocket
import stat
import subprocess
import sys
import tarfile
import threading
import time
import zipfile
import multiprocessing

try:  # Python 3
    from urllib.error import HTTPError, URLError
    from urllib.request import Request, urlopen
except ImportError:  # Python 2
    from urllib2 import HTTPError, Request, URLError, urlopen

try:
    # On windows, sometimes the python install is missing certificates
    import certifi  # noqa: F401
except ImportError:
    pass

try:
    from typing import Literal  # noqa: F401
except ImportError:
    pass

# Clean environment for subprocesses
clean_env = {
    key: value
    for key, value in os.environ.items()
    if key not in ("PYTHONPATH", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH")
}

devnull = open(os.devnull, "wb")  # to redirect unwanted subprocess output
allowed_ssh_connections = {}
concurrent_git_connection_limit = threading.Semaphore(5)


def make_executable(filepath):
    if os.name == "posix":
        mode = os.stat(filepath).st_mode
        mode |= (mode & 0o444) >> 2  # copy R bits to X
        # r--r--r-- => 0o444
        os.chmod(filepath, mode)


def _run_conda_retry(command_list):
    for retry in range(5):
        retry += 1
        try:
            run_command(
                command=command_list,
                workdir=".",
            )
        except Exception:
            print(
                """
*******************************************************************************
There was a failure in constructing the conda environment.
Attempt {retry} of 5 will start {retry} minute(s) from {t}.
*******************************************************************************
""".format(retry=retry, t=time.asctime())
            )
            time.sleep(retry * 60)
        else:
            break
    else:
        sys.exit(
            """
The conda environment could not be constructed. Please check that there is a
working network connection for downloading conda packages.
"""
        )


def get_requirements(conda_platform, conda_arch, python_version, is_cmake, extra_deps):
    # type: (str, Literal["linux", "macos", "windows"], str, bool, list[str] | None) -> str

    """
    Find or create a file of platform-specific dependencies
    """
    # Searching for definitions first looks for platform/version locked:
    # modules/dials/.conda-envs/{conda_arch}-py{python_version}.txt: Full release locked files
    #
    # Otherwise, look for and merge the new requirements format. It is
    # an error if these are not present (e.g. before the update step):
    #
    # - modules/dials/dependencies.yaml
    # - modules/dxtbx/dependencies.yaml
    # - modules/xia2/dependencies.yaml
    #
    # And, if not using a prebuild cctbx:
    # - modules/dials/.conda_envs/cctbx-dependencies.yaml

    # Identify packages required for environment
    env_dir = os.path.join("modules", "dials", ".conda-envs")
    # First, check to see if we have an architecture-and-python-specific environment file
    filename = os.path.join(env_dir, "{0}_py{1}.txt".format(conda_arch, python_version))
    if os.path.isfile(filename):
        return filename

    expected_dependency_lists = [
        "modules/dials/dependencies.yaml",
        "modules/dxtbx/dependencies.yaml",
        # "modules/xia2/dependencies.yaml",
    ]
    if not is_cmake:
        expected_dependency_lists.append(
            "modules/dials/.conda-envs/cctbx-dependencies.yaml"
        )
    expected_dependency_lists.extend(extra_deps or [])

    for reqfile in expected_dependency_lists:
        if not os.path.isfile(reqfile):
            # We need this...
            sys.exit(
                "Error: Could not find dependency list {}. Have you run 'bootstrap update'?".format(
                    reqfile
                )
            )

    # Run the dependency merger
    prebuilt = ["--prebuilt-cctbx"] if is_cmake else []
    platform_selectors = {"linux": "linux", "macos": "osx", "windows": "win"}
    results = subprocess.check_output(
        [
            sys.executable,
            "modules/dials/util/parse_dependency_selectors.py",
            "-p",
            platform_selectors[conda_platform],
        ]
        + prebuilt
        + expected_dependency_lists
    )
    filename = "modules/dials/.conda-envs/requirements.txt"
    with open(filename, "wb") as f:
        f.write(results)
    return filename


def install_micromamba(python, cmake, extra_deps):
    # type: (str, bool, list[str] | None) -> None
    """Download and install Micromamba"""
    if sys.platform.startswith("linux"):
        conda_platform = "linux"
        conda_arch = "linux-64"
        member = "bin/micromamba"
    elif sys.platform == "darwin":
        conda_platform = "macos"
        member = "bin/micromamba"
        if platform.machine() == "arm64":
            conda_arch = "osx-arm64"
        else:
            conda_arch = "osx-64"
    elif os.name == "nt":
        conda_platform = "windows"
        member = "Library/bin/micromamba.exe"
        conda_arch = "win-64"
    else:
        raise NotImplementedError(
            "Unsupported platform %s / %s" % (os.name, sys.platform)
        )
    url = "https://micro.mamba.pm/api/micromamba/{0}/latest".format(conda_arch)
    mamba_prefix = os.path.realpath("micromamba")
    clean_env["MAMBA_ROOT_PREFIX"] = mamba_prefix
    mamba = os.path.join(mamba_prefix, member.split("/")[-1])
    print("Downloading {url}:".format(url=url), end=" ")
    result = download_to_file(url, os.path.join(mamba_prefix, "micromamba.tar.bz2"))
    if result in (0, -1):
        sys.exit("Micromamba download failed")
    with tarfile.open(
        os.path.join(mamba_prefix, "micromamba.tar.bz2"), "r:bz2"
    ) as tar, open(mamba, "wb") as fh:
        fh.write(tar.extractfile(member).read())
    make_executable(mamba)

    # verify micromamba works and check version
    conda_info = subprocess.check_output([mamba, "--version"], env=clean_env)
    if sys.version_info.major > 2:
        conda_info = conda_info.decode("latin-1")
    print("Using Micromamba version", conda_info.strip())

    # Find or generate the requirements list
    filename = get_requirements(
        conda_platform=conda_platform,
        conda_arch=conda_arch,
        python_version=python,
        is_cmake=cmake,
        extra_deps=extra_deps,
    )
    # install a new environment or update an existing one
    prefix = os.path.realpath("conda_base")
    if os.path.exists(prefix):
        command = "install"
        text_messages = ["Updating", "update of"]
    else:
        command = "create"
        text_messages = ["Installing", "installation into"]
    python_requirement = "conda-forge::python=%s.*" % python

    command_list = [
        mamba,
        "--no-env",
        "--no-rc",
        "--prefix",
        prefix,
        "--root-prefix",
        mamba_prefix,
        command,
        "--file",
        filename,
        "--yes",
        "--channel",
        "conda-forge",
        "--override-channels",
        "mamba",
        python_requirement,
    ]
    extra_deps = []

    if cmake:
        # This should no longer be necessary... except for now we also want to
        # continue to support the DIALS-release-style explicit lock file
        # with CMake options. So look for an explicit lockfile, then if
        # present just have a hardcoded list of extra packages.
        #
        # This is a bad solution but will keep us going until we move
        # to a better way of doing releases.
        with open(filename) as dep_file:
            is_explicit_reqs = "@EXPLICIT" in dep_file.read()

        if is_explicit_reqs:
            extra_deps = [
                "cctbx-base",
                "pycbf",
                "cmake",
                "pre-commit",
                "libboost-python-devel",  # Have seen this removed by explicit update
            ]

    print(
        "{text} dials environment from {filename} with Python {python}".format(
            text=text_messages[0], filename=filename, python=python
        )
    )

    _run_conda_retry(command_list)
    # If we wanted extra dependencies, and couldn't install them directly, run again
    if extra_deps:
        command_list = [
            mamba,
            "--no-env",
            "--no-rc",
            "--prefix",
            prefix,
            "--root-prefix",
            mamba_prefix,
            "install",
            "--yes",
            "--channel",
            "conda-forge",
            "--override-channels",
            "mamba",
        ] + extra_deps
        _run_conda_retry(command_list)

    print("Completed {text}:\n  {prefix}".format(text=text_messages[1], prefix=prefix))
    with open(os.path.join(prefix, ".condarc"), "w") as fh:
        fh.write(
            """
changeps1: False
channels:
  - conda-forge
""".lstrip()
        )


def run_command(command, workdir):
    print("Running %s (in %s)" % (" ".join(command), workdir))
    workdir = os.path.abspath(workdir)
    try:
        os.makedirs(workdir)
    except OSError:
        pass
    try:
        p = subprocess.Popen(args=command, cwd=workdir, env=clean_env)
    except Exception as e:
        if isinstance(e, OSError):
            if e.errno == 2:
                executable = os.path.normpath(os.path.join(workdir, command[0]))
                raise RuntimeError("Could not run %s: File not found" % executable)
        if "child_traceback" in dir(e):
            print("Calling subprocess resulted in error; ", e.child_traceback)
        raise e

    try:
        p.wait()
    except KeyboardInterrupt:
        print("\nReceived CTRL+C, trying to stop subprocess...\n")
        p.terminate()
        raise
    if p.returncode:
        sys.exit("Process failed with return code %s" % p.returncode)


def run_indirect_command(command, args):
    print("(via conda environment) " + command)
    try:
        os.mkdir("build")
    except OSError:
        pass
    if os.name == "nt":
        filename = os.path.join("build", "indirection.cmd")
        with open(filename, "w") as fh:
            fh.write("call %s\\conda_base\\condabin\\activate.bat\r\n" % os.getcwd())
            fh.write("shift\r\n")
            fh.write("%*\r\n")
        if not command.endswith((".bat", ".cmd", ".exe")):
            command = command + ".bat"
        indirection = ["cmd.exe", "/C", "indirection.cmd"]
    else:
        filename = os.path.join("build", "indirection.sh")
        with open(filename, "w") as fh:
            fh.write("#!/bin/bash\n")
            fh.write("source %s/conda_base/etc/profile.d/conda.sh\n" % os.getcwd())
            fh.write("conda activate %s/conda_base\n" % os.getcwd())
            fh.write('"$@"\n')
        make_executable(filename)
        indirection = ["./indirection.sh"]
    run_command(
        command=indirection + [command] + args,
        workdir="build",
    )


def download_to_file(url, file, quiet=False):
    """Downloads a URL to file. Returns the file size.
    Returns -1 if the downloaded file size does not match the expected file
    size
    Returns -2 if the download is skipped due to the file at the URL not
    being newer than the local copy (identified by matching timestamp and
    size)
    """

    # Create directory structure if necessary
    if os.path.dirname(file):
        try:
            os.makedirs(os.path.dirname(file))
        except Exception:
            pass

    localcopy = os.path.isfile(file)

    try:
        from ssl import SSLError
    except ImportError:
        SSLError = None

    # Open connection to remote server
    try:
        url_request = Request(url)
        if localcopy:
            # Shorten timeout to 7 seconds if a copy of the file is already present
            socket = urlopen(url_request, None, 7)
        else:
            socket = urlopen(url_request)
    except SSLError as e:
        # This could be a timeout
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            if not quiet:
                print(str(e))
            return -2
        # otherwise pass on the error message
        raise
    except (pysocket.timeout, HTTPError) as e:
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            if not quiet:
                print(str(e))
            return -2
        # otherwise pass on the error message
        raise
    except URLError as e:
        if localcopy:
            # Download failed for some reason, but a valid local copy of
            # the file exists, so use that one instead.
            if not quiet:
                print(str(e))
            return -2
        # if url fails to open, try using curl
        # temporary fix for old OpenSSL in system Python on macOS
        # https://github.com/cctbx/cctbx_project/issues/33
        command = ["/usr/bin/curl", "-fLo", file, "--retry", "5", url]
        subprocess.call(command)
        socket = None  # prevent later socket code from being run
        try:
            received = os.path.getsize(file)
        except OSError:
            raise RuntimeError("Download failed")

    if socket is not None:
        try:
            file_size = int(socket.info().get("Content-Length"))
        except Exception:
            file_size = 0

        remote_mtime = 0
        try:
            remote_mtime = time.mktime(socket.info().getdate("last-modified"))
        except Exception:
            pass

        if file_size > 0:
            if remote_mtime > 0:
                # check if existing file matches remote size and timestamp
                try:
                    (
                        mode,
                        ino,
                        dev,
                        nlink,
                        uid,
                        gid,
                        size,
                        atime,
                        mtime,
                        ctime,
                    ) = os.stat(file)
                    if (size == file_size) and (remote_mtime == mtime):
                        if not quiet:
                            print("local copy is current")
                        socket.close()
                        return -2
                except Exception:
                    # proceed with download if timestamp/size check fails for any reason
                    pass

            hr_size = (file_size, "B")
            if hr_size[0] > 500:
                hr_size = (hr_size[0] / 1024, "kB")
            if hr_size[0] > 500:
                hr_size = (hr_size[0] / 1024, "MB")
            if not quiet:
                print("%.1f %s" % hr_size)
                print("    [0%", end="")
                sys.stdout.flush()  # becomes print(flush=True) when we move to 3.3+

        received = 0
        block_size = 8192
        progress = 1
        # Write to the file immediately so we can empty the buffer
        tmpfile = file + ".tmp"

        with open(tmpfile, "wb") as fh:
            while True:
                block = socket.read(block_size)
                received += len(block)
                fh.write(block)
                if file_size > 0 and not quiet:
                    while (100 * received / file_size) > progress:
                        progress += 1
                        if (progress % 20) == 0:
                            print(progress, end="%")
                            sys.stdout.flush()  # becomes print(flush=True) when we move to 3.3+
                        elif (progress % 2) == 0:
                            print(".", end="")
                            sys.stdout.flush()  # becomes print(flush=True) when we move to 3.3+
                if not block:
                    break
        socket.close()

        if not quiet:
            if file_size > 0:
                print("]")
            else:
                print("%d kB" % (received / 1024))
            sys.stdout.flush()  # becomes print(flush=True) when we move to 3.3+

        # Do not overwrite file during the download. If a download temporarily fails we
        # may still have a clean, working (yet older) copy of the file.
        shutil.move(tmpfile, file)

        if (file_size > 0) and (file_size != received):
            return -1

        if remote_mtime > 0:
            # set file timestamp if timestamp information is available
            st = os.stat(file)
            atime = st[stat.ST_ATIME]  # current access time
            os.utime(file, (atime, remote_mtime))

    return received


def unzip(archive, directory, trim_directory=0):
    """unzip a file into a directory."""
    if not zipfile.is_zipfile(archive):
        raise Exception(
            "Cannot install %s: %s is not a valid .zip file" % (directory, archive)
        )
    z = zipfile.ZipFile(archive, "r")
    for member in z.infolist():
        is_directory = member.filename.endswith("/")
        filename = os.path.join(*member.filename.split("/")[trim_directory:])
        if filename != "":
            filename = os.path.normpath(filename)
            if "../" in filename:
                raise Exception(
                    "Archive %s contains invalid filename %s" % (archive, filename)
                )
            filename = os.path.join(directory, filename)
            upperdirs = os.path.dirname(filename)
            try:
                if is_directory and not os.path.exists(filename):
                    os.makedirs(filename)
                elif upperdirs and not os.path.exists(upperdirs):
                    os.makedirs(upperdirs)
            except Exception:
                pass
            if not is_directory:
                with z.open(member) as source:
                    with open(filename, "wb") as target:
                        shutil.copyfileobj(source, target)

                # Preserve executable permission, if set
                unix_executable = member.external_attr >> 16 & 0o111
                # rwxrwxrwx => --x--x--x => 0o111
                if unix_executable:
                    make_executable(filename)
    z.close()


def set_git_repository_config_to_rebase(config):
    with open(config, "r") as fh:
        cfg = fh.readlines()

    branch, remote, rebase = False, False, False
    insertions = []
    for n, line in enumerate(cfg):
        if line.startswith("["):
            if branch and remote and not rebase:
                insertions.insert(0, (n, branch))
            if line.startswith("[branch"):
                branch = line.split('"')[1]
            else:
                branch = False
            remote, rebase = False, False
        if re.match(r"remote\s*=", line.strip()):
            remote = True
        if re.match(r"rebase\s*=", line.strip()):
            rebase = True
    if branch and remote and not rebase:
        insertions.insert(0, (n + 1, branch))
    for n, branch in insertions:
        cfg.insert(n, "\trebase = true\n")
    with open(config, "w") as fh:
        fh.write("".join(cfg))


def git(module, git_available, ssh_available, reference_base, settings):
    """Retrieve a git repository, either by running git directly
    or by downloading and unpacking an archive.
    """
    destination = os.path.join("modules", module)

    if os.path.exists(destination):
        if os.path.isfile(os.path.join(destination, ".git")):
            return module, "WARNING", "Existing git worktree directory -- skipping"
        if not os.path.isdir(os.path.join(destination, ".git")):
            return module, "WARNING", "Existing non-git directory -- skipping"
        if not git_available:
            return module, "WARNING", "Cannot update module, git command not found"

        with open(os.path.join(destination, ".git", "HEAD"), "r") as fh:
            if fh.read(4) != "ref:":
                return (
                    module,
                    "WARNING",
                    "Cannot update existing git repository! You are not on a branch.\n"
                    "This may be legitimate when run eg. via Jenkins, but be aware that you cannot commit any changes",
                )

        with concurrent_git_connection_limit:
            p = subprocess.Popen(
                args=["git", "pull", "--rebase"],
                cwd=destination,
                env=clean_env,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            # This may fail for unclean trees and merge problems. In this case manual
            # user intervention will be required.
            # For the record, you can clean up the tree and *discard ALL changes* with
            #   git reset --hard origin/master
            #   git clean -dffx
            try:
                output, _ = p.communicate()
                output = output.decode("latin-1")
            except KeyboardInterrupt:
                print("\nReceived CTRL+C, trying to terminate subprocess...\n")
                p.terminate()
                raise
        if p.returncode:
            return (
                module,
                "WARNING",
                "Cannot update existing git repository! Unclean tree or merge problems.\n"
                + output,
            )
        # Show the hash for the checked out commit for debugging purposes
        p = subprocess.Popen(
            args=["git", "rev-parse", "HEAD", "--abbrev-ref", "HEAD"],
            cwd=destination,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        output, _ = p.communicate()
        output = output.decode("latin-1")
        if p.returncode:
            return module, "WARNING", "Cannot get git repository revision\n" + output
        output = output.split()
        if len(output) == 2:
            return module, "OK", "Checked out revision %s (%s)" % (output[0], output[1])
        return module, "OK", "Checked out revision " + output[0].strip()

    try:
        os.makedirs("modules")
    except OSError:
        pass

    remote_branch = settings.get("branch-remote", settings["branch-local"])

    if not git_available:
        # Fall back to downloading a static archive
        url = "https://github.com/%s/archive/%s.zip" % (
            settings.get("effective-repository", settings.get("base-repository")),
            remote_branch,
        )
        filename = os.path.join("modules", "%s-%s.zip" % (module, remote_branch))
        try:
            download_to_file(url, filename, quiet=True)
        except Exception:
            print("Error downloading", url)
            raise
        unzip(filename, destination, trim_directory=1)
        return module, "OK", "Downloaded branch %s from static archive" % remote_branch

    if ssh_available:
        remote_pattern = "git@github.com:%s.git"
    else:
        remote_pattern = "https://github.com/%s.git"

    if git_available:
        # Determine the git version
        p = subprocess.Popen(
            ["git", "--version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        output, _ = p.communicate()
        output = output.decode("latin-1")
        parts = output.split(" ", 2)
        if p.returncode or not parts[:2] == ["git", "version"]:
            raise RuntimeError("Could not determine git version")
        # Version comes in:
        #    "git version x.y.z"
        # or "git version x.y.z.windows.n"
        git_version = tuple(int(x) if x.isnumeric() else x for x in parts[2].split("."))

    secondary_remote = settings.get("effective-repository") and (
        settings["effective-repository"] != settings.get("base-repository")
    )
    direct_branch_checkout = []
    if not secondary_remote and remote_branch == settings["branch-local"]:
        direct_branch_checkout = ["-b", remote_branch]
    reference_parameters = []
    if reference_base and os.path.exists(os.path.join(reference_base, module, ".git")):
        # Use -if-able so that we don't have errors over unreferenced submodules
        reference_type = "--reference-if-able"
        if git_version < (2, 11, 0):
            # As a fallback, use the old parameter. This will fail if
            # there are submodules and the reference does not have all
            # the required submodules
            reference_type = "--reference"
        reference_parameters = [
            reference_type,
            os.path.join(reference_base, module),
        ]

    with concurrent_git_connection_limit:
        p = subprocess.Popen(
            args=["git", "clone", "--recursive"]
            + direct_branch_checkout
            + [
                remote_pattern
                % settings.get("base-repository", settings.get("effective-repository")),
                module,
            ]
            + reference_parameters,
            cwd="modules",
            env=clean_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        try:
            output, _ = p.communicate()
            output = output.decode("latin-1")
        except KeyboardInterrupt:
            print("\nReceived CTRL+C, trying to terminate subprocess...\n")
            p.terminate()
            raise
        if p.returncode:
            return (module, "ERROR", "Cannot checkout git repository\n" + output)

    if reference_parameters:
        # Sever the link between checked out and reference repository
        returncode = subprocess.call(
            ["git", "repack", "-a", "-d"],
            cwd=destination,
            env=clean_env,
            stdout=devnull,
            stderr=devnull,
        )
        if returncode:
            return (
                module,
                "ERROR",
                "Could not detach git repository from reference. Repository may be in invalid state!\n"
                "Run 'git repack -a -d' in the repository, or delete and recreate it",
            )
        alternates = os.path.join(destination, ".git", "objects", "info", "alternates")
        if os.path.exists(alternates):
            os.remove(alternates)

    if secondary_remote:
        returncode = subprocess.call(
            [
                "git",
                "remote",
                "add",
                "upstream",
                remote_pattern % settings["effective-repository"],
            ],
            cwd=destination,
            env=clean_env,
            stdout=devnull,
            stderr=devnull,
        )
        if returncode:
            return (
                module,
                "ERROR",
                "Could not add upstream remote to repository. Repository may be in invalid state!",
            )
        with concurrent_git_connection_limit:
            returncode = subprocess.call(
                ["git", "fetch", "upstream"],
                cwd=destination,
                env=clean_env,
                stdout=devnull,
                stderr=devnull,
            )
        if returncode:
            return (
                module,
                "ERROR",
                "Could not fetch upstream repository %s. Repository may be in invalid state!"
                % settings["effective-repository"],
            )

    set_git_repository_config_to_rebase(os.path.join(destination, ".git", "config"))

    if not direct_branch_checkout:
        # set up the local branch with tracking
        returncode = subprocess.call(
            [
                "git",
                "checkout",
                "-B",
                settings["branch-local"],
                "--track",
                "%s/%s" % ("upstream" if secondary_remote else "origin", remote_branch),
            ],
            cwd=destination,
            env=clean_env,
            stdout=devnull,
            stderr=devnull,
        )
        if returncode:
            return (
                module,
                "ERROR",
                "Could not check out alternate branch %s. Repository may be in invalid state!"
                % remote_branch,
            )

        if secondary_remote:
            # When checking out a branch from a secondary repository under a
            # different name set up a git pre-commit hook to write protect
            # this branch.
            hookfile = os.path.join(destination, ".git", "hooks", "pre-commit")
            with open(hookfile, "w") as fh:
                fh.write(
                    """
#!/bin/sh

#
# Reject commits to '{branch}' to avoid accidental
# commits to the alternate repository {repository}
#

if [ "$(git rev-parse --abbrev-ref HEAD)" == "{branch}" ]
then
    echo "Please do not commit to the '{branch}' branch."
    echo "You can create a new branch or commit directly to master."
    echo
    exit 1
fi
""".format(
                        branch=settings["branch-local"],
                        repository=settings["effective-repository"],
                    )
                )
            make_executable(hookfile)

    # Show the hash for the checked out commit for debugging purposes
    p = subprocess.Popen(
        args=["git", "rev-parse", "HEAD"],
        cwd=destination,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    output, _ = p.communicate()
    output = output.decode("latin-1")
    if p.returncode:
        return (
            module,
            "WARNING",
            "Cannot get git repository revision\n" + output,
        )
    git_status = settings["branch-local"]
    if settings["branch-local"] != remote_branch:
        git_status += " tracking " + remote_branch
    if secondary_remote:
        git_status += " at " + settings["effective-repository"]
    return module, "OK", "Checked out revision %s (%s)" % (output.strip(), git_status)


def update_sources(options):
    if options.git_reference:
        reference_base = os.path.abspath(os.path.expanduser(options.git_reference))
    else:
        if os.name == "posix" and pysocket.gethostname().endswith(".diamond.ac.uk"):
            reference_base = (
                "/dls/science/groups/scisoft/DIALS/repositories/git-reference"
            )
        else:
            reference_base = None
    try:
        git_available = not subprocess.call(
            ["git", "--version"], stdout=devnull, stderr=devnull
        )
    except OSError:
        git_available = False
    ssh_available = False
    if git_available:
        try:
            returncode = subprocess.call(
                [
                    "ssh",
                    "-oBatchMode=yes",
                    "-oStrictHostKeyChecking=no",
                    "-T",
                    "git@github.com",
                ],
                stdout=devnull,
                stderr=devnull,
            )
            # SSH errors lead to 255
            ssh_available = returncode in (0, 1)
        except OSError:
            pass

    if not options.cmake:
        repositories = {
            source.split("/")[1]: {"base-repository": source, "branch-local": branch}
            for source, branch in (
                ("cctbx/annlib_adaptbx", "master"),
                ("cctbx/cctbx_project", "master"),
                ("cctbx/dxtbx", "main"),
                ("dials/annlib", "master"),
                ("dials/cbflib", "main"),
                ("dials/ccp4io", "master"),
                ("dials/ccp4io_adaptbx", "master"),
                ("dials/dials", "main"),
                ("dials/gui_resources", "master"),
                ("xia2/xia2", "main"),
            )
        }
        repositories["cctbx_project"] = {
            "base-repository": "cctbx/cctbx_project",
            "effective-repository": "dials/cctbx",
            "branch-remote": "master",
            "branch-local": "stable",
        }
    else:
        # Only what we need for CMake
        repositories = {
            source.split("/")[1]: {"base-repository": source, "branch-local": branch}
            for source, branch in (
                ("cctbx/dxtbx", "main"),
                ("dials/dials", "main"),
                ("xia2/xia2", "main"),
            )
        }

    for source, setting in options.branch:
        if source not in repositories:
            sys.exit("Unknown repository %s" % source)
        setting = re.match(
            r"^(?:(\w+/\w+)@?)?([a-zA-Z0-9._\-]+)?(?::([a-zA-Z0-9._\-]+))?$", setting
        )
        if not setting:
            sys.exit("Could not parse the branch setting for repository %s" % source)
        _repository, _branch_remote, _branch_local = setting.groups()
        if _repository:
            repositories[source] = {
                "base-repository": _repository,
                "branch-remote": _branch_remote or "master",
                "branch-local": _branch_local or _branch_remote or "master",
            }
        elif _branch_remote:
            repositories[source]["branch-remote"] = _branch_remote
            repositories[source]["branch-local"] = _branch_local or _branch_remote
        elif _branch_local:
            repositories[source]["branch-local"] = _branch_local

    def _git_fn(repository):
        return git(
            repository,
            git_available,
            ssh_available,
            reference_base,
            repositories[repository],
        )

    success = True
    update_pool = multiprocessing.pool.ThreadPool(20)
    try:
        for result in update_pool.imap_unordered(_git_fn, repositories):
            module, result, output = result
            output = (result + " - " + output).replace(
                "\n", "\n" + " " * (len(module + result) + 5)
            )
            if result == "ERROR":
                success = False
            if os.name == "posix" and sys.stdout.isatty():
                if result == "OK":
                    output = "\x1b[32m" + output + "\x1b[0m"
                elif result == "WARNING":
                    output = "\x1b[33m" + output + "\x1b[0m"
                elif result == "ERROR":
                    output = "\x1b[31m" + output + "\x1b[0m"
            print(module + ": " + output)
    except KeyboardInterrupt:
        update_pool.terminate()
        sys.exit("\naborted with Ctrl+C")
    except Exception:
        update_pool.terminate()
        raise
    update_pool.close()
    update_pool.join()
    if not success:
        sys.exit("\nFailed to update one or more repositories")


def run_tests():
    dispatch_extension = ".bat" if os.name == "nt" else ""
    print("Running dxtbx tests")
    run_command(
        [
            os.path.join(
                "..", "..", "build", "bin", "libtbx.pytest" + dispatch_extension
            ),
            "--regression",
            "-n",
            "auto",
        ],
        workdir=os.path.join("modules", "dxtbx"),
    )
    print("Running dials tests")
    run_command(
        [
            os.path.join(
                "..", "..", "build", "bin", "libtbx.pytest" + dispatch_extension
            ),
            "--regression",
            "-n",
            "auto",
        ],
        workdir=os.path.join("modules", "dials"),
    )


def refresh_build():
    print("Running libtbx.refresh")
    run_indirect_command(
        os.path.join("bin", "libtbx.refresh"),
        args=[],
    )


def install_precommit(cmake=False):
    print("Installing precommits")
    if cmake:
        # Just find all repository directories under modules, since we
        # don't have a central registry
        subdirs = []
        for sub in os.listdir("modules"):
            if os.path.isdir(os.path.join("modules", sub, ".git")):
                subdirs.append(os.path.abspath(os.path.join("modules", sub)))

        conda_base_root = os.path.join(os.path.abspath("."), "conda_base")
        if os.name == "nt":
            precommit_command = os.path.join(
                conda_base_root, "Scripts", "libtbx.precommit.exe"
            )
        else:
            precommit_command = os.path.join(conda_base_root, "bin", "libtbx.precommit")
        run_indirect_command(
            precommit_command,
            args=["install"] + subdirs,
        )
    else:
        run_indirect_command(
            os.path.join("bin", "libtbx.precommit"),
            args=["install"],
        )


def _get_base_python():
    if os.name == "nt":
        conda_python = os.path.join(os.getcwd(), "conda_base", "python.exe")
    elif sys.platform.startswith("darwin"):
        conda_python = os.path.join(
            "conda_base", "python.app", "Contents", "MacOS", "python"
        )
    else:
        conda_python = os.path.join("conda_base", "bin", "python")
    return os.path.abspath(conda_python)


def _get_cmake_exe():
    if os.name == "nt":
        return os.path.abspath(
            os.path.join("conda_base", "Library", "bin", "cmake.exe")
        )
    else:
        return os.path.abspath(os.path.join("conda_base", "bin", "cmake"))


def refresh_build_cmake():
    conda_python = _get_base_python()
    run_indirect_command(
        conda_python,
        [
            "-mpip",
            "install",
            "--no-deps",
            "-e",
            "../modules/dxtbx",
            "-e",
            "../modules/dials",
            "-e",
            "../modules/xia2",
        ],
    )


def configure_build_cmake(extra_args):
    # type: (list[str] | None) -> None
    cmake_exe = _get_cmake_exe()
    python_exe = _get_base_python()

    # Get the location of site-packages
    site_path = subprocess.check_output(
        [
            python_exe,
            "-c",
            "import os, site; print(os.path.abspath(site.getsitepackages()[0]))",
        ]
    ).strip()
    if not isinstance(site_path, str):
        site_path = site_path.decode()

    # Write a .pth file here pointing to the build/lib folder. This
    # is a somewhat "clever" convenience such that the PYTHONPATH doesn't
    # need to be explicitly set before running code.
    lib_pth = os.path.join(site_path, "__bootstrap__.dials.pth")
    print("Writing libdir .pth file to: " + lib_pth)
    # Best-guess work out where the libraries are put by the build
    if os.name == "nt":
        build_lib_dir = os.path.join(os.getcwd(), "build", "lib", "RelWithDebInfo")
    else:
        build_lib_dir = os.path.join(os.getcwd(), "build", "lib")

    with open(lib_pth, "w") as f:
        f.write(build_lib_dir)

    # write a new-style environment setup script
    if os.name == "nt":
        activate = os.path.join(os.getcwd(), "conda_base", "condabin", "activate.bat")
        with open("dials.bat", "w") as f:
            f.write(
                """\
rem enable conda environment
call {}
conda activate {}
""".format(
                    activate,
                    os.path.join(os.getcwd(), "conda_base"),
                )
            )
    else:
        with open("dials", "w") as f:
            f.write(
                """\
# enable conda environment
source {dist_root}/conda_base/etc/profile.d/conda.sh
conda activate {dist_root}/conda_base
""".format(
                    dist_root=os.getcwd(),
                )
            )

    # Write a compound CMakeLists.txt, if one doesn't exist
    if not os.path.isfile("modules/CMakeLists.txt"):
        with open("modules/CMakeLists.txt", "w") as f:
            f.write(
                """\
cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(dials)

if (CMAKE_UNITY_BUILD AND MSVC)
    # Windows can fail in this scenario because too many objects
    add_compile_options(/bigobj)
endif()

add_subdirectory(dxtbx)
add_subdirectory(dials)
"""
            )

    # run_indirect runs inside the build folder with an activated environment
    conda_base_root = os.path.join(os.path.abspath("."), "conda_base")
    assert os.path.isdir(conda_base_root)
    extra_args = extra_args or []
    if os.name == "nt":
        extra_args.append("-DPython_ROOT_DIR=" + conda_base_root)
    sys.stdout.flush()
    run_indirect_command(
        cmake_exe,
        [
            "../modules",
            "-DCMAKE_INSTALL_PREFIX=" + conda_base_root,
            "-DHDF5_DIR=" + conda_base_root,
            "-DPython_ROOT_DIR=" + conda_base_root,
        ]
        + extra_args,
    )


def configure_build(config_flags):
    conda_python = _get_base_python()

    # write a new-style environment setup script
    with open(("dials.bat" if os.name == "nt" else "dials"), "w"):
        pass  # ensure we write a new-style environment setup script

    if os.name != "nt" and not any(
        flag.startswith("--compiler=") for flag in config_flags
    ):
        config_flags.append("--compiler=conda")
    if "--use_conda" not in config_flags:
        config_flags.append("--use_conda")
    # Default to C++14 if otherwise unspecified
    if not any(x.startswith("--cxxstd") for x in config_flags):
        config_flags.append("--cxxstd=c++14")

    print("Setting up build directory")

    configcmd = [
        os.path.join("..", "modules", "cctbx_project", "libtbx", "configure.py"),
        "cctbx",
        "scitbx",
        "libtbx",
        "iotbx",
        "mmtbx",
        "smtbx",
        "gltbx",
        "wxtbx",
        "cbflib",
        "dxtbx",
        "xfel",
        "dials",
        "xia2",
        "prime",
        "--skip_phenix_dispatchers",
        "--use_environment",
    ] + config_flags

    run_indirect_command(
        command=conda_python,
        args=configcmd,
    )


def make_build():
    try:
        nproc = len(os.sched_getaffinity(0))
    except AttributeError:
        nproc = multiprocessing.cpu_count()
    command = os.path.join("bin", "libtbx.scons")

    run_indirect_command(command, args=["-j", str(nproc)])
    # run build again to make sure everything is built
    run_indirect_command(command, args=["-j", str(nproc)])


def make_build_cmake():
    cmake_exe = _get_cmake_exe()
    if os.name == "nt":
        run_indirect_command(cmake_exe, ["--build", ".", "--config", "RelWithDebInfo"])
    else:
        parallel = []
        if "CMAKE_GENERATOR" not in os.environ:
            if hasattr(os, "sched_getaffinity"):
                cpu = os.sched_getaffinity(0)
            else:
                cpu = multiprocessing.cpu_count()
            if isinstance(cpu, int):
                parallel = ["--parallel", str(cpu)]
        run_indirect_command(cmake_exe, ["--build", "."] + parallel)


def repository_at_tag(string):
    try:
        repository, tag = string.split("@", 1)
        return (repository, tag)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "%s does not follow the repository@branch format" % string
        )


class Choices(tuple):
    # Python bug https://bugs.python.org/issue27227, https://bugs.python.org/issue9625
    def __new__(cls, *args, **kwargs):
        x = tuple.__new__(cls, *args, **kwargs)
        Choices.__init__(x, *args, **kwargs)
        return x

    def __init__(self, *args, **kwargs):
        self.default = []

    def __contains__(self, item):
        return tuple.__contains__(self, item) or item is self.default


def run():
    prog = os.environ.get("LIBTBX_DISPATCHER_NAME")
    if prog is None or prog.startswith("python") or prog.endswith("python"):
        prog = os.path.basename(sys.argv[0])

    description = """
  You may specify one or more actions:
    update - Update source repositories (cctbx, cbflib, etc.)
    base - Create conda environment
    build - Build
    test - Run tests

  The default actions are: update, base, build

  Example:

    python bootstrap.py update base build test
  """

    parser = argparse.ArgumentParser(
        prog=prog,
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    action_choices = Choices(("update", "base", "build", "test"))
    action_choices.default = ["update", "base", "build"]
    parser.add_argument(
        "actions",
        nargs="*",
        help="actions for building as described above",
        choices=action_choices,
        default=action_choices.default,
    )
    parser.add_argument(
        "--git-reference",
        help="Path to a directory containing reference copies of repositories for faster checkouts.",
    )
    parser.add_argument(
        "--config-flags",
        help="""Pass flags to the configuration step (CMake or libtbx). Flags should
be passed separately with quotes to avoid confusion (e.g
--config-flags="--build=debug" --config-flags="--another_flag")""",
        action="append",
        default=[],
    )
    parser.add_argument(
        "--python",
        help="Install this minor version of Python (default: %(default)s)",
        default="3.12",
        choices=("3.10", "3.11", "3.12"),
    )
    parser.add_argument(
        "--branch",
        type=repository_at_tag,
        action="append",
        default=[],
        help=(
            "during 'update' step when a repository is newly cloned set it to a given branch."
            "Specify as repository@branch, eg. 'dials@dials-next'"
        ),
    )
    parser.add_argument(
        "--clean",
        help="Remove temporary conda environments and package caches after installation",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "--libtbx",
        help="Use the libtbx build system, compiling cctbx from scratch.",
        action="store_false",
        dest="cmake",
    )
    parser.add_argument(
        "--cmake",
        action="store_true",
        dest="removed_cmake",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--extra-dependencies",
        action="append",
        help=argparse.SUPPRESS,
    )

    options = parser.parse_args()
    if options.removed_cmake:
        # User passed the obsolete parameter
        sys.exit("Error: --cmake is now the default, please remove --cmake.")

    print("Performing actions:", " ".join(options.actions))

    # Add sources
    if "update" in options.actions:
        update_sources(options)

    # Build base packages
    if "base" in options.actions:
        install_micromamba(
            options.python,
            cmake=options.cmake,
            extra_deps=options.extra_dependencies,
        )
        if options.clean:
            shutil.rmtree(os.path.realpath("micromamba"))

    # Configure, make, get revision numbers
    if "build" in options.actions:
        if options.cmake:
            refresh_build_cmake()
            configure_build_cmake(options.config_flags)
            make_build_cmake()
        else:
            configure_build(options.config_flags)
            make_build()
            refresh_build()
        install_precommit(options.cmake)

    # Tests, tests
    if "test" in options.actions:
        run_tests()

    print("\nBootstrap success: %s" % ", ".join(options.actions))


if __name__ == "__main__":
    run()
