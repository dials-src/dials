from __future__ import absolute_import, division, print_function

import collections
import os
import re

from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from libtbx.utils import Sorry
import libtbx.phil


FilenameDataWrapper = collections.namedtuple("FilenameDataWrapper", "filename, data")


class ExperimentListConverters(object):
    """A phil converter for the experiment list class."""

    phil_type = "experiment_list"

    def __init__(self, check_format=True):
        self._check_format = check_format

    def __str__(self):
        return self.phil_type

    def from_words(self, words, master):
        s = libtbx.phil.str_from_words(words=words)
        if s is None:
            return None
        if not os.path.exists(s):
            raise Sorry("File %s does not exist" % s)
        return FilenameDataWrapper(
            s, ExperimentListFactory.from_json_file(s, check_format=self._check_format),
        )

    def as_words(self, python_object, master):
        if python_object is None:
            value = "None"
        else:
            value = python_object.filename
        return [libtbx.phil.tokenizer.word(value=value)]


class ReflectionTableConverters(object):
    """A phil converter for the reflection table class."""

    phil_type = "reflection_table"

    def __str__(self):
        return self.phil_type

    def from_words(self, words, master):
        s = libtbx.phil.str_from_words(words=words)
        if s is None:
            return None
        if not os.path.exists(s):
            raise Sorry("File %s does not exist" % s)
        return FilenameDataWrapper(s, flex.reflection_table.from_file(s))

    def as_words(self, python_object, master):
        if python_object is None:
            value = "None"
        else:
            value = python_object.filename
        return [libtbx.phil.tokenizer.word(value=value)]


class ReflectionTableSelectorConverters(object):
    """A phil converter for the reflection table selector class."""

    phil_type = "reflection_table_selector"

    def __str__(self):
        return self.phil_type

    def from_words(self, words, master):
        s = libtbx.phil.str_from_words(words=words)
        if s is None:
            return None
        regex = r"^\s*([\w\.]+)\s*(<=|!=|==|>=|<|>|&)\s*(.+)\s*$"
        matches = re.findall(regex, s)
        if len(matches) == 0:
            raise RuntimeError("%s is not a valid selection" % s)
        elif len(matches) > 1:
            raise RuntimeError("%s is not a single selection" % s)
        else:
            match = matches[0]
        assert len(match) == 3
        col, op, value = match
        return flex.reflection_table_selector(col, op, value)

    def as_words(self, python_object, master):
        if python_object is None:
            value = "None"
        else:
            value = "%s%s%s" % (
                python_object.column,
                python_object.op_string,
                python_object.value,
            )
        return [libtbx.phil.tokenizer.word(value=value)]


# Create the default converter registry with the extract converters
default_converter_registry = libtbx.phil.extended_converter_registry(
    additional_converters=[
        ExperimentListConverters,
        ReflectionTableConverters,
        ReflectionTableSelectorConverters,
    ]
)


def parse(
    input_string=None,
    source_info=None,
    file_name=None,
    converter_registry=None,
    process_includes=False,
):
    """Redefinition of the parse function."""
    if converter_registry is None:
        converter_registry = default_converter_registry
    return libtbx.phil.parse(
        input_string=input_string,
        source_info=source_info,
        file_name=file_name,
        converter_registry=converter_registry,
        process_includes=process_includes,
    )
