"""
A datastructure for summing over groups of symmetry equivalent reflections.

This module defines a blocked datastructures for summing over groups of
symmetry equivalent reflections, as required for scaling.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from ordered_set import OrderedSet
from scipy.sparse import csc_matrix

from cctbx import crystal, miller, sgtbx, uctbx
from dxtbx import flumpy
from scitbx import sparse

from dials.algorithms.scaling.error_model.error_model import BasicErrorModel
from dials.array_family import flex


def map_indices_to_asu(miller_indices, space_group, anomalous=False):
    """Map the indices to the asymmetric unit."""
    crystal_symmetry = crystal.symmetry(space_group=space_group)
    miller_set = miller.set(
        crystal_symmetry=crystal_symmetry,
        indices=miller_indices,
        anomalous_flag=anomalous,
    )
    miller_set_in_asu = miller_set.map_to_asu()
    return miller_set_in_asu.indices()


def get_sorted_asu_indices(asu_indices, space_group, anomalous=False):
    """Return the sorted asu indices and the permutation selection."""
    crystal_symmetry = crystal.symmetry(space_group=space_group)
    miller_set_in_asu = miller.set(
        crystal_symmetry=crystal_symmetry, indices=asu_indices, anomalous_flag=anomalous
    )
    permuted = miller_set_in_asu.sort_permutation(by_value="packed_indices")
    sorted_asu_miller_index = asu_indices.select(permuted)
    return sorted_asu_miller_index, permuted


class IhTable:
    """
    A class to manage access to Ih_table blocks.

    The idea here is to split the data into blocks to allow parallelized
    computations, but within the blocks the data are sorted by dataset.
    In each block, there exists a block_selection_list which contains the indices
    for each dataset from the input reflection table.

    This class acts as a 'master' to setup the block structure and control access
    to the underlying blocks - only metadata is kept in this class after
    initialisation, the reflections etc are all contained in the blocks.

    Attributes:
        space_group: The space group for the dataset.
        Ih_table_blocks (list): A list of IhTableBlock instances. All symmetry
            equivalent reflections are recorded in the same block, to allow
            splitting of the dataset for parallelized computations.
        nblocks (int): The number of blocks in the Ih_table_blocks list.
        blocked_selection_list (list): A list of lists. bsl[i][j] is the selection
            list for block i, dataset j.
        n_datasets: The number of input reflection tables used to make the Ih_table.
        size: The number of reflections across all blocks
        asu_index_dict (dict): A dictionary, key: asu_miller_index, value tuple
            containing group_id and block_id (where group id is the group index
            within its block).
    """

    id_ = "IhTable"

    def __init__(
        self,
        reflection_tables: list[flex.reflection_table],
        space_group: sgtbx.space_group,
        indices_lists: list[flex.size_t] | None = None,
        nblocks: int = 1,
        free_set_percentage: float = 0,
        free_set_offset: int = 0,
        additional_cols: list[str] | None = None,
        anomalous: bool = False,
    ):
        """
        Distribute the input data into the required structure.

        The reflection data can be split into blocks, while the relevant
        metadata is also generated.

        A list of flex.size_t indices can be provided - this allows the
        reflection table data to maintain a reference to a dataset from which
        it was selected; these will be used when making the block selections.
        e.g selection = flex.bool([True, False, True])
            r_1 = r_master.select(selection)
            indices_list = selection.iselection() = flex.size_t([0, 2])
            then the block selection will contain 0 and 2 to refer back
            to the location of the data in r_master.
        """
        if indices_lists:
            assert len(indices_lists) == len(reflection_tables)
        self.anomalous = anomalous
        self._asu_index_dict = {}
        self._free_asu_index_dict = {}
        self.space_group = space_group
        self.n_work_blocks = nblocks
        self.n_datasets = len(reflection_tables)
        self.Ih_table_blocks = []
        self.blocked_selection_list = []
        self.properties_dict = {
            "n_unique_in_each_block": [],
            "n_reflections_in_each_block": {},
            "miller_index_boundaries": [],
        }
        self.free_set_percentage = free_set_percentage
        self._determine_required_block_structures(
            reflection_tables, free_set_percentage, free_set_offset
        )
        self._create_empty_Ih_table_blocks()
        for i, table in enumerate(reflection_tables):
            if indices_lists:
                self._add_dataset_to_blocks(
                    i, table, indices_lists[i], additional_cols=additional_cols
                )
            else:
                self._add_dataset_to_blocks(i, table, additional_cols=additional_cols)
        self.generate_block_selections()
        self.free_Ih_table = None
        if free_set_percentage > 0:
            self.extract_free_set()
            self.free_Ih_table = True
        self.calc_Ih()

    def update_data_in_blocks(
        self, data: flex.double, dataset_id: int, column: str = "intensity"
    ) -> None:
        """
        Update a given column across all blocks for a given dataset.

        Given an array of data (of the same size as the input reflection
        table) and the name of the column, use the internal data to split
        this up and set in individual blocks.
        """
        assert column in ["intensity", "variance", "inverse_scale_factor"]
        assert dataset_id in range(self.n_datasets)
        # split up data for blocks
        data = flumpy.to_numpy(data)
        for block in self.blocked_data_list:
            data_for_block = data[block.block_selections[dataset_id]]
            start = block.dataset_info[dataset_id]["start_index"]
            end = block.dataset_info[dataset_id]["end_index"]
            block.Ih_table.loc[np.arange(start=start, stop=end), column] = (
                data_for_block
            )

    def get_block_selections_for_dataset(self, dataset: int) -> list[flex.size_t]:
        """Generate the block selection list for a given dataset."""
        assert dataset in range(self.n_datasets)
        if self.free_Ih_table:
            return [
                self.blocked_selection_list[i][dataset]
                for i in range(self.n_work_blocks + 1)
            ]
        return [
            self.blocked_selection_list[i][dataset] for i in range(self.n_work_blocks)
        ]

    @property
    def size(self) -> int:
        """Sum the sizes of all work blocks to give the total number of reflections."""
        if self.free_Ih_table:
            return sum(block.size for block in self.Ih_table_blocks[:-1])
        return sum(block.size for block in self.Ih_table_blocks)

    def generate_block_selections(self) -> None:
        """Generate and set an updated blocked_selection_list."""
        self.blocked_selection_list = [
            block.block_selections for block in self.Ih_table_blocks
        ]

    def update_weights(
        self, error_model: BasicErrorModel | None = None, dataset_id: int = None
    ) -> None:
        """Update the error model in the blocks."""
        for block in self.Ih_table_blocks:
            block.update_weights(error_model, dataset_id)

    @property
    def blocked_data_list(self) -> list[IhTableBlock]:
        """Return the list of IhTableBlock instances."""
        return self.Ih_table_blocks

    def set_derivatives(self, derivatives: sparse.matrix, block_id: int) -> None:
        """Set the derivatives matrix for a given block."""
        self.Ih_table_blocks[block_id].derivatives = derivatives

    def set_inverse_scale_factors(self, new_scales: np.array, block_id: int) -> None:
        """Set the inverse scale factors for a given block."""
        self.Ih_table_blocks[block_id].inverse_scale_factors = new_scales

    def calc_Ih(self, block_id: int = None) -> None:
        """Calculate the latest value of Ih, for a given block or for all blocks."""
        if block_id is not None:
            self.Ih_table_blocks[block_id].calc_Ih()
        else:
            for block in self.Ih_table_blocks:
                block.calc_Ih()

    def _determine_required_block_structures(
        self,
        reflection_tables: list[flex.reflection_table],
        free_set_percentage: float = 0,
        free_set_offset: int = 0,
    ) -> None:
        """
        Inspect the input to determine how to split into blocks.

        Extract the asu miller indices from the reflection table and
        add data to the asu_index_dict and properties dict.
        """
        joint_asu_indices = flex.miller_index()
        for table in reflection_tables:
            if "asu_miller_index" not in table:
                table["asu_miller_index"] = map_indices_to_asu(
                    table["miller_index"], self.space_group, self.anomalous
                )
            joint_asu_indices.extend(table["asu_miller_index"])
        sorted_joint_asu_indices, _ = get_sorted_asu_indices(
            joint_asu_indices, self.space_group, self.anomalous
        )
        if not sorted_joint_asu_indices:
            raise ValueError("No data found in input file(s)")

        asu_index_set = OrderedSet(sorted_joint_asu_indices)
        n_unique_groups = len(asu_index_set)
        n_free_groups = None
        interval_between_free_groups = None
        if free_set_percentage:
            n_free_groups = int(free_set_percentage * n_unique_groups / 100.0)
            n_work_groups = n_unique_groups - n_free_groups
            interval_between_free_groups = int(100 / free_set_percentage)
        else:
            n_work_groups = n_unique_groups
        self.n_work_blocks = min(self.n_work_blocks, n_work_groups)
        # first remove the free set groups
        if free_set_percentage:
            groups_for_free_set = np.full(n_unique_groups, False, dtype=bool)
            for_free = np.arange(
                0 + free_set_offset, n_unique_groups, interval_between_free_groups
            )
            groups_for_free_set[for_free] = True
            asu_index_set = np.array(list(asu_index_set))
            # work_asu_index_set = asu_index_set[~groups_for_free_set]
            free_asu_index_set = asu_index_set[groups_for_free_set]
        else:
            # work_asu_index_set = asu_index_set
            free_asu_index_set = None

        # also record how many unique groups go into each block
        group_boundaries = [
            int(i * n_unique_groups / self.n_work_blocks)
            for i in range(self.n_work_blocks)
        ]
        group_boundaries.append(n_unique_groups)

        next_boundary = group_boundaries[1]
        block_id = 0
        group_id_in_block_i = 0
        for i, index in enumerate(asu_index_set):
            if i == next_boundary:
                self.properties_dict["n_unique_in_each_block"].append(
                    group_id_in_block_i
                )
                self.properties_dict["miller_index_boundaries"].append(tuple(index))
                block_id += 1
                next_boundary = group_boundaries[block_id + 1]
                group_id_in_block_i = 0
            self._asu_index_dict[tuple(index)] = group_id_in_block_i
            group_id_in_block_i += 1
        # record the number in the last work block
        self.properties_dict["n_unique_in_each_block"].append(group_id_in_block_i)
        self.properties_dict["miller_index_boundaries"].append((10000, 10000, 10000))

        block_id += 1
        group_id_in_block_i = 0
        if free_asu_index_set is not None:
            for index in free_asu_index_set:
                # no boundaries as all go into the final block
                self._free_asu_index_dict[tuple(index)] = group_id_in_block_i
                group_id_in_block_i += 1
        # record the number in the free block
        self.properties_dict["n_unique_in_each_block"].append(group_id_in_block_i)
        self.properties_dict["miller_index_boundaries"].append((10000, 10000, 10000))
        # ^ to avoid bounds checking when in last group
        # need to know how many reflections will be in each block also

        block_id = 0
        idx_prev = 0
        boundary = self.properties_dict["miller_index_boundaries"][0]
        for i, index in enumerate(sorted_joint_asu_indices):
            if index == boundary:
                n_in_prev_group = i - idx_prev
                self.properties_dict["n_reflections_in_each_block"][block_id] = (
                    n_in_prev_group
                )
                block_id += 1
                boundary = self.properties_dict["miller_index_boundaries"][block_id]
                idx_prev = i
        self.properties_dict["n_reflections_in_each_block"][block_id] = (
            len(sorted_joint_asu_indices) - idx_prev
        )

    def _create_empty_Ih_table_blocks(self) -> None:
        for n in range(self.n_work_blocks):
            n_refl_in_block = self.properties_dict["n_reflections_in_each_block"][n]
            n_groups_in_block = self.properties_dict["n_unique_in_each_block"][n]
            self.Ih_table_blocks.append(
                IhTableBlock(
                    n_groups=n_groups_in_block,
                    n_refl=n_refl_in_block,
                    n_datasets=self.n_datasets,
                )
            )

    def _add_dataset_to_blocks(
        self,
        dataset_id: int,
        reflections: flex.reflection_table,
        indices_array: flex.size_t | None = None,
        additional_cols: list[str] | None = None,
    ) -> None:
        sorted_asu_indices, perm = get_sorted_asu_indices(
            reflections["asu_miller_index"], self.space_group, self.anomalous
        )
        hkl = reflections["asu_miller_index"]
        df = pd.DataFrame()
        df["intensity"] = flumpy.to_numpy(reflections["intensity"])
        df["variance"] = flumpy.to_numpy(reflections["variance"])
        df["inverse_scale_factor"] = flumpy.to_numpy(
            reflections["inverse_scale_factor"]
        )
        if isinstance(additional_cols, list):
            for col in additional_cols:
                if col in reflections:
                    df[col] = flumpy.to_numpy(reflections[col])
        if indices_array:
            df["loc_indices"] = flumpy.to_numpy(indices_array)
        else:
            df["loc_indices"] = np.arange(df.shape[0], dtype=np.uint64)
        df = df.iloc[flumpy.to_numpy(perm)]
        hkl = hkl.select(perm)
        df["dataset_id"] = np.full(df.shape[0], dataset_id, dtype=np.uint64)
        # if data are sorted by asu_index, then up until boundary, should be in same
        # block (still need to read group_id though)
        # sort data, get group ids and block_ids
        group_ids = np.zeros(sorted_asu_indices.size(), dtype=np.uint64)
        boundary = self.properties_dict["miller_index_boundaries"][0]
        boundary_id = 0
        boundaries_for_this_datset = [0]  # use to slice
        # make this a c++ method for speed?
        prev = (0, 0, 0)
        group_id = -1
        for i, index in enumerate(sorted_asu_indices):
            if index != prev:
                while index >= boundary:
                    boundaries_for_this_datset.append(i)
                    boundary_id += 1
                    boundary = self.properties_dict["miller_index_boundaries"][
                        boundary_id
                    ]
                group_id = self._asu_index_dict[tuple(index)]
                prev = index
            group_ids[i] = group_id
        while len(boundaries_for_this_datset) < self.n_work_blocks + 1:
            # catch case where last boundaries aren't reached
            boundaries_for_this_datset.append(len(sorted_asu_indices))
        # so now have group ids as well for individual dataset
        if self.n_work_blocks == 1:
            self.Ih_table_blocks[0].add_data(dataset_id, group_ids, df, hkl)
        else:
            for i, val in enumerate(boundaries_for_this_datset[:-1]):
                start = val
                end = boundaries_for_this_datset[i + 1]
                self.Ih_table_blocks[i].add_data(
                    dataset_id, group_ids[start:end], df[start:end], hkl[start:end]
                )

    def extract_free_set(self) -> None:
        """Extract a free set from all blocks."""
        assert not self.free_Ih_table
        free_reflection_table = pd.DataFrame()
        free_indices = np.array([], dtype=int).reshape((0,))
        free_hkl = flex.miller_index([])
        # for each block, remove a fraction of the groups
        for j, block in enumerate(self.Ih_table_blocks):
            n_groups = block.n_groups
            groups_for_free_set = np.full(n_groups, False, dtype=bool)
            for_free = np.array(
                [
                    tuple(i) in self._free_asu_index_dict
                    for i in OrderedSet(block.asu_miller_index)
                ]
            )
            groups_for_free_set[for_free] = True
            free_block = block.select_on_groups(groups_for_free_set)
            free_reflection_table = pd.concat(
                [free_reflection_table, free_block.Ih_table]
            )
            free_hkl.extend(free_block.asu_miller_index)
            for sel in free_block.block_selections:
                free_indices = np.concatenate([free_indices, sel])
            self.Ih_table_blocks[j] = block.select_on_groups(~groups_for_free_set)
            # Now need to update dataset_info dict.
            removed_from_each_dataset = [
                np.count_nonzero(free_block.Ih_table["dataset_id"].to_numpy() == i)
                for i in range(0, block.n_datasets)
            ]
            n_removed = 0
            for i in range(0, self.Ih_table_blocks[j].n_datasets):
                self.Ih_table_blocks[j].dataset_info[i]["start_index"] -= n_removed
                n_removed += removed_from_each_dataset[i]
                self.Ih_table_blocks[j].dataset_info[i]["end_index"] -= n_removed
        self.blocked_selection_list = [
            block.block_selections for block in self.Ih_table_blocks
        ]
        # now split by dataset and use to instantiate another Ih_table
        datasets = set(free_reflection_table["dataset_id"])
        tables = []
        indices_lists = []
        n_refl = 0
        for id_ in datasets:
            dataset_sel = free_reflection_table["dataset_id"].to_numpy() == id_
            n_refl += np.count_nonzero(dataset_sel)
            tables.append(free_reflection_table[dataset_sel])
            indices_lists.append(free_indices[dataset_sel])
        free_block = IhTableBlock(
            n_groups=len(set(free_hkl)), n_refl=n_refl, n_datasets=len(datasets)
        )
        group_ids = np.array(
            [self._free_asu_index_dict[tuple(index)] for index in free_hkl],
            dtype=np.uint64,
        )
        for id_, t in zip(datasets, tables):
            dataset_sel = free_reflection_table["dataset_id"].to_numpy() == id_
            group_id_this = group_ids[dataset_sel]
            hkl_this = free_hkl.select(flumpy.from_numpy(dataset_sel))
            free_block.add_data(id_, group_id_this, t, hkl_this)

        self.Ih_table_blocks.append(free_block)
        self.blocked_selection_list.append(free_block.block_selections)

    def as_miller_array(
        self, unit_cell: uctbx.unit_cell, return_free_set_data: bool = False
    ) -> miller.array:
        """Get a scaled miller array from the Ih_table and an experiment."""
        blocked_data_list = self.blocked_data_list
        joint_table = flex.reflection_table([])
        if self.free_Ih_table:
            if return_free_set_data:
                blocked_data_list = [blocked_data_list[-1]]
            else:
                blocked_data_list = blocked_data_list[:-1]
        if len(blocked_data_list) > 1:
            for block in blocked_data_list:
                joint_table.extend(block.as_reflection_table())
        else:
            joint_table = blocked_data_list[0].as_reflection_table()
        # Filter out negative scale factors to avoid merging statistics errors.
        return _reflection_table_to_iobs(joint_table, unit_cell, self.space_group)


class TargetAsuDictCache:
    instances = {}

    def __new__(cls, target_Ih_table):
        id_ = id(target_Ih_table)
        if id_ not in cls.instances:
            cls.instances[id_] = dict(
                zip(
                    target_Ih_table.blocked_data_list[0].asu_miller_index,
                    target_Ih_table.blocked_data_list[0].Ih_values,
                )
            )
        return cls.instances[id_]


class IhTableBlock:
    """
    A datastructure for efficient summations over symmetry equivalent reflections.

    This contains a reflection table, sorted by dataset, called the Ih_table,
    a h_index_matrix (sparse) for efficiently calculating sums over symmetry
    equivalent reflections as well as 'block_selections' which relate the order
    of the data to the initial reflection tables used to initialise the (master)
    IhTable.

    Attributes:
        Ih_table: A reflection table, containing I, g, w, var, Ih,
            asu_miller_index, loc_indices and dataset_id.
        block_selections: A list of flex.size_t arrays of indices, that can be
            used to select and reorder data from the input reflection tables to
            match the order in the Ih_table.
        h_index_matrix: A sparse matrix used to sum over groups of equivalent
            reflections by multiplication. Sum_h I = I * h_index_matrix. The
            dimension is n_refl by n_groups; each row has a single nonzero
            entry with a value of 1.
        h_expand_matrix: The transpose of the h_index_matrix, used to expand an
            array of values for symmetry groups into an array of size n_refl.
        derivatives: A matrix of derivatives of the reflections wrt the model
            parameters.
    """

    def __init__(self, n_groups: int, n_refl: int, n_datasets: int = 1):
        """Create empty datastructures to which data can later be added."""
        self.Ih_table = pd.DataFrame()
        self.block_selections = [None] * n_datasets
        self.h_index_matrix = sparse.matrix(n_refl, n_groups)
        self._setup_info = {"next_row": 0, "next_dataset": 0, "setup_complete": False}
        self.dataset_info = {}
        self.n_datasets = n_datasets
        self.h_expand_matrix = None
        self.derivatives = None
        self.binner = None
        self._csc_rows = np.array([], dtype=np.uint64).reshape((0,))
        self._csc_cols = np.array([], dtype=np.uint64).reshape((0,))
        self._csc_h_index_matrix = None
        self._csc_h_expand_matrix = None
        self._hkl = flex.miller_index([])

    def add_data(
        self,
        dataset_id: int,
        group_ids: np.array,
        reflections: pd.DataFrame,
        hkl: flex.miller_index,
    ) -> None:
        """
        Add data to all blocks for a given dataset.

        Add data to the Ih_table, write data to the h_index_matrix and
        add the loc indices to the block_selections list.
        """
        assert not self._setup_info["setup_complete"], """
No further data can be added to the IhTableBlock as setup marked complete."""
        assert (
            self._setup_info["next_row"] + len(group_ids) <= self.h_index_matrix.n_rows
        ), """
Not enough space left to add this data, please check for correct block initialisation."""
        assert dataset_id == self._setup_info["next_dataset"], """
Datasets must be added in correct order: expected: {}, this dataset: {}""".format(
            self._setup_info["next_dataset"],
            dataset_id,
        )
        for i, id_ in enumerate(group_ids):
            rowidx = i + self._setup_info["next_row"]
            self.h_index_matrix[rowidx, int(id_)] = 1.0
        cols = group_ids
        rows = np.arange(
            start=self._setup_info["next_row"],
            stop=self._setup_info["next_row"] + group_ids.size,
            dtype=np.uint64,
        )

        self._csc_cols = np.concatenate([self._csc_cols, cols])
        self._csc_rows = np.concatenate([self._csc_rows, rows])
        self._hkl.extend(hkl)

        self.dataset_info[dataset_id] = {"start_index": self._setup_info["next_row"]}
        self._setup_info["next_row"] += len(group_ids)
        self._setup_info["next_dataset"] += 1
        self.dataset_info[dataset_id]["end_index"] = self._setup_info["next_row"]
        self.Ih_table = pd.concat([self.Ih_table, reflections], ignore_index=True)
        if "loc_indices" in reflections:
            self.block_selections[dataset_id] = reflections["loc_indices"].to_numpy()
        else:
            self.block_selections[dataset_id] = np.arange(
                reflections.shape[0], dtype=np.uint64
            )
        if self._setup_info["next_dataset"] == len(self.block_selections):
            self._complete_setup()

    def _complete_setup(self) -> None:
        """Finish the setup of the Ih_table once all data has been added."""
        self.h_index_matrix.compact()
        assert self._setup_info["next_row"] == self.h_index_matrix.n_rows, """
Not all rows of h_index_matrix appear to be filled in IhTableBlock setup."""
        self.h_expand_matrix = self.h_index_matrix.transpose()
        data = np.full(self._csc_cols.size, 1.0)
        self._csc_h_index_matrix = csc_matrix((data, (self._csc_rows, self._csc_cols)))
        self._csc_h_expand_matrix = self._csc_h_index_matrix.transpose()
        self.weights = 1.0 / self.variances
        self._setup_info["setup_complete"] = True

    def group_multiplicities(self, output: str = "per_group") -> np.array:
        """Return the multiplicities of the symmetry groups."""
        return self.sum_in_groups(np.full(self.size, 1.0), output=output)

    def select(self, sel: np.array) -> IhTableBlock:
        """Select a subset of the data, returning a new IhTableBlock object."""
        Ih_table = self.Ih_table[sel]
        Ih_table.reset_index(drop=True, inplace=True)
        h_idx_sel = self.h_expand_matrix.select_columns(
            flumpy.from_numpy(sel).iselection()
        )
        reduced_h_idx = h_idx_sel.transpose()
        unity = flex.double(int(Ih_table.size), 1.0)
        nz_col_sel = (unity * reduced_h_idx) > 0
        h_index_matrix = reduced_h_idx.select_columns(nz_col_sel.iselection())
        h_expand = h_index_matrix.transpose()
        csc_h_idx_sel = self._csc_h_expand_matrix[:, sel]
        csc_h_index_matrix = csc_h_idx_sel.transpose()[:, flumpy.to_numpy(nz_col_sel)]
        csc_h_expand_matrix = csc_h_index_matrix.transpose()
        newtable = IhTableBlock(n_groups=0, n_refl=0, n_datasets=self.n_datasets)
        newtable.Ih_table = Ih_table
        newtable._hkl = self._hkl.select(flumpy.from_numpy(sel))
        newtable.h_expand_matrix = h_expand
        newtable.h_index_matrix = h_index_matrix
        newtable._csc_h_index_matrix = csc_h_index_matrix
        newtable._csc_h_expand_matrix = csc_h_expand_matrix
        newtable.block_selections = []
        offset = 0
        for i in range(newtable.n_datasets):
            newtable.dataset_info[i] = {"start_index": offset}
            block_sel_i = self.block_selections[i]
            n_in_dataset_i = len(block_sel_i)
            newtable.block_selections.append(
                block_sel_i[sel[offset : offset + n_in_dataset_i]]
            )
            offset += n_in_dataset_i
            newtable.dataset_info[i]["end_index"] = offset
        return newtable

    def select_on_groups(self, sel: np.array) -> IhTableBlock:
        """Select a subset of the unique groups, returning a new IhTableBlock."""
        reduced_h_idx = self._csc_h_index_matrix[:, sel]
        unity = np.full(reduced_h_idx.shape[1], 1.0)
        nz_row_sel = (unity * reduced_h_idx.transpose()) > 0
        return self.select(nz_row_sel)

    def calc_Ih(self) -> None:
        """Calculate the current best estimate for Ih for each reflection group."""
        scale_factors = self.inverse_scale_factors
        sumgsq = self.sum_in_groups(np.square(scale_factors) * self.weights)
        sumgI = self.sum_in_groups(scale_factors * self.intensities * self.weights)
        Ih = sumgI / sumgsq
        self.Ih_table.loc[:, "Ih_values"] = Ih @ self._csc_h_expand_matrix

    def update_weights(
        self,
        error_model: BasicErrorModel | None = None,
        dataset_id: int | None = None,
    ) -> None:
        """Update the scaling weights based on an error model."""
        if error_model:
            if dataset_id is not None:  # note the first dataset has an id of 0
                sel = self.Ih_table["dataset_id"].to_numpy() == dataset_id
                sigmaprimesq = error_model.update_variances(
                    self.variances[sel], self.intensities[sel]
                )
                self.Ih_table.loc[sel, "weights"] = 1.0 / sigmaprimesq
            else:
                sigmaprimesq = error_model.update_variances(
                    self.variances, self.intensities
                )
                self.Ih_table.loc[:, "weights"] = 1.0 / sigmaprimesq
        else:
            if dataset_id is not None:  # note the first dataset has an id of 0
                sel = self.Ih_table["dataset_id"].to_numpy() == dataset_id
                self.Ih_table.loc[sel, "weights"] = 1.0 / self.variances[sel]
            else:
                self.Ih_table.loc[:, "weights"] = 1.0 / self.variances

    def calc_nh(self) -> np.array:
        """Calculate the number of refls in the group to which the reflection belongs.

        This is a vector of length n_refl."""
        return self.sum_in_groups(np.full(self.size, 1.0), output="per_refl")

    def match_Ih_values_to_target(self, target_Ih_table: IhTable) -> None:
        """
        Use an Ih_table as a target to set Ih values in this table.

        Given an Ih table as a target, the common reflections across the tables
        are determined and the Ih_values are set to those of the target. If no
        matching reflection is found, then the values are removed from the table.
        """
        assert target_Ih_table.n_work_blocks == 1
        target_asu_Ih_dict = TargetAsuDictCache(target_Ih_table)
        new_Ih_values = np.zeros(self.size, dtype=float)
        location_in_unscaled_array = 0
        sorted_asu_indices, permuted = get_sorted_asu_indices(
            self.asu_miller_index,
            target_Ih_table.space_group,
            anomalous=target_Ih_table.anomalous,
        )
        for j, miller_idx in enumerate(OrderedSet(sorted_asu_indices)):
            n_in_group = self._csc_h_index_matrix.getcol(j).count_nonzero()
            if miller_idx in target_asu_Ih_dict:
                i = location_in_unscaled_array
                new_Ih_values[np.arange(i, i + n_in_group, dtype=np.uint64)] = np.full(
                    n_in_group, target_asu_Ih_dict[miller_idx]
                )
            location_in_unscaled_array += n_in_group
        self.Ih_table.loc[flumpy.to_numpy(permuted), "Ih_values"] = new_Ih_values
        sel = self.Ih_values != 0.0
        new_table = self.select(sel)
        # now set attributes to update object
        self.Ih_table = new_table.Ih_table
        self.h_index_matrix = new_table.h_index_matrix
        self.h_expand_matrix = new_table.h_expand_matrix
        self.block_selections = new_table.block_selections
        self._csc_h_expand_matrix = new_table._csc_h_expand_matrix
        self._csc_h_index_matrix = new_table._csc_h_index_matrix

    @property
    def inverse_scale_factors(self) -> np.array:
        """The inverse scale factors of the reflections."""
        return self.Ih_table["inverse_scale_factor"].to_numpy()

    @inverse_scale_factors.setter
    def inverse_scale_factors(self, new_scales: np.array) -> None:
        if new_scales.size != self.size:
            assert 0, f"""attempting to set a new set of scale factors of different
      length than previous assignment: was {self.inverse_scale_factors.size}, attempting {new_scales.size}"""
        else:
            self.Ih_table.loc[:, "inverse_scale_factor"] = new_scales

    @property
    def variances(self) -> np.array:
        """The variances of the reflections."""
        return self.Ih_table["variance"].to_numpy()

    @variances.setter
    def variances(self, new_variances: np.array) -> None:
        assert new_variances.size == self.size
        self.Ih_table.loc[:, "variance"] = new_variances

    @property
    def intensities(self) -> np.array:
        """The unscaled reflection intensities."""
        return self.Ih_table["intensity"].to_numpy()

    @intensities.setter
    def intensities(self, new_intensities):
        assert new_intensities.size == self.size
        self.Ih_table.loc[:, "intensity"] = new_intensities

    @property
    def Ih_values(self) -> np.array:
        """The bset-estimated intensities of symmetry equivalent reflections."""
        return self.Ih_table["Ih_values"].to_numpy()

    @property
    def weights(self) -> np.array:
        """The weights that will be used in scaling."""
        return self.Ih_table["weights"].to_numpy()

    @weights.setter
    def weights(self, new_weights):
        if new_weights.size != self.size:
            assert 0, f"""attempting to set a new set of weights of different
      length than previous assignment: was {self.size}, attempting {new_weights.size}"""
        self.Ih_table.loc[:, "weights"] = new_weights

    @property
    def size(self) -> int:
        """Return the length of the stored Ih_table (a reflection table)."""
        return self.Ih_table.shape[0]

    @property
    def n_groups(self) -> int:
        """Return the length of the stored Ih_table (a reflection table)."""
        return self._csc_h_index_matrix.shape[1]

    @property
    def asu_miller_index(self) -> flex.miller_index:
        """Return the miller indices in the asymmetric unit."""
        return self._hkl

    def setup_binner(
        self,
        unit_cell: uctbx.unit_cell,
        space_group: sgtbx.space_group,
        n_resolution_bins: int,
    ) -> None:
        """Create a binner for the reflections contained in the table."""
        ma = _reflection_table_to_iobs(
            self.as_reflection_table(), unit_cell, space_group
        )
        # need d star sq step
        d_star_sq = ma.d_star_sq().data()
        d_star_sq_min = flex.min(d_star_sq)
        d_star_sq_max = flex.max(d_star_sq)
        span = d_star_sq_max - d_star_sq_min
        relative_tolerance = 1e-6
        d_star_sq_max += span * relative_tolerance
        d_star_sq_min -= span * relative_tolerance
        # Avoid a zero-size step that would otherwise anger the d_star_sq_step binner.
        step = max((d_star_sq_max - d_star_sq_min) / n_resolution_bins, 0.004)

        self.binner = ma.setup_binner_d_star_sq_step(
            auto_binning=False,
            d_max=uctbx.d_star_sq_as_d(d_star_sq_max),
            d_min=uctbx.d_star_sq_as_d(d_star_sq_min),
            d_star_sq_step=step,
        )

    def sum_in_groups(
        self, array: csc_matrix | np.array, output: str = "per_group"
    ) -> np.array:
        """
        Sums an array object over the symmetry equivalent groups.
        The array's final dimension must equal the size of the Ih_table.
        """
        if output == "per_group":
            return array @ self._csc_h_index_matrix
        elif output == "per_refl":  # return the summed quantity per reflection
            return (array @ self._csc_h_index_matrix) @ self._csc_h_expand_matrix
        else:
            raise ValueError(
                f"""Bad value for output= parameter
(value={output}, allowed values: per_group, per_refl)"""
            )

    def as_reflection_table(self) -> flex.reflection_table:
        """Return the data in flex reflection table format"""
        table = flex.reflection_table()
        table["asu_miller_index"] = self.asu_miller_index
        for k, v in self.Ih_table.items():
            table[k] = flumpy.from_numpy(v.to_numpy())
        return table


def _reflection_table_to_iobs(
    table: flex.reflection_table,
    unit_cell: uctbx.unit_cell,
    space_group: sgtbx.space_group,
) -> miller.array:
    miller_set = miller.set(
        crystal_symmetry=crystal.symmetry(
            unit_cell=unit_cell,
            space_group=space_group,
            assert_is_compatible_unit_cell=False,
        ),
        indices=table["asu_miller_index"],
        anomalous_flag=False,
    )
    i_obs = miller.array(
        miller_set, data=table["intensity"] / table["inverse_scale_factor"]
    )
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(flex.sqrt(table["variance"]) / table["inverse_scale_factor"])
    i_obs.set_info(miller.array_info(source="DIALS", source_type="reflection_tables"))
    return i_obs
