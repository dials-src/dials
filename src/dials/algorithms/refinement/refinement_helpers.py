"""Auxiliary functions for the refinement package"""

from __future__ import annotations

import logging
import math
import random

import scitbx.matrix
from scitbx.array_family import flex

from dials_refinement_helpers_ext import CrystalOrientationCompose as xloc_cpp
from dials_refinement_helpers_ext import PanelGroupCompose as pgc_cpp
from dials_refinement_helpers_ext import dR_from_axis_and_angle as dR_cpp

logger = logging.getLogger(__name__)


def ordinal_number(array_index=None, cardinal_number=None):
    """Return a string representing the ordinal number for the input integer. One
    of array_index or cardinal_number must be set, depending on whether the
    input is from a 0-based or 1-based sequence.

    Based on Thad Guidry's post at
    https://groups.google.com/forum/#!topic/openrefine/G7_PSdUeno0"""
    if [array_index, cardinal_number].count(None) != 1:
        raise ValueError("One of array_index or cardinal_number should be set")
    if array_index is not None:
        i = int(array_index) + 1
    if cardinal_number is not None:
        i = int(cardinal_number)
    return str(i) + {1: "st", 2: "nd", 3: "rd"}.get(
        4 if 10 <= i % 100 < 20 else i % 10, "th"
    )


class PanelGroupCompose(pgc_cpp):
    """Wrapper for the C++ PanelGroupCompose class with accessors that
    return scitbx matrix values."""

    def d1(self):
        return scitbx.matrix.col(super().d1())

    def d2(self):
        return scitbx.matrix.col(super().d2())

    def origin(self):
        return scitbx.matrix.col(super().origin())

    def derivatives_for_panel(self, offset, dir1_new_basis, dir2_new_basis):
        d = super().derivatives_for_panel(offset, dir1_new_basis, dir2_new_basis)
        return [scitbx.matrix.sqr(e) for e in d]


class CrystalOrientationCompose(xloc_cpp):
    """Wrapper for the C++ CrystalOrientationCompose class with accessors that
    return matrix.sqr values."""

    def U(self):
        return scitbx.matrix.sqr(super().U())

    def dU_dphi1(self):
        return scitbx.matrix.sqr(super().dU_dphi1())

    def dU_dphi2(self):
        return scitbx.matrix.sqr(super().dU_dphi2())

    def dU_dphi3(self):
        return scitbx.matrix.sqr(super().dU_dphi3())


def dR_from_axis_and_angle(axis, angle, deg=False):
    """Wrapper for C++ version of dR_from_axis_and_angle returning a matrix.sqr"""
    return scitbx.matrix.sqr(dR_cpp(axis, angle, deg))


def dR_from_axis_and_angle_py(axis, angle, deg=False):
    """return the first derivative of a rotation matrix specified by its
    axis and angle"""

    # NB it is inefficient to do this separately from the calculation of
    # the rotation matrix itself, but it seems the Python interface to
    # scitbx does not have a suitable function. It might perhaps be
    # useful to write one, which could come straight from David Thomas'
    # RTMATS (present in Mosflm and MADNES).

    # NB RTMATS does calculation for a clockwise rotation of a vector
    # whereas axis_and_angle_as_r3_rotation_matrix does anticlockwise
    # rotation. Therefore flip the axis in here compared with
    # RTMATS in order to match the axis_and_angle_as_r3_rotation_matrix
    # convention

    # See also axis_and_angle_as_r3_derivative_wrt_angle, which does the same
    # as this function, but this function is faster.

    assert axis.n in ((3, 1), (1, 3))
    if deg:
        angle *= math.pi / 180
    axis = -1.0 * axis.normalize()
    ca, sa = math.cos(angle), math.sin(angle)

    return scitbx.matrix.sqr(
        (
            sa * axis[0] * axis[0] - sa,
            sa * axis[0] * axis[1] + ca * axis[2],
            sa * axis[0] * axis[2] - ca * axis[1],
            sa * axis[1] * axis[0] - ca * axis[2],
            sa * axis[1] * axis[1] - sa,
            sa * axis[1] * axis[2] + ca * axis[0],
            sa * axis[2] * axis[0] + ca * axis[1],
            sa * axis[2] * axis[1] - ca * axis[0],
            sa * axis[2] * axis[2] - sa,
        )
    )


def skew_symm(v):
    """Make matrix [v]_x from v. Essentially multiply vector by SO(3) basis
    set Lx, Ly, Lz. Equation (2) from Gallego & Yezzi paper.

    NB a C++ version exists in gallego_yezzi.h."""

    L1 = scitbx.matrix.sqr((0, 0, 0, 0, 0, -1, 0, 1, 0))
    L2 = scitbx.matrix.sqr((0, 0, 1, 0, 0, 0, -1, 0, 0))
    L3 = scitbx.matrix.sqr((0, -1, 0, 1, 0, 0, 0, 0, 0))

    v1, v2, v3 = v.elems

    return v1 * L1 + v2 * L2 + v3 * L3


def dRq_de(theta, e, q):
    """Calculate derivative of rotated vector r = R*q with respect to the elements
    of the rotation axis e, where the angle of rotation is theta.

    Implementation of Equation (8) from Gallego & Yezzi.

    NB a C++ version exists in gallego_yezzi.h."""

    # ensure e is unit vector
    e = e.normalize()

    # rotation matrix
    R = e.axis_and_angle_as_r3_rotation_matrix(theta, deg=False)

    # rotation vector v
    v = theta * e

    qx = skew_symm(q)
    vx = skew_symm(v)
    vvt = v * v.transpose()
    Rt = R.transpose()
    I3 = scitbx.matrix.identity(3)

    return (-1.0 / theta) * R * qx * (vvt + (Rt - I3) * vx)


def random_param_shift(vals, sigmas):
    """Add a random (normal) shift to a parameter set, for testing"""

    assert len(vals) == len(sigmas)
    shifts = [random.gauss(0, sd) for sd in sigmas]
    newvals = [(x + y) for x, y in zip(vals, shifts)]

    return newvals


def get_fd_gradients(mp, deltas, multi_state_elt=None):
    """Calculate centered finite difference gradients for each of the
    parameters of the model parameterisation mp.

    "deltas" must be a sequence of the same length as the parameter list, and
    contains the step size for the difference calculations for each parameter.

    "multi_state_elt" selects a particular state for use when mp is a multi-
    state parameterisation.
    """

    p_vals = mp.get_param_vals()
    assert len(deltas) == len(p_vals)
    fd_grad = []

    for i in range(len(deltas)):
        val = p_vals[i]

        p_vals[i] -= deltas[i] / 2.0
        mp.set_param_vals(p_vals)
        if multi_state_elt is None:
            rev_state = mp.get_state()
        else:
            rev_state = mp.get_state(multi_state_elt=multi_state_elt)

        p_vals[i] += deltas[i]
        mp.set_param_vals(p_vals)
        if multi_state_elt is None:
            fwd_state = mp.get_state()
        else:
            fwd_state = mp.get_state(multi_state_elt=multi_state_elt)

        fd_grad.append((fwd_state - rev_state) / deltas[i])

        p_vals[i] = val

    # return to the initial state
    mp.set_param_vals(p_vals)

    return fd_grad


def get_panel_groups_at_depth(group, depth=0):
    """Return a list of the panel groups at a certain depth below the node group"""
    assert depth >= 0
    if depth == 0:
        return [group]
    else:
        assert group.is_group()
        return [
            p
            for gp in group.children()
            for p in get_panel_groups_at_depth(gp, depth - 1)
        ]


def get_panel_ids_at_root(panel_list, group):
    """Get the sequential panel IDs for a set of panels belonging to a group"""
    if group.is_group():
        return [
            p for gp in group.children() for p in get_panel_ids_at_root(panel_list, gp)
        ]
    else:  # we got down to Panels
        return [panel_list.index(group)]


def string_sel(l, full_names, prefix=""):
    """Provide flexible matching between a list of input strings, l,
    consisting either of indices or partial names, and a list of full names,
    with an optional shared prefix. The input list l may arrive from PHIL
    conversion of the strings type. In that case, comma-separated values will
    require splitting, and bracket characters will be removed. The values in
    the processed list l should consist of integers or partial names. Integers
    will be treated as 0-based indices and partial names will be matched to
    as many full names as possible. The created selection is returned as a
    boolean list."""

    sel = [False] * len(full_names)
    full_names = [prefix + s for s in full_names]

    # expand elements of the list that are comma separated strings and remove
    # braces/brackets
    l = (s.strip("(){}[]") for e in l for s in str(e).split(","))
    l = (e for e in l if e != "")
    for e in l:
        try:
            i = int(e)
            sel[i] = True
            continue
        except ValueError:
            pass
        except IndexError:
            pass
        sel = [True if e in name else s for (name, s) in zip(full_names, sel)]

    return sel


def calculate_frame_numbers(reflections, experiments):
    """calculate observed frame numbers for all reflections, if not already
    set"""

    # Only do this if we have to
    if "xyzobs.px.value" in reflections:
        return reflections

    # Ok, frames are not set, so set them, with dummy observed pixel values
    frames = flex.double(len(reflections), 0.0)
    for iexp, exp in enumerate(experiments):
        scan = exp.scan
        if not scan:
            continue
        sel = reflections["id"] == iexp
        xyzobs = reflections["xyzobs.mm.value"].select(sel)
        angles = xyzobs.parts()[2]
        to_update = scan.get_array_index_from_angle(angles, deg=False)
        frames.set_selected(sel, to_update)
    reflections["xyzobs.px.value"] = flex.vec3_double(
        flex.double(len(reflections), 0.0), flex.double(len(reflections), 0.0), frames
    )

    return reflections


def set_obs_s1(reflections, experiments):
    """Set observed s1 vectors for reflections if required, return the number
    of reflections that have been set."""

    refs_wo_s1_sel = reflections["s1"].norms() < 1.0e-6
    nrefs_wo_s1 = refs_wo_s1_sel.count(True)
    if nrefs_wo_s1 == 0:
        return nrefs_wo_s1

    for i_expt, expt in enumerate(experiments):
        detector = expt.detector
        beam = expt.beam
        expt_sel = reflections["id"] == i_expt
        for i_panel, panel in enumerate(detector):
            panel_sel = reflections["panel"] == i_panel
            isel = (expt_sel & panel_sel & refs_wo_s1_sel).iselection()
            spots = reflections.select(isel)
            x, y, rot_angle = spots["xyzobs.mm.value"].parts()
            s1 = panel.get_lab_coord(flex.vec2_double(x, y))
            s1 = s1 / s1.norms() * (1 / beam.get_wavelength())
            reflections["s1"].set_selected(isel, s1)
    return nrefs_wo_s1


def compute_radial_and_transverse_residuals(experiments, reflections, two_theta=False):
    """Compute radial and transverse residuals, relative to the beam center"""
    beam_centre_lab = flex.vec3_double(len(reflections))
    obs_lab_coords = flex.vec3_double(len(reflections))
    calc_lab_coords = flex.vec3_double(len(reflections))
    delta_lab_coords = flex.vec3_double(len(reflections))
    delta_xy_panel = reflections["xyzcal.mm"] - reflections["xyzobs.mm.value"]
    panel_ori = flex.vec3_double(len(reflections))

    all_s0 = flex.vec3_double(len(reflections))
    for expt_id, expt in enumerate(experiments):
        all_s0.set_selected(reflections["id"] == expt_id, expt.beam.get_s0())

    # faster to do this row by row than using selections for multi-panel, multi-experiment data
    for i in range(len(reflections)):
        expt_id = reflections["id"][i]
        panel = experiments[expt_id].detector[reflections["panel"][i]]

        # Compute the beam center in lab space (a vector pointing from the origin to where the beam would intersect
        # the panel, if it did intersect the panel)
        beam_centre_lab[i] = panel.get_beam_centre_lab(all_s0[i])

        # Compute obs in lab space
        obs_lab_coords[i] = panel.get_lab_coord(reflections["xyzobs.mm.value"][i][0:2])

        # Compute deltaXY in lab space
        delta_lab_coords[i] = panel.get_lab_coord(delta_xy_panel[i][0:2])
        panel_ori[i] = panel.get_origin()

        # Compute calc in lab space
        calc_lab_coords[i] = panel.get_lab_coord(reflections["xyzcal.mm"][i][0:2])

    # subtract off the panel origin to center the deltas relative to the lab origin
    delta_lab_coords = delta_lab_coords - panel_ori

    if two_theta:
        tto = obs_lab_coords.angle(all_s0)
        ttc = calc_lab_coords.angle(all_s0)

        d = all_s0.norms() / 2 / flex.sin(ttc / 2)
        reflections["delpsidstar"] = reflections["delpsical.rad"] * flex.sqrt(4.0 / d)

    # The radial vector points from the center of the reflection to the beam center
    radial_vectors = (obs_lab_coords - beam_centre_lab).each_normalize()
    # The transverse vector is orthogonal to the radial vector and the beam vector
    transverse_vectors = radial_vectors.cross(beam_centre_lab).each_normalize()
    # Compute the raidal and transverse components of each deltaXY
    reflections["t_resid"] = delta_lab_coords.dot(transverse_vectors)
    if two_theta:
        reflections["twotheta_resid"] = tto - ttc
    else:
        reflections["r_resid"] = delta_lab_coords.dot(radial_vectors)
