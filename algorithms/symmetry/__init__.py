from __future__ import division, absolute_import, print_function

import logging
logger = logging.getLogger(__name__)

from cStringIO import StringIO

import libtbx
from scitbx.array_family import flex
from cctbx import adptbx
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from mmtbx import scaling
from mmtbx.scaling import absolute_scaling
from mmtbx.scaling import matthews


class symmetry_base(object):
  def __init__(self, intensities,
               normalisation='ml_aniso',
               lattice_symmetry_max_delta=2.0,
               d_min=libtbx.Auto,
               min_i_mean_over_sigma_mean=4,
               min_cc_half=0.6,
               relative_length_tolerance=None,
               absolute_angle_tolerance=None):

    self.input_intensities = intensities

    uc_params = [flex.double() for i in range(6)]
    for d in self.input_intensities:
      for i, p in enumerate(d.unit_cell().parameters()):
        uc_params[i].append(p)
    self.median_unit_cell = uctbx.unit_cell(parameters=[flex.median(p) for p in uc_params])
    for d in self.input_intensities:
      if (relative_length_tolerance is not None and
          absolute_angle_tolerance is not None):
        assert d.unit_cell().is_similar_to(
          self.median_unit_cell, relative_length_tolerance,
          absolute_angle_tolerance), (
            str(d.unit_cell()), str(self.median_unit_cell))

    self.intensities = self.input_intensities[0]
    self.dataset_ids = flex.double(self.intensities.size(), 0)
    for i, d in enumerate(self.input_intensities[1:]):
      self.intensities = self.intensities.concatenate(
        d, assert_is_similar_symmetry=False)
      self.dataset_ids.extend(flex.double(d.size(), i+1))
    self.intensities = self.intensities.customized_copy(
      unit_cell=self.median_unit_cell).set_info(
          self.intensities.info())
    self.intensities.set_observation_type_xray_intensity()

    self.cb_op_inp_min = self.intensities.change_of_basis_op_to_niggli_cell()
    self.intensities = self.intensities.change_basis(
      self.cb_op_inp_min).customized_copy(
        space_group_info=sgtbx.space_group_info('P1')).map_to_asu().set_info(
          self.intensities.info())

    self.lattice_symmetry_max_delta = lattice_symmetry_max_delta
    self.subgroups = metric_subgroups(
      self.intensities.crystal_symmetry(),
      max_delta=self.lattice_symmetry_max_delta,
      bravais_types_only=False)
    self.cb_op_min_best = self.subgroups.result_groups[0]['cb_op_inp_best']
    self.lattice_group = self.subgroups.result_groups[0]['best_subsym'].space_group()
    self.lattice_group = self.lattice_group.change_basis(self.cb_op_min_best.inverse())
    self.patterson_group = self.lattice_group.build_derived_patterson_group()

    sel = self.patterson_group.epsilon(self.intensities.indices()) == 1
    self.intensities = self.intensities.select(sel).set_info(
      self.intensities.info())
    self.dataset_ids = self.dataset_ids.select(sel)

    # Correct SDs by "typical" SD factors
    self.correct_sigmas(sd_fac=2.0, sd_b=0.0, sd_add=0.03)

    if normalisation is not None:
      if normalisation == 'kernel':
        normalise = self.kernel_normalisation
      elif normalisation == 'quasi':
        normalise = self.quasi_normalisation
      elif normalisation == 'ml_iso':
        normalise = self.ml_iso_normalisation
      elif normalisation == 'ml_aniso':
        normalise = self.ml_aniso_normalisation

      for i in range(int(flex.max(self.dataset_ids)+1)):
        logger.info('Normalising intensities for dataset %i' % (i+1))
        sel = self.dataset_ids == i
        intensities = self.intensities.select(self.dataset_ids == i)
        if i == 0:
          normalised_intensities = normalise(intensities)
        else:
          normalised_intensities = normalised_intensities.concatenate(
            normalise(intensities))
      self.intensities = normalised_intensities.set_info(
        self.intensities.info()).set_observation_type_xray_intensity()

    if d_min is not None or d_min is libtbx.Auto:
      self.resolution_filter(d_min, min_i_mean_over_sigma_mean, min_cc_half)

  def correct_sigmas(self, sd_fac, sd_b, sd_add):
    # sd' = SDfac * Sqrt(sd^2 + SdB * I + (SDadd * I)^2)
    sigmas = sd_fac * flex.sqrt(
      flex.pow2(self.intensities.sigmas() + (sd_b * self.intensities.data()) + flex.pow2(sd_add * self.intensities.data())))
    variance = flex.pow2(self.intensities.sigmas())
    si2 = flex.pow2(sd_add * self.intensities.data())
    ssc = variance + sd_b * self.intensities.data() + si2
    MINVARINFRAC = 0.1
    ssc.set_selected(ssc < MINVARINFRAC * variance, MINVARINFRAC * variance)
    sd = sd_fac * flex.sqrt(ssc)
    self.intensities = self.intensities.customized_copy(
      sigmas=sd).set_info(self.intensities.info())

  @staticmethod
  def kernel_normalisation(intensities):
    normalisation = absolute_scaling.kernel_normalisation(
      intensities, auto_kernel=True)
    return normalisation.normalised_miller.deep_copy().set_info(
      intensities.info())

  @staticmethod
  def quasi_normalisation(intensities):
    #handle negative reflections to minimise effect on mean I values.
    intensities.data().set_selected(intensities.data() < 0.0, 0.0)

    #set up binning objects
    if intensities.size() > 20000:
      n_refl_shells = 20
    elif intensities.size() > 15000:
      n_refl_shells = 15
    else:
      n_refl_shells = 10
    d_star_sq = intensities.d_star_sq().data()
    step = (flex.max(d_star_sq) - flex.min(d_star_sq) + 1e-8) / n_refl_shells
    binner = intensities.setup_binner_d_star_sq_step(d_star_sq_step=step)

    normalisations = intensities.intensity_quasi_normalisations()
    return intensities.customized_copy(
      data=(intensities.data()/normalisations.data()),
      sigmas=(intensities.sigmas()/normalisations.data()))

  @staticmethod
  def ml_aniso_normalisation(intensities):
    return symmetry_base._ml_normalisation(intensities, aniso=True)

  @staticmethod
  def ml_iso_normalisation(intensities):
    return symmetry_base._ml_normalisation(intensities, aniso=False)

  @staticmethod
  def _ml_normalisation(intensities, aniso):
    # estimate number of residues per unit cell
    mr = matthews.matthews_rupp(intensities.crystal_symmetry())
    n_residues = mr.n_residues

    # estimate B-factor and scale factors for normalisation
    if aniso:
      normalisation = absolute_scaling.ml_aniso_absolute_scaling(
        intensities, n_residues=n_residues)
      u_star = normalisation.u_star
    else:
      normalisation = absolute_scaling.ml_iso_absolute_scaling(
        intensities, n_residues=n_residues)
      u_star = adptbx.b_as_u(
        adptbx.u_iso_as_u_star(
          intensities.unit_cell(), normalisation.b_wilson))

    # record output in log file
    if aniso:
      b_cart = normalisation.b_cart
      logger.info('ML estimate of overall B_cart value:')
      logger.info('''\
  %5.2f, %5.2f, %5.2f
  %12.2f, %5.2f
  %19.2f''' % (b_cart[0], b_cart[3], b_cart[4],
                          b_cart[1], b_cart[5],
                                     b_cart[2]))
    else:
      logger.info('ML estimate of overall B value:')
      logger.info('   %5.2f A**2' %normalisation.b_wilson)
    logger.info('ML estimate of  -log of scale factor:')
    logger.info('  %5.2f' %(normalisation.p_scale))

    s = StringIO()
    mr.show(out=s)
    normalisation.show(out=s)
    logger.debug(s.getvalue())

    # apply scales
    return intensities.customized_copy(
      data=scaling.ml_normalise_aniso(
        intensities.indices(), intensities.data(),
        normalisation.p_scale, intensities.unit_cell(),
        u_star),
      sigmas=scaling.ml_normalise_aniso(
        intensities.indices(), intensities.sigmas(),
        normalisation.p_scale, intensities.unit_cell(),
        u_star))

  def resolution_filter(self, d_min, min_i_mean_over_sigma_mean, min_cc_half):
    if d_min is libtbx.Auto and (
        min_i_mean_over_sigma_mean is not None or min_cc_half is not None):
      from dials.util import Resolutionizer
      rparams = Resolutionizer.phil_defaults.extract().resolutionizer
      rparams.nbins = 20
      resolutionizer = Resolutionizer.resolutionizer(self.intensities, rparams)
      d_min_isigi = 0
      d_min_cc_half = 0
      if min_i_mean_over_sigma_mean is not None:
        d_min_isigi = resolutionizer.resolution_i_mean_over_sigma_mean(min_i_mean_over_sigma_mean)
        logger.info('Resolution estimate from <I>/<sigI> > %.1f : %.2f' % (
          min_i_mean_over_sigma_mean, d_min_isigi))
      if min_cc_half is not None:
        d_min_cc_half = resolutionizer.resolution_cc_half(min_cc_half)
        logger.info('Resolution estimate from CC1/2 > %.2f: %.2f' % (
          min_cc_half, d_min_cc_half))
      d_min = min(d_min_isigi, d_min_cc_half)
      logger.info('High resolution limit set to: %.2f' % d_min)
    if d_min is not None:
      sel = self.intensities.resolution_filter_selection(d_min=d_min)
      self.intensities = self.intensities.select(sel).set_info(
        self.intensities.info())
      self.dataset_ids = self.dataset_ids.select(sel)
      logger.info('Selecting %i reflections with d > %.2f' % (self.intensities.size(), d_min))

