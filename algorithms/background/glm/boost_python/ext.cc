/*
 * ext.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/background/glm/robust_poisson_mean.h>
#include <dials/algorithms/background/glm/creator.h>

namespace dials { namespace algorithms { namespace background {
  namespace boost_python {

  using namespace boost::python;

  BOOST_PYTHON_MODULE(dials_algorithms_background_glm_ext)
  {
    class_<RobustPoissonMean>("RobustPoissonMean", no_init)
      .def(init<const af::const_ref<double>&,
                double,
                double,
                double,
                std::size_t>((
                  arg("Y"),
                  arg("mean0"),
                  arg("c")=1.345,
                  arg("tolerance")=1e-3,
                  arg("max_iter")=100)))
      .def("mean", &RobustPoissonMean::mean)
      .def("niter", &RobustPoissonMean::niter)
      .def("error", &RobustPoissonMean::error)
      .def("converged", &RobustPoissonMean::converged)
      ;


    class_<Creator> creator("Creator", no_init);
    creator
      .def(init<
          Creator::Model,
          double,
          std::size_t>((
              arg("model"),
              arg("tuning_constant"),
              arg("max_iter"))))
      .def("__call__", &Creator::shoebox)
      .def("__call__", &Creator::volume)
      ;

    scope in_creator = creator;

    enum_<Creator::Model>("model")
      .value("constant2d", Creator::Constant2d)
      .value("constant3d", Creator::Constant3d)
      .value("loglinear2d", Creator::LogLinear2d)
      .value("loglinear3d", Creator::LogLinear3d)
      ;
  }

}}}} // namespace = dials::algorithms::background::boost_python
