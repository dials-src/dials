/*
 * fitting.cc
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
#include <boost/python/iterator.hpp>
#include <dials/algorithms/integration/profile/fitting.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_fitting()
  {
    class_<ProfileModel>("ProfileModel", no_init)
      .def(init<const flex_double&,
                const flex_double&,
                const flex_double&>((
        arg("profile"),
        arg("contents"),
        arg("background"))))
      .def("__call__", &ProfileModel::operator(), (arg("I")))
      .def("variance", &ProfileModel::variance, (arg("I")));
      
    class_<ProfileFitting>("ProfileFitting", no_init)
      .def(init<const flex_double&,
                const flex_double&,
                const flex_double&,
                int,
                std::size_t>((
        arg("profile"),
        arg("contents"),
        arg("background"),
        arg("bits") = 16,
        arg("max_iter") = 50)))
      .def("intensity", &ProfileFitting::intensity)
      .def("variance", &ProfileFitting::variance);
      
    class_<ProfileFitting2>("ProfileFitting2", no_init)
      .def(init<const flex_double&,
                const flex_double&,
                const flex_double&,
                double,
                std::size_t>((
        arg("profile"),
        arg("contents"),
        arg("background"),
        arg("bits") = 1e-3,
        arg("max_iter") = 10)))
      .def("intensity", &ProfileFitting2::intensity)
      .def("variance", &ProfileFitting2::variance)
      .def("niter", &ProfileFitting2::niter)
      .def("error", &ProfileFitting2::error);
  }

}}} // namespace = dials::algorithms::boost_python
