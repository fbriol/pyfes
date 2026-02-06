// Copyright (c) 2026 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#include "fes/harmonic_analysis.hpp"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace fes {

inline auto init_harmonic_analysis(py::module& m) -> void {
  m.def(
      "harmonic_analysis",
      [](const Eigen::Ref<const Eigen::VectorXd>& h,
         const DynamicRef<const Eigen::MatrixXd>& f,
         const DynamicRef<const Eigen::MatrixXd>& vu) -> Eigen::VectorXcd {
        py::gil_scoped_release release;
        return harmonic_analysis(h, f, vu);
      },
      py::arg("h"), py::arg("f"), py::arg("vu"),
      R"__doc__(
Harmonic analysis of the tide.

The harmonic analysis method consists in expressing the ocean tidal
variations as a sum of independent constituents accordingly to the tidal
potential spectrum. Then the sea surface elevation at a point
:math:`(x, y)` and time :math:`(t)` can be expressed as a linear sum as
follow:

.. math::

  S_{ap} = S_{0}(x, y) + \sum_{k=0}^n f_{k}(t)S_{k}(x, y)
    	imes \cos [\omega_{k}t + v_{k}(t) + u_{k}(t) - G_{k}(x,y)]

where:
  - :math:`(n)` is the number of constituents,
  - :math:`S_{0}(x, y)` is the mean sea level,
  - :math:`S_{k}(x, y)` is the amplitude of the constituent of index
    :math:`(k)`,
  - :math:`G_{k}(x, y)` is the phase lag relative to Greenwich time,
  - :math:`\omega_{k}` is the angular frequency of the constituent of index
    :math:`(k)`,
  - :math:`v_{k}` is the astronomical argument at time :math:`(t)`,
  - :math:`f_{k}(t)` is the nodal correction coefficient applied to
    the amplitude of the constituent of index :math:`(k)`,
  - :math:`u_{k}(t)` is the nodal correction coefficient applied to
    the phase of the constituent of index :math:`(k)`.

The a priori analysis spectrum includes the most important astronomical
constituents in the Darwin development, completed by Schureman in 1958,
and many non-linear waves. The definition of tidal constants and
astronomical arguments is taken from FES2014 tidal prediction software
and a complete definition of waves is also available in Schureman (1958).
This spectrum is the most commonly used for harmonic analysis due the
simplification given by the nodal correction concept (:math:`(f)` and
:math:`(u)` coefficients above) which allows dealing with slow motions of
the lunar ascending node and reducing the number of constituents in the
tidal spectrum. More details about this harmonic analysis method can be
found in Ponchaut et al. 1999.

Args:
  h: Sea level.
  f: Nodal correction coefficient applied to the amplitude of the
    constituents analyzed.
  vu: Astronomical argument at time :math:`(t)` + the nodal correction
    coefficient applied to the phase of the constituents analyzed.

Returns:
  The complex number representing the different reconstructed waves.
)__doc__");
}

}  // namespace fes

auto init_harmonic_analysis(py::module& m) -> void {
  fes::init_harmonic_analysis(m);
}