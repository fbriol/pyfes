// Copyright (c) 2026 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#include "fes/perth/wave.hpp"

#include <pybind11/complex.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "fes/perth/constituent.hpp"
#include "fes/perth/wave_table.hpp"

namespace py = pybind11;

namespace fes {
namespace perth {

inline auto init_wave(py::module& m) -> void {
  py::class_<Wave, WaveInterface, std::unique_ptr<Wave>>(
      m, "Wave",
      R"__doc__(
A tidal wave using Doodson's notation system.

Perth waves are described by a 7-element Doodson number vector that directly
encodes the tidal argument in terms of astronomical variables.
)__doc__")
      .def_property_readonly("name", &Wave::name,
                             "The name of the tidal wave.")
      .def_property_readonly("ident", &Wave::ident,
                             "The constituent identifier.")
      .def_property_readonly("type", &Wave::type,
                             "The type of tidal wave.")
      .def(
          "frequency",
          [](const Wave& self, const FrequencyUnit unit) -> double {
            return unit == FrequencyUnit::kRadianPerHour
                       ? self.frequency()
                       : detail::math::degrees(self.frequency());
          },
          py::arg("unit") = FrequencyUnit::kRadianPerHour,
          R"__doc__(
Get the frequency of the tidal wave.

Args:
  unit: The frequency unit. Default is radians per hour.

Returns:
  The frequency of the tidal wave.
)__doc__")
      .def_property_readonly("period", &Wave::period,
                             "The period of the wave in hours.")
      .def("doodson_numbers", &Wave::doodson_numbers,
           "Get the Doodson number of the wave.")
      .def("xdo_numerical", &Wave::xdo_numerical,
           "Get the XDO numerical representation of the wave.")
      .def("xdo_alphabetical", &Wave::xdo_alphabetical,
           "Get the XDO alphabetical representation of the wave.")
      .def("__repr__", [](const Wave& self) -> std::string {
        return std::string("<perth.Wave '") + self.name() + "'>";
      });
}

inline auto init_wave_table(py::module& m) -> void {
  py::class_<WaveTable, WaveTableInterface, std::unique_ptr<WaveTable>>(
      m, "WaveTable",
      R"__doc__(
Table of tidal constituents using Doodson's notation system.

This table manages 80 tidal constituents used by GOT/Perth models. Each
constituent is represented as a Perth wave with its Doodson number encoding.
)__doc__")
      .def(py::init<>(),
           "Default constructor with all known Perth constituents.")
      .def(py::init<const std::vector<std::string>&>(), py::arg("names"),
           R"__doc__(
Constructor with a list of constituent names.

Args:
  names: List of constituent names to include in the table.
)__doc__")
      .def("__repr__", [](const WaveTable& self) -> std::string {
        return "<perth.WaveTable with " + std::to_string(self.size()) +
               " constituents>";
      });
}

inline auto init_constituents(py::module& m) -> void {
  auto constituents_mod =
      m.def_submodule("constituents", "Perth constituent utilities.");

  constituents_mod.def("parse", &constituents::parse, py::arg("name"),
                       R"__doc__(
Parse a tidal constituent name to its identifier.

Parsing is case insensitive. So ``Mm``, ``MM`` and ``mm`` are equivalent.

Args:
  name: The constituent name.

Returns:
  The constituent identifier.

Raises:
  ValueError: If the constituent name is not recognized by the Perth model.
)__doc__");

  constituents_mod.def(
      "name",
      [](const ConstituentId ident) -> std::string {
        return constituents::name(ident);
      },
      py::arg("ident"),
      R"__doc__(
Get the name of a tidal constituent from its identifier.

Args:
  ident: The constituent identifier.

Returns:
  The constituent name.
)__doc__");

  constituents_mod.def(
      "known",
      []() -> std::vector<std::string> {
        auto names = constituents::known();
        return {names.begin(), names.end()};
      },
      R"__doc__(
Get all tidal constituent names handled by the Perth model.

Returns:
  List of all 80 Perth constituent names.
)__doc__");
}

}  // namespace perth
}  // namespace fes

void init_perth(py::module& m) {
  auto perth_mod = m.def_submodule("perth", "Perth tidal model components.");
  fes::perth::init_wave(perth_mod);
  fes::perth::init_wave_table(perth_mod);
  fes::perth::init_constituents(perth_mod);
}
