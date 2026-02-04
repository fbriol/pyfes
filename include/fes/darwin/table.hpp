// Copyright (c) 2026 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
/// @file include/fes/darwin/table.hpp
/// @brief Table of tidal constituents handled by FES models.
#pragma once

#include <memory>
#include <string>
#include <vector>
#include <boost/optional.hpp>

#include "fes/interface.hpp"

namespace fes {
namespace darwin {

/// Properties of tide waves handled by FES models.
class WaveTable : public WaveTableInterface<WaveTable> {
 public:
  /// @brief Build a table from a list of constituent names
  explicit WaveTable(
      const boost::optional<std::vector<std::string>>& waves = {});

  /// @brief Generates a wave properties from its identifier
  /// @param[in] ident Wave identifier
  /// @return Wave properties
  static auto wave_factory_impl(ConstituentId ident) 
        -> std::unique_ptr<WaveInterface>;

  /// @brief Computes the nodal corrections for all constituents in the table.
  ///
  /// @param[in] angles Astronomic angles used to compute the nodal corrections.
  auto compute_nodal_corrections(const angle::Astronomic& angles) -> void final;

  /// @brief Computes the missing constituents by inference.
  auto process_inference() -> void final;

  /// @brief Clones the wave table.
  /// @return A unique pointer to the cloned wave table.
  auto clone() const -> std::unique_ptr<WaveTableInterface> final {
    return std::make_unique<WaveTable>(*this);
  }
};

}  // namespace darwin
}  // namespace fes