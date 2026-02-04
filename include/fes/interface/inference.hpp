// Copyright (c) 2026 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
/// @file include/fes/interface/inference.hpp
/// @brief Inference interface.
#pragma once

#include <vector>

#include "fes/constituent.hpp"

namespace fes {

/// Forward declaration of WaveTableInterface
/// @tparam Derived The derived type of the wave table interface.
template <typename Derived>
class WaveTableInterface;

/// @brief Inference interface.
/// @tparam Derived The derived type of the inference.
/// @details This class uses the Curiously Recurring Template Pattern (CRTP)
/// to allow static polymorphism for different inference implementations.
/// The derived class must implement the `apply_impl` method to define the
/// specific inference logic.
template <typename Derived>
class Inference {
 public:
  /// @brief Destructor.
  virtual ~Inference() = default;

  /// @brief Apply inference to compute minor constituents.
  /// @tparam WaveTableInterfaceDerived The derived type of the wave table
  /// interface.
  /// @param[in,out] wave_table The wave table to process.
  /// @param[in] lat Latitude in degrees.
  template <typename WaveTableInterfaceDerived>
  auto apply(WaveTableInterface<WaveTableInterfaceDerived>& wave_table,
             const double lat) -> void {
    static_cast<Derived*>(this)->apply_impl(wave_table, lat);
  }

  /// @brief Returns the list of the tidal constituents inferred by the model.
  /// @return A vector of constituent identifiers.
  virtual auto inferred_constituents() const -> std::vector<ConstituentId> = 0;
};

}  // namespace fes
