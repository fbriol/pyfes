// Copyright (c) 2026 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
/// @file include/fes/map.hpp
/// @brief Map using binary search for small set.
#pragma once

#include <array>
#include <cstddef>

namespace fes {

/// @brief Map using binary search for small set.
/// @tparam KeyT Type of the keys.
/// @tparam ValueT Type of the values.
/// @tparam N Number of elements in the map.
template <typename KeyT, typename ValueT, size_t N>
class Map {
  /// Array of key-value pairs.
  std::array<std::pair<KeyT, ValueT>, N> data_{};

 public:
  /// @brief Constructor.
  /// @param[in] data Array of key-value pairs.
  explicit constexpr Map(std::array<std::pair<KeyT, ValueT>, N>&& data) noexcept
      : data_(std::move(data)) {}

  /// @brief Default constructor.
  constexpr Map() noexcept = default;

  /// @brief Gets the value associated with a key.
  /// @param[in] key The key.
  /// @return Pointer to the associated value, or nullptr if the key is not
  /// present.
  constexpr auto get(KeyT key) const noexcept -> const ValueT* {
    for (const auto& pair : data_) {
      if (pair.first == key) {
        return &pair.second;
      }
    }
    return nullptr;
  }

  /// @brief Gets the value associated with a key.
  /// @param[in] key The key.
  /// @return Reference to the associated value.
  constexpr auto operator[](KeyT key) const noexcept -> const ValueT& {
    return *get(key);
  }

  /// @brief Checks if a key is present in the map.
  /// @param[in] key The key.
  /// @return true if the key is present, false otherwise.
  constexpr auto contains(KeyT key) const noexcept -> bool {
    return get(key) != nullptr;
  }

  /// @brief Gets the number of elements in the map.
  /// @return Number of elements.
  constexpr auto size() const noexcept -> size_t { return N; }

  /// @brief Gets an iterator to the beginning.
  /// @return Iterator to the beginning.
  constexpr auto begin() const noexcept { return data_.begin(); }

  /// @brief Gets an iterator to the end.
  /// @return Iterator to the end.
  constexpr auto end() const noexcept { return data_.end(); }
};

}  // namespace fes
