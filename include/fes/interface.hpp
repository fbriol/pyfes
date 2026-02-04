// Copyright (c) 2026 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
/// @file include/fes/wave.hpp
/// @brief Tidal wave interface.
#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "fes/angle/astronomic.hpp"
#include "fes/constituent.hpp"
#include "fes/detail/math.hpp"
#include "fes/detail/parallel_for.hpp"
#include "fes/enum_map.hpp"
#include "fes/types.hpp"

namespace fes {

/// @brief Frequency unit enumeration.
enum FrequencyUnit : uint8_t {
  kRadian,  ///< Radian per hour
  kDegree   ///< Degree per hour
};

/// @brief Possible type of tidal wave.
enum WaveType : uint8_t {
  kLongPeriod,  //!< Long period tidal waves
  kShortPeriod  //!< Short period tidal waves
};

/// @brief Tidal wave interface.
class WaveInterface {
 public:
  /// @brief Constructor.
  WaveInterface(const ConstituentId ident, WaveType type) noexcept
      : id_(ident), type_(type) {}

  /// @brief Destructor.
  virtual ~WaveInterface() = default;

  /// @brief Clones the wave.
  /// @return A unique pointer to the cloned wave.
  virtual auto clone() const -> std::unique_ptr<WaveInterface> = 0;

  /// @brief Indicates whether the wave is inferred.
  /// @return `true` if the wave is inferred, `false` otherwise.
  constexpr auto is_inferred() const noexcept -> bool { return is_inferred_; }

  /// @brief Sets whether the wave is inferred.
  /// @param inferred `true` if the wave is inferred, `false` otherwise.
  constexpr void set_is_inferred(bool is_inferred) noexcept {
    is_inferred_ = is_inferred;
  }

  /// @brief Indicates whether the tide is provided by the model.
  /// @return `true` if the tide is provided by the model, `false` otherwise.
  constexpr auto is_modeled() const noexcept -> bool { return is_modeled_; }

  /// @brief Sets whether the tide is provided by the model.
  /// @param is_modeled `true` if the tide is provided by the model, `false`
  /// otherwise.
  constexpr void set_is_modeled(bool is_modeled) noexcept {
    is_modeled_ = is_modeled;
  }

  /// @brief Gets the constituent identifier.
  /// @return The constituent identifier.
  constexpr auto ident() const noexcept -> ConstituentId { return id_; }

  /// @brief Gets the name of the tidal wave.
  /// @return The name of the tidal wave.
  inline auto name() const -> const char* { return constituents::name(id_); }

  /// @brief Gets the type of tidal wave.
  /// @return The type of tidal wave.
  constexpr auto type() const noexcept -> WaveType { return type_; }

  /// @brief Gets the frequency.
  /// @param unit The frequency unit.
  /// @return The frequency in radians per hour or degrees per hour.
  virtual auto frequency(const FrequencyUnit& unit) const noexcept
      -> double = 0;

  /// Gets the period of the wave (hours)
  constexpr auto period() const noexcept -> double {
    return detail::math::two_pi<double>() / frequency(kRadian);
  }

  /// @brief Gets the tide value.
  /// @return The tide value.
  constexpr auto tide() const noexcept -> const Complex& { return tide_; }

  /// @brief Sets the tide value.
  /// @param tide The tide value.
  constexpr void set_tide(const Complex& tide) noexcept { tide_ = tide; }

  /// @brief Gets the Greenwich argument, in radians.
  /// @return The Greenwich argument.
  constexpr auto v() const noexcept -> double { return v_; }

  /// @brief Gets the nodal correction for phase, in radians.
  /// @return The nodal correction for phase.
  constexpr auto u() const noexcept -> double { return u_; }

  /// @brief Gets the nodal correction for amplitude, in radians.
  /// @return The nodal correction for amplitude.
  constexpr auto f() const noexcept -> double { return f_; }

  /// @brief Gets the sum of the Greenwich argument and the nodal phase
  /// correction, in radians.
  /// @return The sum of the Greenwich argument and the nodal phase
  /// correction, normalized to the interval [0, 2π).
  inline auto vu() const noexcept -> double {
    return detail::math::normalize_angle(v_ + u_, 0.0,
                                         detail::math::two_pi<double>());
  }

  /// @brief Computes the nodal corrections for the wave.
  /// @param[in] angles Astronomic angles used to compute the nodal corrections.
  virtual auto compute_nodal_corrections(const angle::Astronomic& angles)
      -> void = 0;

  /// Gets the XDO numerical representation of the wave
  virtual auto xdo_numerical() const -> std::string = 0;

  /// Gets the XDO alphabetical representation of the wave
  virtual auto xdo_alphabetical() const -> std::string = 0;

  /// Gets the Doodson number of the wave
  /// @note The 7th number follows the convention established in Doodson &
  /// Warburg's 1941 book. This number can be 0, 1, 2, or -1, representing
  /// multiples of 90 degrees added to the tidal argument when using cosine
  /// functions
  virtual auto doodson_numbers() const -> Vector7b = 0;

 protected:
  Complex tide_{};      ///< Tide value
  double v_{};          ///< Greenwich argument
  double f_{};          ///< Nodal correction for amplitude
  double u_{};          ///< Nodal correction for phase
  ConstituentId id_;    ///< Constituent identifier
  WaveType type_;       ///< Type of tidal wave
  bool is_inferred_{};  ///< Indicates whether the wave is inferred
  bool is_modeled_{};   ///< Indicates whether the tide is provided by the model
};

/// @brief Map of constituents to wave interfaces.
using ConstituentMap =
    EnumMap<ConstituentId, std::unique_ptr<WaveInterface>, kKnownConstituents>;

/// @brief Tidal wave table interface.
template <typename Derived>
class WaveTableInterface {
 public:
  /// @brief Default constructor.
  /// @param constituents List of constituent identifiers to include in the
  /// table
  inline explicit WaveTableInterface(
      const std::vector<ConstituentId>& constituents) {
    for (const auto& ident : constituents) {
      map_.set(ident, WaveTableInterface::wave_factory(ident));
    }
  }

  /// @brief Copy constructor
  WaveTableInterface(const WaveTableInterface& other) {
    for (const auto& item : other.map_) {
      map_.set(item.key(), item.value()->clone());
    }
  }

  /// @brief Move constructor
  WaveTableInterface(WaveTableInterface&& other) noexcept = default;

  /// @brief Copy assignment operator
  auto operator=(const WaveTableInterface& other) -> WaveTableInterface& {
    if (this != &other) {
      map_.clear();
      for (const auto& item : other.map_) {
        map_.set(item.key(), item.value()->clone());
      }
    }
    return *this;
  }

  /// @brief Move assignment operator
  auto operator=(WaveTableInterface&& other) noexcept
      -> WaveTableInterface& = default;

  /// @brief Destructor.
  virtual ~WaveTableInterface() = default;

  /// @brief Clones the wave table.
  /// @return A unique pointer to the cloned wave table.
  virtual auto clone() const -> std::unique_ptr<WaveTableInterface> = 0;

  /// Create a wave properties from its identifier
  /// @param[in] ident Wave identifier
  /// @return Wave properties
  static auto wave_factory(ConstituentId ident)
      -> std::unique_ptr<WaveInterface> {
    return Derived::wave_factory_impl(ident);
  }

  /// Set the tide of a constituent
  /// @param[in] ident The constituent identifier
  /// @param[out] value The tide value
  void set_tide(ConstituentId ident, const Complex& value) {
    const auto* ptr = map_.get(ident);
    if (ptr == nullptr) {
      throw out_of_range(ident);
    }
    set_tide(static_cast<ConstituentId>(ident), value);
  }

  /// @brief Computes the nodal corrections for all constituents in the table.
  ///
  /// @param[in] angles Astronomic angles used to compute the nodal corrections.
  virtual auto compute_nodal_corrections(const angle::Astronomic& angles)
      -> void = 0;

  /// @brief Computes the missing constituents by inference.
  virtual auto process_inference() -> void = 0;

  /// @brief Get the wave corresponding at the given index.
  /// @param index The index of the constituent, must be in the range
  /// [0, size()).
  /// @return A constant reference to the unique pointer to the wave.
  inline auto operator[](const size_t index) const
      -> const std::unique_ptr<WaveInterface>& {
    return map_[index];
  }

  /// @brief Get the wave corresponding at the given index.
  /// @param index The index of the constituent, must be in the range
  /// [0, size()).
  /// @return A reference to the unique pointer to the wave.
  inline auto operator[](const size_t index)
      -> std::unique_ptr<WaveInterface>& {
    return map_[index];
  }

  /// @brief Get the wave corresponding to the given constituent identifier.
  /// @param ident The constituent identifier.
  /// @return A constant reference to the unique pointer to the wave.
  inline auto operator[](ConstituentId ident) const
      -> const std::unique_ptr<WaveInterface>& {
    const auto* ptr = map_.get(ident);
    if (ptr == nullptr) {
      throw out_of_range(ident);
    }
    return map_[ident];
  }

  /// @brief Get the wave corresponding to the given constituent identifier.
  /// @param ident The constituent identifier.
  /// @return A reference to the unique pointer to the wave.
  inline auto operator[](ConstituentId ident)
      -> std::unique_ptr<WaveInterface>& {
    const auto* ptr = map_.get(ident);
    if (ptr == nullptr) {
      throw out_of_range(ident);
    }
    return map_[ident];
  }

  /// @brief Returns an iterator to the beginning of the wave table
  /// @return An iterator to the beginning of the wave table
  inline auto begin() const noexcept -> ConstituentMap::const_iterator {
    return map_.begin();
  }

  /// @brief Returns an iterator to the end of the wave table
  /// @return An iterator to the end of the wave table
  inline auto end() const noexcept -> ConstituentMap::const_iterator {
    return map_.end();
  }

  /// @brief Returns an iterator to the beginning of the wave table
  /// @return An iterator to the beginning of the wave table
  inline auto begin() noexcept -> ConstituentMap::iterator {
    return map_.begin();
  }

  /// @brief Returns an iterator to the end of the wave table
  /// @return An iterator to the end of the wave table
  inline auto end() noexcept -> ConstituentMap::iterator { return map_.end(); }

  /// @brief Gets the list of constituent names in the table
  /// @return The list of constituent names in the table
  inline auto constituents() const -> std::vector<std::string> {
    auto names = std::vector<std::string>();
    names.reserve(map_.size());
    for (const auto& item : map_) {
      names.emplace_back(item.value()->name());
    }
    return names;
  }

  /// @brief Returns the size of the table
  inline auto size() const noexcept -> size_t { return map_.size(); }

  /// @brief Return the list of tidal waves such that their period is more than
  /// twice the duration of the time series analyzed
  ///
  /// @param[in] duration Duration of the time series analyzed in seconds
  /// @param[in] f Number of times the period of the wave is greater than
  /// the duration of the time series analyzed
  /// @return List of selected tidal waves.
  auto select_waves_for_analysis(double duration, double f = 2.0)
      -> std::vector<std::string>;

  /// @brief Calculate the tide of a given time series.
  ///
  /// @param[in] epoch Desired UTC time expressed in number of seconds elapsed
  /// since 1970-01-01T00:00:00.
  /// @param[in] wave Tidal wave properties computed by an harmonic analysis.
  /// @return the tide at the given time.
  /// @param[in] formulae The formulae used to compute the astronomical angles.
  /// @return the tide at the given time.
  auto tide_from_tide_series(
      const Eigen::Ref<const Eigen::VectorXd>& epoch,
      const Eigen::Ref<const Eigen::VectorXcd>& wave,
      const angle::Formulae& formulae = angle::Formulae::kSchuremanOrder3) const
      -> Eigen::VectorXd;

  /// @brief Calculate the tide for a given date from a grid describing the
  /// properties of tidal waves over an area.
  ///
  /// @param[in] epoch Desired UTC time expressed in number of seconds elapsed
  /// since 1970-01-01T00:00:00.
  /// @param[in] wave Tidal wave properties computed by an harmonic analysis.
  /// @param[in] formulae The formulae used to compute the astronomical angles.
  /// @param[in] num_threads Number of threads to use for the computation. If
  /// set to 0, the number of threads is automatically determined.
  auto tide_from_mapping(
      double epoch, const DynamicRef<const Eigen::MatrixXcd>& wave,
      const angle::Formulae& formulae = angle::Formulae::kSchuremanOrder3,
      size_t num_threads = 0) const -> Eigen::MatrixXd;

  /// @brief Compute nodal modulations for amplitude and phase.
  ///
  /// @param[in] epoch: Desired UTC time expressed in number of seconds elapsed
  /// since 1970-01-01T00:00:00.
  /// @param[in] formulae The formulae used to compute the astronomical angles.
  /// @return The nodal correction for amplitude, v greenwich argument) + u
  /// (nodal correction for phase).
  /// @throw std::invalid_argument if the size of the epoch vector is not
  /// equal to the size of the leap seconds vector.
  auto compute_nodal_modulations(const Eigen::Ref<const Eigen::VectorXd>& epoch,
                                 const angle::Formulae& formulae) const
      -> std::tuple<Eigen::MatrixXd, Eigen::MatrixXd>;

  /// @brief Check if a constituent is in the table
  /// @param[in] ident The constituent identifier
  /// @return true if the constituent is in the table
  inline auto contains(const ConstituentId ident) const noexcept -> bool {
    return map_.contains(ident);
  }

 private:
  ConstituentMap map_;  ///< Map of constituents

  /// @brief Generates an out_of_range exception for a missing constituent.
  /// @param ident The constituent identifier.
  /// @return The out_of_range exception.
  auto out_of_range(ConstituentId ident) const -> std::out_of_range {
    return std::out_of_range("Constituent ID#" +
                             std::to_string(static_cast<uint8_t>(ident)) +
                             " not found in the wave table.");
  }
};

// ============================================================================
// Implementation
// ============================================================================

template <typename Derived>
auto WaveTableInterface<Derived>::select_waves_for_analysis(
    const double duration, const double f) -> std::vector<std::string> {
  auto result = std::vector<std::string>();
  for (auto& item : *this) {
    auto& wave = item.value();
    if (wave->period() < f * (duration / 3600.0)) {
      result.emplace_back(wave->name());
    }
  }
  return result;
}

// ============================================================================

template <typename Derived>
auto WaveTableInterface<Derived>::tide_from_tide_series(
    const Eigen::Ref<const Eigen::VectorXd>& epoch,
    const Eigen::Ref<const Eigen::VectorXcd>& wave,
    const angle::Formulae& formulae) const -> Eigen::VectorXd {
  if (static_cast<size_t>(wave.rows()) != size()) {
    throw std::invalid_argument(
        "wave must contain as many elements as the number of waves in the "
        "table");
  }
  auto result = Eigen::VectorXd(epoch.rows());

  /// The object responsible for the calculation of astronomical angles.
  auto angles = angle::Astronomic(formulae);

  // The wave properties of the object must be immutable for the provided
  // instance.
  auto wt = clone();

  for (auto ix = 0; ix < epoch.rows(); ++ix) {
    double tide = 0;
    angles.update(epoch(ix));
    wt->compute_nodal_corrections(angles);

    for (size_t jx = 0; jx < wt->size(); ++jx) {
      const auto& item = (*wt)[jx];
      double phi = item->vu();

      tide += item->f() * (wave(jx).real() * std::cos(phi) +
                           wave(jx).imag() * std::sin(phi));
    }
    result(ix) = tide;
  }
  return result;
}

// ============================================================================

template <typename Derived>
auto WaveTableInterface<Derived>::tide_from_mapping(
    const double epoch, const DynamicRef<const Eigen::MatrixXcd>& wave,
    const angle::Formulae& formulae, const size_t num_threads) const
    -> Eigen::MatrixXd {
  if (static_cast<size_t>(wave.rows()) != size()) {
    throw std::invalid_argument(
        "wave must contain as many elements as the number of waves in the "
        "table");
  }
  auto result = Eigen::MatrixXd(wave.cols(), wave.rows());
  auto worker = [&](const int64_t start, const int64_t end) {
    // The wave properties of the object must be immutable for the provided
    // instance.
    auto wt = clone();
    wt->compute_nodal_corrections(angle::Astronomic(formulae, epoch));

    for (auto ix = start; ix < end; ++ix) {
      for (size_t jx = 0; jx < wt->size(); ++jx) {
        const auto& item = (*wt)[jx];
        double phi = item->vu();

        result(ix, jx) += item->f() * (wave(jx, ix).real() * std::cos(phi) +
                                       wave(jx, ix).imag() * std::sin(phi));
      }
    }
  };
  detail::parallel_for(worker, wave.cols(), num_threads);
  return result;
}

// ============================================================================

template <typename Derived>
auto WaveTableInterface<Derived>::compute_nodal_modulations(
    const Eigen::Ref<const Eigen::VectorXd>& epoch,
    const angle::Formulae& formulae) const
    -> std::tuple<Eigen::MatrixXd, Eigen::MatrixXd> {
  auto f = Eigen::MatrixXd(size(), epoch.size());
  auto vu = Eigen::MatrixXd(size(), epoch.size());

  /// The object responsible for the calculation of astronomical angles.
  auto angles = angle::Astronomic(formulae);

  // The wave properties of the object must be immutable for the provided
  // instance.
  auto wt = clone();

  for (auto ix = 0; ix < epoch.size(); ++ix) {
    angles.update(epoch(ix));
    wt->compute_nodal_corrections(angles);
    for (size_t jx = 0; jx < wt->size(); ++jx) {
      const auto& wave = (*wt)[jx];
      f(jx, ix) = wave->f();
      vu(jx, ix) = wave->vu();
    }
  }
  return std::make_tuple(f, vu);
}

}  // namespace fes
