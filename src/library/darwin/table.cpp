// Copyright (c) 2026 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#include "fes/darwin/table.hpp"

#include <memory>

#include "fes/constituent.hpp"
#include "fes/darwin/constituent.hpp"
#include "fes/darwin/wave.hpp"
#include "fes/interface.hpp"

namespace fes {
namespace darwin {

auto WaveTable::wave_factory_impl(ConstituentId ident)
    -> std::unique_ptr<WaveInterface> {
  switch (ident) {
    case kO1:
      return std::unique_ptr<WaveInterface>(new wave::O1());
    case kP1:
      return std::unique_ptr<WaveInterface>(new wave::P1());
    case kK1:
      return std::unique_ptr<WaveInterface>(new wave::K1());
    case k2N2:
      return std::unique_ptr<WaveInterface>(new wave::_2N2());
    case kMu2:
      return std::unique_ptr<WaveInterface>(new wave::Mu2());
    case kN2:
      return std::unique_ptr<WaveInterface>(new wave::N2());
    case kNu2:
      return std::unique_ptr<WaveInterface>(new wave::Nu2());
    case kM2:
      return std::unique_ptr<WaveInterface>(new wave::M2());
    case kL2:
      return std::unique_ptr<WaveInterface>(new wave::L2());
    case kT2:
      return std::unique_ptr<WaveInterface>(new wave::T2());
    case kS2:
      return std::unique_ptr<WaveInterface>(new wave::S2());
    case kK2:
      return std::unique_ptr<WaveInterface>(new wave::K2());
    case kM4:
      return std::unique_ptr<WaveInterface>(new wave::M4());
    case kS1:
      return std::unique_ptr<WaveInterface>(new wave::S1());
    case kQ1:
      return std::unique_ptr<WaveInterface>(new wave::Q1());
    case kMm:
      return std::unique_ptr<WaveInterface>(new wave::Mm());
    case kMf:
      return std::unique_ptr<WaveInterface>(new wave::Mf());
    case kMtm:
      return std::unique_ptr<WaveInterface>(new wave::Mtm());
    case kMSqm:
      return std::unique_ptr<WaveInterface>(new wave::MSqm());
    case kEps2:
      return std::unique_ptr<WaveInterface>(new wave::Eps2());
    case kLambda2:
      return std::unique_ptr<WaveInterface>(new wave::Lambda2());
    case kEta2:
      return std::unique_ptr<WaveInterface>(new wave::Eta2());
    case k2Q1:
      return std::unique_ptr<WaveInterface>(new wave::_2Q1());
    case kSigma1:
      return std::unique_ptr<WaveInterface>(new wave::Sigma1());
    case kRho1:
      return std::unique_ptr<WaveInterface>(new wave::Rho1());
    case kM1:
      return std::unique_ptr<WaveInterface>(new wave::M1());
    case kM11:
      return std::unique_ptr<WaveInterface>(new wave::M11());
    case kM12:
      return std::unique_ptr<WaveInterface>(new wave::M12());
    case kM13:
      return std::unique_ptr<WaveInterface>(new wave::M13());
    case kChi1:
      return std::unique_ptr<WaveInterface>(new wave::Chi1());
    case kPi1:
      return std::unique_ptr<WaveInterface>(new wave::Pi1());
    case kPhi1:
      return std::unique_ptr<WaveInterface>(new wave::Phi1());
    case kTheta1:
      return std::unique_ptr<WaveInterface>(new wave::Theta1());
    case kJ1:
      return std::unique_ptr<WaveInterface>(new wave::J1());
    case kOO1:
      return std::unique_ptr<WaveInterface>(new wave::OO1());
    case kM3:
      return std::unique_ptr<WaveInterface>(new wave::M3());
    case kM6:
      return std::unique_ptr<WaveInterface>(new wave::M6());
    case kMN4:
      return std::unique_ptr<WaveInterface>(new wave::MN4());
    case kMS4:
      return std::unique_ptr<WaveInterface>(new wave::MS4());
    case kN4:
      return std::unique_ptr<WaveInterface>(new wave::N4());
    case kR2:
      return std::unique_ptr<WaveInterface>(new wave::R2());
    case kR4:
      return std::unique_ptr<WaveInterface>(new wave::R4());
    case kS4:
      return std::unique_ptr<WaveInterface>(new wave::S4());
    case kMNS2:
      return std::unique_ptr<WaveInterface>(new wave::MNS2());
    case kMK4:
      return std::unique_ptr<WaveInterface>(new wave::MK4());
    case kSN4:
      return std::unique_ptr<WaveInterface>(new wave::SN4());
    case kSK4:
      return std::unique_ptr<WaveInterface>(new wave::SK4());
    case k2MN6:
      return std::unique_ptr<WaveInterface>(new wave::_2MN6());
    case k2MS6:
      return std::unique_ptr<WaveInterface>(new wave::_2MS6());
    case k2MK6:
      return std::unique_ptr<WaveInterface>(new wave::_2MK6());
    case kMSN6:
      return std::unique_ptr<WaveInterface>(new wave::MSN6());
    case k2SM6:
      return std::unique_ptr<WaveInterface>(new wave::_2SM6());
    case kMSK6:
      return std::unique_ptr<WaveInterface>(new wave::MSK6());
    case kMP1:
      return std::unique_ptr<WaveInterface>(new wave::MP1());
    case k2SM2:
      return std::unique_ptr<WaveInterface>(new wave::_2SM2());
    case kPsi1:
      return std::unique_ptr<WaveInterface>(new wave::Psi1());
    case k2MS2:
      return std::unique_ptr<WaveInterface>(new wave::_2MS2());
    case kMKS2:
      return std::unique_ptr<WaveInterface>(new wave::MKS2());
    case k2MN2:
      return std::unique_ptr<WaveInterface>(new wave::_2MN2());
    case kMSN2:
      return std::unique_ptr<WaveInterface>(new wave::MSN2());
    case kMO3:
      return std::unique_ptr<WaveInterface>(new wave::MO3());
    case k2MK3:
      return std::unique_ptr<WaveInterface>(new wave::_2MK3());
    case kMK3:
      return std::unique_ptr<WaveInterface>(new wave::MK3());
    case kS6:
      return std::unique_ptr<WaveInterface>(new wave::S6());
    case kM8:
      return std::unique_ptr<WaveInterface>(new wave::M8());
    case kMSf:
      return std::unique_ptr<WaveInterface>(new wave::MSf());
    case kSsa:
      return std::unique_ptr<WaveInterface>(new wave::Ssa());
    case kSa:
      return std::unique_ptr<WaveInterface>(new wave::Sa());
    case kA5:
      return std::unique_ptr<WaveInterface>(new wave::A5());
    case kSa1:
      return std::unique_ptr<WaveInterface>(new wave::Sa1());
    case kSta:
      return std::unique_ptr<WaveInterface>(new wave::Sta());
    case kMm2:
      return std::unique_ptr<WaveInterface>(new wave::Mm2());
    case kMm1:
      return std::unique_ptr<WaveInterface>(new wave::Mm1());
    case kMf1:
      return std::unique_ptr<WaveInterface>(new wave::Mf1());
    case kMf2:
      return std::unique_ptr<WaveInterface>(new wave::Mf2());
    case kM0:
      return std::unique_ptr<WaveInterface>(new wave::M0());
    case kL2P:
      return std::unique_ptr<WaveInterface>(new wave::L2P());
    case kN2P:
      return std::unique_ptr<WaveInterface>(new wave::N2P());
    case kMSK2:
      return std::unique_ptr<WaveInterface>(new wave::MSK2());
    case kSKM2:
      return std::unique_ptr<WaveInterface>(new wave::SKM2());
    case kOQ2:
      return std::unique_ptr<WaveInterface>(new wave::OQ2());
    case k3MS4:
      return std::unique_ptr<WaveInterface>(new wave::_3MS4());
    case kMNu4:
      return std::unique_ptr<WaveInterface>(new wave::MNu4());
    case k2MSN4:
      return std::unique_ptr<WaveInterface>(new wave::_2MSN4());
    case k2NS2:
      return std::unique_ptr<WaveInterface>(new wave::_2NS2());
    case kMNuS2:
      return std::unique_ptr<WaveInterface>(new wave::MNuS2());
    case k2MK2:
      return std::unique_ptr<WaveInterface>(new wave::_2MK2());
    case kNKM2:
      return std::unique_ptr<WaveInterface>(new wave::NKM2());
    case kML4:
      return std::unique_ptr<WaveInterface>(new wave::ML4());
    case kSO1:
      return std::unique_ptr<WaveInterface>(new wave::SO1());
    case kSO3:
      return std::unique_ptr<WaveInterface>(new wave::SO3());
    case kNK4:
      return std::unique_ptr<WaveInterface>(new wave::NK4());
    case kMNK6:
      return std::unique_ptr<WaveInterface>(new wave::MNK6());
    case k2NM6:
      return std::unique_ptr<WaveInterface>(new wave::_2NM6());
    case k3MS8:
      return std::unique_ptr<WaveInterface>(new wave::_3MS8());
    case kSK3:
      return std::unique_ptr<WaveInterface>(new wave::SK3());
    case k2MNS4:
      return std::unique_ptr<WaveInterface>(new wave::_2MNS4());
    case k2SMu2:
      return std::unique_ptr<WaveInterface>(new wave::_2SMu2());
    case k2MP5:
      return std::unique_ptr<WaveInterface>(new wave::_2MP5());
    default:
      throw std::invalid_argument("wave identifier not recognized: " +
                                  std::to_string(ident));
  }
}

// Builds the list of constituent identifiers from the optional list of names
static auto build_constituent_ids(
    const boost::optional<std::vector<std::string>>& waves)
    -> std::vector<ConstituentId> {
  auto to_ids = [](const auto& names) {
    std::vector<ConstituentId> result;
    result.reserve(names.size());
    std::transform(names.begin(), names.end(), std::back_inserter(result),
                   darwin::constituents::parse);
    return result;
  };

  return waves ? to_ids(*waves) : to_ids(darwin::constituents::known());
}

WaveTable::WaveTable(const boost::optional<std::vector<std::string>>& waves)
    : WaveTableInterface(build_constituent_ids(waves)) {}

auto WaveTable::compute_nodal_corrections(const angle::Astronomic& angles)
    -> void {
  for (auto& item : *this) {
    item.value()->compute_nodal_corrections(angles);
  }
}

void WaveTable::process_inference() {
  // Arrays who contains the spline coefficients needed to compute MU2, NU2,
  // L2, T2 and Lambda2 by admittance.
  constexpr auto mu2 =
      std::array<double, 3>{0.069439968323, 0.351535557706, -0.046278307672};
  constexpr auto nu2 =
      std::array<double, 3>{-0.006104695053, 0.156878802427, 0.006755704028};
  constexpr auto l2 =
      std::array<double, 3>{0.077137765667, -0.051653455134, 0.027869916824};
  constexpr auto t2 =
      std::array<double, 3>{0.180480173707, -0.020101177502, 0.008331518844};
  constexpr auto lda2 =
      std::array<double, 3>{0.016503557465, -0.013307812292, 0.007753383202};

  // infer additional constituents by admittance DIURNALS (from Richard Ray
  // perth2 program)

  // from Q1 and O1 (0-1)

  auto* x = (*this)[kQ1].get();
  auto* y = (*this)[kO1].get();
  auto* z = (*this)[kK1].get();

  auto set_tide = [](WaveInterface* wave, const std::complex<double>& value) {
    if (wave->is_inferred()) {
      wave->set_tide(value);
    }
  };

  // 2Q1
  set_tide((*this)[k2Q1].get(), 0.263 * x->tide() - 0.0252 * y->tide());

  // Sigma1
  set_tide((*this)[kSigma1].get(), 0.297 * x->tide() - 0.0264 * y->tide());

  // rho1
  set_tide((*this)[kRho1].get(), 0.164 * x->tide() + 0.0048 * y->tide());

  // from O1 and K1  (1-2)

  // M11
  set_tide((*this)[kM11].get(), 0.0140 * y->tide() + 0.0101 * z->tide());

  // M12
  set_tide((*this)[kM12].get(), 0.0389 * y->tide() + 0.0282 * z->tide());

  // CHI1
  set_tide((*this)[kChi1].get(), 0.0064 * y->tide() + 0.0060 * z->tide());

  // pi1
  set_tide((*this)[kPi1].get(), 0.0030 * y->tide() + 0.0171 * z->tide());

  // phi1
  set_tide((*this)[kPhi1].get(), -0.0015 * y->tide() + 0.0152 * z->tide());

  // theta1
  set_tide((*this)[kTheta1].get(), -0.0065 * y->tide() + 0.0155 * z->tide());

  // J1
  set_tide((*this)[kJ1].get(), -0.0389 * y->tide() + 0.0836 * z->tide());

  // OO1
  set_tide((*this)[kOO1].get(), -0.0431 * y->tide() + 0.0613 * z->tide());

  // infer additional constituents by admittance SEMI-DIURNALS
  // (from Richard Ray perth3 program)

  // from M2 - N2
  x = (*this)[kN2].get();
  y = (*this)[kM2].get();

  // 2N2
  set_tide((*this)[k2N2].get(), 0.264 * x->tide() - 0.0253 * y->tide());

  // SEMI-DIURNAL (from Grenoble to take advantage of 2N2)

  // from 2N2 -N2 (3-4)
  x = (*this)[k2N2].get();
  y = (*this)[kN2].get();

  // eps2
  set_tide((*this)[kEps2].get(), 0.53285 * x->tide() - 0.03304 * y->tide());
  // from M2 - K2 [5-6]
  x = (*this)[kN2].get();
  y = (*this)[kM2].get();
  z = (*this)[kK2].get();

  // eta2
  set_tide((*this)[kEta2].get(),
           -0.0034925 * y->tide() + 0.0831707 * z->tide());

  // from N2 -M2- K2 by spline admittances [see GRL 18[5]:845-848,1991]

  // mu2
  set_tide((*this)[kMu2].get(),
           mu2[0] * z->tide() + mu2[1] * x->tide() + mu2[2] * y->tide());
  // nu2
  set_tide((*this)[kNu2].get(),
           nu2[0] * z->tide() + nu2[1] * x->tide() + nu2[2] * y->tide());
  // lambda2
  set_tide((*this)[kLambda2].get(),
           lda2[0] * z->tide() + lda2[1] * x->tide() + lda2[2] * y->tide());
  // L2
  set_tide((*this)[kL2].get(),
           l2[0] * z->tide() + l2[1] * x->tide() + l2[2] * y->tide());

  // T2
  set_tide((*this)[kT2].get(),
           t2[0] * z->tide() + t2[1] * x->tide() + t2[2] * y->tide());
}

}  // namespace darwin
}  // namespace fes