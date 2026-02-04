// Copyright (c) 2026 CNES
//
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.
#include "fes/constituent.hpp"

#include <array>
#include <stdexcept>

#include "fes/detail/string.hpp"

namespace fes {
namespace constituents {

constexpr auto kConstituentNames =
    std::array<const char* const, kNumConstituents>{
        "2mk2",  "2mk3",  "2mk6",  "2mn2",   "2mn6", "2mns4",   "2mp5",
        "2ms2",  "2ms6",  "2msn4", "2n2",    "2nm6", "2ns2",    "2q1",
        "2sm2",  "2sm6",  "2smu2", "3ms4",   "3ms8", "a5",      "alpha2",
        "beta1", "beta2", "chi1",  "delta2", "eps2", "eta2",    "gamma2",
        "j1",    "k1",    "k2",    "l2",     "l2p",  "lambda2", "m0",
        "m1",    "m11",   "m12",   "m13",    "m2",   "m3",      "m4",
        "m6",    "m8",    "mf",    "mf1",    "mf2",  "mk3",     "mk4",
        "mks2",  "ml4",   "mm",    "mm1",    "mm2",  "mn4",     "mnk6",
        "mns2",  "mnu4",  "mnus2", "mo3",    "mp1",  "mqm",     "ms4",
        "msf",   "msk2",  "msk6",  "msm",    "msn2", "msn6",    "msqm",
        "mstm",  "mtm",   "mu2",   "n2",     "n2p",  "n4",      "nk4",
        "nkm2",  "node",  "nu2",   "o1",     "oo1",  "oq2",     "p1",
        "phi1",  "pi1",   "psi1",  "q1",     "r2",   "r4",      "rho1",
        "s1",    "s2",    "s4",    "s6",     "sa",   "sa1",     "sigma1",
        "sk3",   "sk4",   "skm2",  "sn4",    "so1",  "so3",     "ssa",
        "sta",   "t2",    "tau1",  "theta1", "ups1",
    };

auto parse(const std::string& constituent_name) -> ConstituentId {
  for (size_t ix = 0; ix < kNumConstituents; ++ix) {
    if (detail::iequals(constituent_name, std::string{kConstituentNames[ix]})) {
      return static_cast<ConstituentId>(ix);
    }
  }
  throw std::invalid_argument("unknown constituent name: " + constituent_name);
}

auto name(ConstituentId constituent) -> const char* {
  const auto ix = static_cast<std::size_t>(constituent);
  if (ix >= kNumConstituents) {
    throw std::invalid_argument("unknown constituent: " +
                                std::to_string(static_cast<int>(constituent)));
  }
  return kConstituentNames[ix];
}

}  // namespace constituents
}  // namespace fes
