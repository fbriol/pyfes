#include "fes/inference.hpp"

#include <gtest/gtest.h>

#include "fes/darwin//wave_table.hpp"
#include "fes/interface/wave_table.hpp"
#include "fes/perth/wave_table.hpp"

namespace fes {

static void admittance(const double r, const ConstituentId ident,
                       WaveTableInterface& table) {
  EXPECT_NEAR(table[ident]->tide().real(), r, 1e-6);
  EXPECT_NEAR(table[ident]->tide().imag(), r, 1e-6);

//   table[ident]->admittance(false);
//   table[ident]->tide({1, 1});
//   table.admittance();
//   EXPECT_EQ(table[ident]->tide(), std::complex<double>(1, 1));

//   table[ident]->admittance(true);
//   table[ident]->tide({1, 1});
}

TEST(InferenceTest, SplineInference) {
  auto wt = fes::darwin::WaveTable();
  for (auto& item : wt) {
    item.value()->set_tide({1, 1});
  }

  auto inference = inference_factory(wt, InferenceType::kSpline);
  inference->apply(wt, 45.0);
  admittance(0.2378, k2Q1, wt);
  admittance(0.2706, kSigma1, wt);
  admittance(0.1688, kRho1, wt);
  admittance(0.0241, kM11, wt);
  admittance(0.0671, kM12, wt);
  admittance(0.0124, kChi1, wt);
  admittance(0.0201, kPi1, wt);
  admittance(0.0137, kPhi1, wt);
  admittance(0.009, kTheta1, wt);
  admittance(0.0447, kJ1, wt);
  admittance(0.0182, kOO1, wt);
  admittance(0.0796782, kEta2, wt);
  admittance(0.374697218357, kMu2, wt);
  admittance(0.157529811402, kNu2, wt);
  admittance(0.010949128375, kLambda2, wt);
  admittance(0.053354227357, kL2, wt);
  admittance(0.16871051505, kT2, wt);
  admittance(0.2387, k2N2, wt);
  admittance(0.094151295, kEps2, wt);
  // disable 2N2 from the admittance interpolation, and check that it falls back
  // to the default value
  wt[k2N2]->set_is_modeled(true);
  wt[k2N2]->set_tide({1, 1});
  inference->apply(wt, 45.0);
  admittance(0.499810, kEps2, wt);
}

}  // namespace fes