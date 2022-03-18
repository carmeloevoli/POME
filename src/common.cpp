#include "pome/common.h"

#include <algorithm>
#include <cmath>

#include "pome/cgs.h"

namespace pome {

double freeExpansionPwnRadius(const double& t, const double& M_ej, const double& rho_0, const double& E_SN,
                              const double& L_pulsar, const double& tau_0) {
  const auto R_ch = std::pow(M_ej / rho_0, 1. / 3.);
  const auto t_ch = std::pow(E_SN, -0.5) * std::pow(M_ej, 5. / 6.) * std::pow(rho_0, -1. / 3.);
  const auto L_ch = E_SN / t_ch;
  const auto L_star = L_pulsar / L_ch;
  const auto tau_star = tau_0 / t_ch;
  const auto t_star = t / t_ch;
  const auto x = t_star / tau_star;

  auto value = 1.911 * std::pow(L_star * tau_star, 1. / 5.) * tau_star;
  value *= std::pow(1. + 0.965 * std::pow(x, 0.719), 1.390) / std::pow(1. + 1.157 * std::pow(x, -0.730), 1.645);
  return std::max(value * R_ch, 1e-3 * cgs::pc);
}

double potentialDropEnergy(const double& Omega_0, const double& BSurface, const double& pulsarRadius) {
  const double factor = 0.5 * BSurface * pow3(pulsarRadius) / pow2(cgs::cLight) * cgs::eV / cgs::volt;
  return factor * pow2(Omega_0);
}

double decayTime(const double& Omega_0, const double& BSurface, const double& pulsarRadius, const double& I) {
  return 3. * pow3(cgs::cLight) * I / pow2(BSurface) / std::pow(pulsarRadius, 6.) / pow2(Omega_0);
}

double luminosityDecay(const double& t, const double& tau_0, const double& n) {
  const double slope = -(n + 1) / (n - 1);
  return std::pow(1. + t / tau_0, slope);
}

}  // namespace pome