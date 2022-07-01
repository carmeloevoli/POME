#include "pome/common.h"

#include <algorithm>
#include <cmath>

#include "pome/cgs.h"

namespace pome {

double freeExpansionPwnRadius(const double& t, const double& M_ej, const double& rho_0, const double& E_SN,
                              const double& L_pulsar, const double& tau_0) {
  const auto R_ch = std::pow(M_ej / rho_0, 1. / 3.);                                             // Eq. 5
  const auto t_ch = std::pow(E_SN, -0.5) * std::pow(M_ej, 5. / 6.) * std::pow(rho_0, -1. / 3.);  // Eq. 6
  const auto L_ch = E_SN / t_ch;                                                                 // Eq. 7
  const auto tau_star = tau_0 / t_ch;                                                            // Eq. 8
  const auto L_star = L_pulsar / L_ch;                                                           // Eq. 9
  const auto t_star = t / t_ch;
  const auto x = t_star / tau_star;

  auto cnst = 1.911 * std::pow(L_star * tau_star, 1. / 5.) * tau_star;
  auto fx = std::pow(1. + 0.965 * std::pow(x, 0.719), 1.390) / std::pow(1. + 1.157 * std::pow(x, -0.730), 1.645);
  return cnst * fx * R_ch;
}

double freeExpansionPwnVelocity(const double& t, const double& M_ej, const double& rho_0, const double& E_SN,
                                const double& L_pulsar, const double& tau_0) {
  const auto R_ch = std::pow(M_ej / rho_0, 1. / 3.);
  const auto t_ch = std::pow(E_SN, -0.5) * std::pow(M_ej, 5. / 6.) * std::pow(rho_0, -1. / 3.);
  const auto L_ch = E_SN / t_ch;
  const auto L_star = L_pulsar / L_ch;
  const auto tau_star = tau_0 / t_ch;
  const auto t_star = t / t_ch;
  const auto x = t_star / tau_star;

  auto cnst = 1.911 * std::pow(L_star * tau_star, 1. / 5.);
  auto dfdx = 1.38938 * std::pow(1. + 0.965 * std::pow(x, 0.719), 1.39) / std::pow(x, 1.73) /
              std::pow(1. + 1.157 / std::pow(x, 0.73), 2.645);
  dfdx += 0.964431 * std::pow(1. + 0.965 * std::pow(x, 0.719), 0.39) / std::pow(x, 0.281) /
          std::pow(1. + 1.157 / std::pow(x, 0.73), 1.645);
  return cnst * dfdx * (R_ch / t_ch);
}

double potentialDropEnergy(const double& Omega_0, const double& BSurface, const double& pulsarRadius) {
  const double factor = 0.5 * BSurface * pow3(pulsarRadius) / pow2(cgs::cLight) * cgs::eV / cgs::volt;
  return factor * pow2(Omega_0);
}

double decayTime(const double& Omega_0, const double& BSurface, const double& pulsarRadius, const double& I) {
  return 3. * pow3(cgs::cLight) * I / pow2(BSurface) / std::pow(pulsarRadius, 6.) / pow2(Omega_0);
}

double computeMagneticField(double magneticFieldEnergy, double radius) {
  auto volume = (4. / 3.) * M_PI * pow3(radius);
  auto w_B = magneticFieldEnergy / volume;
  return std::sqrt(8. * M_PI * w_B);
}

}  // namespace pome