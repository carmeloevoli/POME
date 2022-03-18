#include "pome/spectrum/brokenPowerLaw.h"

#include "pome/common.h"

namespace pome {

double BrokenPowerLaw::get(double E, double t) const {
  const auto x = E / m_state.E_b;
  const auto alpha_1 = m_state.alpha_1;
  const auto alpha_2 = m_state.alpha_2;
  double Q = sourceNormalization(t);
  Q *= (x < 1) ? std::pow(x, -alpha_1) : std::pow(x, -alpha_2);
  Q *= std::exp(-E / E_c(t));
  return Q;
}

double BrokenPowerLaw::sourceNormalization(const double &t) const {
  const double slope = -(m_state.n + 1) / (m_state.n - 1);
  const double L_pairs = (1. - m_state.eta_B) * m_state.L_0;
  const double tau_0 = decayTime(2. * M_PI / m_state.P_0, cgs::pulsarBSurface, cgs::pulsarRadius, cgs::pulsarInertia);
  double value = L_pairs / m_state.E_b / m_state.E_b;
  value *= I(t);
  value *= std::pow(1. + t / tau_0, slope);
  return value;
}

double BrokenPowerLaw::E_c(const double &t) const {
  const double Omega_0 = 2. * M_PI / m_state.P_0;
  const double tau_0 = decayTime(Omega_0, cgs::pulsarBSurface, cgs::pulsarRadius, cgs::pulsarInertia);
  const double V_0 = potentialDropEnergy(Omega_0, cgs::pulsarBSurface, cgs::pulsarRadius);
  double value = V_0 / (1. + t / tau_0);
  return value;
}

double BrokenPowerLaw::I(const double &t) const {
  const double alpha_1 = m_state.alpha_1;
  const double alpha_2 = m_state.alpha_2;
  double value = 1. / (2. - alpha_1);
  value += 1. / (alpha_2 - 2.) * (1. - std::pow(E_c(t) / m_state.E_b, 2. - alpha_2));
  return 1. / value;
}

}  // namespace pome