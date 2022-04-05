#include "pome/model/spectrum.h"

#include "pome/common.h"

namespace pome {

Spectrum::Spectrum(ModelState m) : m_state(m) {
  const auto Omega_0 = 2. * M_PI / m_state.P_0;
  m_tau_0 = decayTime(Omega_0, cgs::pulsarBSurface, cgs::pulsarRadius, cgs::pulsarInertia);
  m_V_0 = potentialDropEnergy(Omega_0, cgs::pulsarBSurface, cgs::pulsarRadius);
  m_L_pairs = (1. - m_state.eta_B) * m_state.L_0;
  m_slope = -(m_state.n + 1) / (m_state.n - 1);
  std::cout << "tau_0 : " << m_tau_0 / cgs::kyr << " kyr\n";
  std::cout << "V_0 : " << m_V_0 / cgs::TeV << " TeV\n";
  std::cout << "L_pairs : " << m_L_pairs / (cgs::erg / cgs::sec) << " erg/s\n";
  std::cout << "slope : " << m_slope << "\n";
}

double Spectrum::getPairLuminosity(double t) const { return m_L_pairs * std::pow(1. + t / m_tau_0, m_slope); }

double Spectrum::getLuminosity(double t) const { return m_state.L_0 * std::pow(1. + t / m_tau_0, m_slope); }

double Spectrum::get(double E, double t) const {
  const auto x = E / m_state.E_b;
  const auto alpha_1 = m_state.alpha_1;
  const auto alpha_2 = m_state.alpha_2;
  double Q = sourceNormalization(t);
  Q *= (x < 1) ? std::pow(x, -alpha_1) : std::pow(x, -alpha_2);
  Q *= std::exp(-E / E_c(t));
  return Q;
}

double Spectrum::sourceNormalization(const double &t) const {
  double value = getPairLuminosity(t) / m_state.E_b / m_state.E_b;
  value *= I(t);
  return value;
}

double Spectrum::E_c(const double &t) const {
  double value = m_V_0 / (1. + t / m_tau_0);
  return value;
}

double Spectrum::I(const double &t) const {
  const double alpha_1 = m_state.alpha_1;
  const double alpha_2 = m_state.alpha_2;
  double value = 1. / (2. - alpha_1);
  value += 1. / (alpha_2 - 2.) * (1. - std::pow(E_c(t) / m_state.E_b, 2. - alpha_2));
  return 1. / value;
}

}  // namespace pome