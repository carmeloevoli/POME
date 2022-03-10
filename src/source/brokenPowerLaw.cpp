#include "pome/source/brokenPowerLaw.h"

namespace pome {

double BrokenPowerLaw::get(double E, double t) const {
  const auto x = E / m_params.E_b;
  const auto alpha_1 = m_params.alpha.first;
  const auto alpha_2 = m_params.alpha.second;
  double Q = sourceNormalization(t);
  Q *= (x < 1) ? std::pow(x, -alpha_1) : std::pow(x, -alpha_2);
  Q *= std::exp(-E / E_c(t));
  return Q;
}

double BrokenPowerLaw::sourceNormalization(const double &t) const {
  const double slope = -(m_params.n + 1) / (m_params.n - 1);
  double value = m_params.L_B / m_params.E_b / m_params.E_b;
  value *= I(t);
  value *= std::pow(1. + t / m_params.tau_0, slope);
  return value;
}

double BrokenPowerLaw::E_c(const double &t) const {
  double value = m_params.V_0;
  value /= (1. + t / m_params.tau_0);
  return value;
}

double BrokenPowerLaw::I(const double &t) const {
  const double alpha_1 = m_params.alpha.first;
  const double alpha_2 = m_params.alpha.second;
  double value = 1. / (2. - alpha_1);
  value += 1. / (alpha_2 - 2.) * (1. - std::pow(E_c(t) / m_params.E_b, 2. - alpha_2));
  return 1. / value;
}

}  // namespace pome