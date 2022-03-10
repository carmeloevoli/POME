#ifndef POME_MODELSTATE_H
#define POME_MODELSTATE_H

#include <iostream>
#include <utility>

#include "cgs.h"

namespace pome {

double potentialDropEnergy(double Omega_0, double BSurface, double pulsarRadius);
double decayTime(double Omega_0, double BSurface, double pulsarRadius, double I);

class ModelState {
 private:
  const double k_BSurface = std::pow(10., 12.65) * cgs::gauss;
  const double k_ESN = 1e51 * cgs::erg;
  const double k_pulsarMass = 1.4 * cgs::sunMass;
  const double k_pulsarRadius = 10. * cgs::km;
  const double k_I = (2. / 5.) * k_pulsarMass * pow2(k_pulsarRadius);

 public:
  ModelState() {
    // Input
    m_alpha_1 = 1.5;
    m_alpha_2 = 2.5;
    m_E_b = 7e5 * cgs::electronMassC2;
    m_eta_B = 0.03;
    m_n = 2.509;
    m_P_0 = 19. * cgs::msec;
    m_L_0 = 3.1e39 * cgs::erg / cgs::sec;
    // Derived parameters
    m_L_B = (1. - m_eta_B) * m_L_0;
    m_Omega_0 = 2. * M_PI / m_P_0;
    m_tau_0 = decayTime(m_Omega_0, k_BSurface, k_pulsarRadius, k_I);
    m_V_0 = potentialDropEnergy(m_Omega_0, k_BSurface, k_pulsarRadius);
  }

  virtual ~ModelState() = default;
  void dump();

 public:
  const double &alpha_1 = m_alpha_1;
  const double &alpha_2 = m_alpha_2;
  const double &E_b = m_E_b;
  const double &eta_B = m_eta_B;
  const double &L_0 = m_L_0;
  const double &L_B = m_L_B;
  const double &n = m_n;
  const double &Omega_0 = m_Omega_0;
  const double &tau_0 = m_tau_0;
  const double &V_0 = m_V_0;

 private:
  double m_alpha_1;
  double m_alpha_2;
  double m_E_b;
  double m_eta_B;
  double m_L_0;
  double m_L_B;
  double m_n;
  double m_Omega_0;
  double m_P_0;
  double m_tau_0;
  double m_V_0;
};

}  // namespace pome

#endif  // POME_MODELSTATE_H