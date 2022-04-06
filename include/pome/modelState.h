#ifndef POME_MODELSTATE_H
#define POME_MODELSTATE_H

#include <iostream>
#include <utility>

#include "cgs.h"

namespace pome {

class ModelState {
 public:
  ModelState() {
    // Input
    m_alpha_1 = 1.25;
    m_alpha_2 = 2.49;
    m_E_b = 5e5 * cgs::electronMassC2;
    m_eta_B = 0.03;
    m_n = 3.0;
    m_P_0 = 0.073 * cgs::sec;
    m_L_0 = 2.26e36 * cgs::erg / cgs::sec;
    m_E_SN = 1e51 * cgs::erg;
    m_M_ejecta = 18.6 * cgs::sunMass;
    m_rho_0 = 0.68 * cgs::protonMass / cgs::cm3;
    // Energy Axis
    m_eSize = 6 * 32;
    m_minEnergy = cgs::GeV;
    m_maxEnergy = cgs::PeV;
    dump();
  }

  virtual ~ModelState() = default;
  void dump();

 public:
  const double &alpha_1 = m_alpha_1;
  const double &alpha_2 = m_alpha_2;
  const double &E_b = m_E_b;
  const double &eta_B = m_eta_B;
  const double &n = m_n;
  const double &P_0 = m_P_0;
  const double &L_0 = m_L_0;
  const double &E_SN = m_E_SN;
  const double &M_ejecta = m_M_ejecta;
  const double &rho_0 = m_rho_0;
  const size_t &eSize = m_eSize;
  const double &minEnergy = m_minEnergy;
  const double &maxEnergy = m_maxEnergy;

 private:
  double m_alpha_1;
  double m_alpha_2;
  double m_E_b;
  double m_eta_B;
  double m_n;
  double m_P_0;
  double m_L_0;
  double m_E_SN;
  double m_M_ejecta;
  double m_rho_0;
  size_t m_eSize;
  double m_minEnergy;
  double m_maxEnergy;
};

}  // namespace pome

#endif  // POME_MODELSTATE_H