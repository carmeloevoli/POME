#include "pome/modelState.h"

#include "pome/common.h"

namespace pome {

void ModelState::dump() {
  std::cout << "alpha 1 : " << m_alpha_1 << "\n";
  std::cout << "alpha 2 : " << m_alpha_2 << "\n";
  std::cout << "E_b : " << m_E_b / cgs::TeV << " TeV\n";
  std::cout << "eta_B : " << m_eta_B << "\n";
  std::cout << "n : " << m_n << "\n";
  std::cout << "P_0 : " << m_P_0 / cgs::msec << " msec\n";
  std::cout << "L_0 : " << m_L_0 / (cgs::erg / cgs::sec) << " erg/s\n";
  // Derived parameters
  std::cout << "L_B : " << m_L_0 / (cgs::erg / cgs::sec) << " erg/s\n";
  std::cout << "Omega_0 : " << 2. * M_PI / m_P_0 * cgs::sec << " 1/s\n";
  std::cout << "tau_0 : "
            << decayTime(2. * M_PI / m_P_0, cgs::pulsarBSurface, cgs::pulsarRadius, cgs::pulsarInertia) / cgs::kyr
            << " kyr\n";
  std::cout << "V_0 : " << potentialDropEnergy(2. * M_PI / m_P_0, cgs::pulsarBSurface, cgs::pulsarRadius) / cgs::PeV
            << " PeV\n";
}

}  // namespace pome