#include "pome/transport/losses.h"

#include "pome/cgs.h"

namespace pome {

double magneticEnergyDensity(double B) { return B * B / 8. / M_PI; }

double KnFactor(const double& T, const double& gamma_e) {
  double value = 45. * pow2(cgs::electronMassC2) / 64. / pow2(M_PI) / pow2(cgs::kBoltzmann * T);
  return value / (value + pow2(gamma_e));
}

double Losses::get(double E, double magneticField, double tauAdiabatic) {
  double gamma_e = E / cgs::electronMassC2;
  double energyDensity = 0.26 * cgs::eV / cgs::cm3 * KnFactor(2.7 * cgs::K, gamma_e);  // CMB
  // energyDensity += 0.30 * cgs::eV / cgs::cm3 * KnFactor(20 * cgs::K, gamma_e);         // IR TODO add this
  // energyDensity += 0.30 * cgs::eV / cgs::cm3 * KnFactor(5000 * cgs::K, gamma_e);       // star
  // energyDensity += 0.10 * cgs::eV / cgs::cm3 * KnFactor(20000 * cgs::K, gamma_e);      // UV
  energyDensity += magneticEnergyDensity(magneticField);
  double dEdt = 4. / 3. * cgs::cLight * cgs::sigmaTh * energyDensity * pow2(gamma_e);
  dEdt += E / tauAdiabatic;
  return dEdt;
}

}  // namespace pome