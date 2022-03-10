#include "pome/modelState.h"

namespace pome {

double potentialDropEnergy(double Omega_0, double BSurface, double pulsarRadius) {
  const double factor = 0.5 * BSurface * pow3(pulsarRadius) / pow2(cgs::cLight) * cgs::eV / cgs::volt;
  return factor * pow2(Omega_0);
}

double decayTime(double Omega_0, double BSurface, double pulsarRadius, double I) {
  return 3. * pow3(cgs::cLight) * I / pow2(BSurface) / std::pow(pulsarRadius, 6.) / pow2(Omega_0);
}

void ModelState::dump() {}

}  // namespace pome