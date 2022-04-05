#ifndef POME_LOSSES_H
#define POME_LOSSES_H

#include <cmath>
#include <utility>

#include "pome/modelState.h"

namespace pome {

class EnergyLosses {
 public:
  EnergyLosses() {}
  EnergyLosses(ModelState m) : m_state(m) {}
  virtual ~EnergyLosses() = default;
  double bLosses(double E, double magneticField) const;
  double bAdiabatic(double E, double rPwn, double vPwn) const;
  double get(double E, double magneticField, double rPwn, double vPwn) const {
    return bLosses(E, magneticField) + bAdiabatic(E, rPwn, vPwn);
  }

 protected:
  ModelState m_state;
};

}  // namespace pome

#endif  // POME_LOSSES_H