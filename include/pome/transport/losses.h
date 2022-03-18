#ifndef POME_LOSSES_H
#define POME_LOSSES_H

#include <cmath>
#include <utility>

#include "pome/modelState.h"

namespace pome {

class Losses {
 public:
  Losses() {}
  Losses(ModelState m) : m_state(m) {}
  virtual ~Losses() = default;
  double get(double E, double magneticField, double tauAdiabatic);

 protected:
 protected:
  ModelState m_state;
};

}  // namespace pome

#endif  // POME_LOSSES_H