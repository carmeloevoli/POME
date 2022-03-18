#ifndef POME_BROKENPOWERLAW_H
#define POME_BROKENPOWERLAW_H

#include <cmath>
#include <utility>

#include "pome/modelState.h"
#include "pome/spectrum/abstractSource.h"

namespace pome {

class BrokenPowerLaw : public AbstractSource {
 public:
  BrokenPowerLaw() {}
  BrokenPowerLaw(ModelState m) : m_state(m) {}
  virtual ~BrokenPowerLaw() = default;
  double get(double E, double t = 0) const override;

 protected:
  double sourceNormalization(const double &t) const;
  double E_c(const double &t) const;
  double I(const double &t) const;

 protected:
  ModelState m_state;
};

}  // namespace pome

#endif  // POME_BROKENPOWERLAW_H