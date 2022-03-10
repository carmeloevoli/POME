#ifndef POME_BROKENPOWERLAW_H
#define POME_BROKENPOWERLAW_H

#include <cmath>
#include <utility>

#include "pome/source/abstractSource.h"

namespace pome {

struct BrokenPowerLawParams {
  std::pair<double, double> alpha = {0, 0};
  double E_b = 0;
  double V_0 = 0;
  double tau_0 = 0;
  double L_B = 0;
  double n = 0;
};

class BrokenPowerLaw : public AbstractSource {
 public:
  BrokenPowerLaw() {}
  BrokenPowerLaw(BrokenPowerLawParams p) : m_params(p) {}
  virtual ~BrokenPowerLaw() = default;
  double get(double E, double t = 0) const override;

 protected:
  double sourceNormalization(const double &t) const;
  double E_c(const double &t) const;
  double I(const double &t) const;

 protected:
  BrokenPowerLawParams m_params;
};

}  // namespace pome

#endif  // POME_BROKENPOWERLAW_H