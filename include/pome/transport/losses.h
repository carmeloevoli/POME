#ifndef POME_LOSSES_H
#define POME_LOSSES_H

#include <cmath>
#include <utility>

namespace pome {

struct LossesParams {
  double magneticField;
};

class Losses {
 public:
  Losses() {}
  Losses(LossesParams p) : m_params(p) {}
  virtual ~Losses() = default;
  double get(double E, double t = 0) const;

 protected:
 protected:
  LossesParams m_params;
};

}  // namespace pome

#endif  // POME_LOSSES_H