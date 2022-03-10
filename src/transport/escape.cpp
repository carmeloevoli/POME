#include "pome/transport/escape.h"

#include <iostream>

#include "pome/cgs.h"

namespace pome {

double larmorRadius(double pc, double B) { return pc / cgs::elementaryCharge / B; }

double Escape::get(double E, double t) const {
  const auto pc = std::sqrt(pow2(E) - pow2(cgs::electronMassC2));
  const auto rL = larmorRadius(pc, m_params.magneticField);
  const auto DBohm = 1. / 3. * cgs::cLight * rL;
  return pow2(m_params.rPwn) / DBohm;
}

}  // namespace pome