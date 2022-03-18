#include "pome/transport/escape.h"

#include <iostream>

#include "pome/cgs.h"

namespace pome {

double larmorRadius(double pc, double B) { return pc / cgs::elementaryCharge / B; }

double Escape::get(double E, double magneticField, double rPwn) const {
  const auto pc = std::sqrt(pow2(E) - pow2(cgs::electronMassC2));
  const auto rL = larmorRadius(pc, magneticField);
  const auto DBohm = 1. / 3. * cgs::cLight * rL;
  return pow2(rPwn) / DBohm;
}

}  // namespace pome