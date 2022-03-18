#ifndef POME_ESCAPE_H
#define POME_ESCAPE_H

#include <cmath>
#include <utility>

#include "pome/modelState.h"

namespace pome {

class Escape {
 public:
  Escape() {}
  Escape(ModelState m) : m_state(m) {}
  virtual ~Escape() = default;
  double get(double E, double magneticField, double rPwn) const;

 protected:
 protected:
  ModelState m_state;
};

}  // namespace pome

#endif  // POME_ESCAPE_H