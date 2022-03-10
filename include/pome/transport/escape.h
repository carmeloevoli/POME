#ifndef POME_ESCAPE_H
#define POME_ESCAPE_H

#include <cmath>
#include <utility>

namespace pome {

struct EscapeParams {
  double magneticField;
  double rPwn;
};

class Escape {
 public:
  Escape() {}
  Escape(EscapeParams p) : m_params(p) {}
  virtual ~Escape() = default;
  double get(double E, double t = 0) const;

 protected:
 protected:
  EscapeParams m_params;
};

}  // namespace pome

#endif  // POME_ESCAPE_H