#ifndef POME_SPECTRUM_H
#define POME_SPECTRUM_H

#include <cmath>
#include <utility>

#include "pome/modelState.h"

namespace pome {

class Spectrum {
 public:
  Spectrum() {}
  Spectrum(ModelState m);
  virtual ~Spectrum() = default;
  inline double getTau0() const { return m_tau_0; }
  inline double getV0() const { return m_V_0; }
  double getPairLuminosity(double t = 0) const;
  double getLuminosity(double t = 0) const;
  double get(double E, double t = 0) const;

 protected:
  double sourceNormalization(const double &t) const;
  double E_c(const double &t) const;
  double I(const double &t) const;

 protected:
  ModelState m_state;
  double m_tau_0 = 0;
  double m_V_0 = 0;
  double m_L_pairs = 0;
  double m_slope = 0;
};

}  // namespace pome

#endif  // POME_SPECTRUM_H