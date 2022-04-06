#ifndef POME_POME_H
#define POME_POME_H

#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

#include "pome/common.h"
#include "pome/model/escape.h"
#include "pome/model/losses.h"
#include "pome/model/spectrum.h"
#include "pome/modelState.h"
#include "pome/tridiag.h"
#include "pome/utilities.h"

namespace pome {

class Pome {
 public:
  Pome(ModelState state);
  virtual ~Pome() = default;

 public:
  void buildEnergyAxis();
  void buildSource();
  void buildLosses();

  void evolve(double dt, double t_dump, double t_max);

 protected:
  ModelState m_state;
  size_t m_eSize;
  std::unique_ptr<Spectrum> m_Q;
  double m_tau_0;
  std::unique_ptr<Escape> m_tauEscape;
  std::unique_ptr<EnergyLosses> m_losses;
  std::vector<double> m_energyAxis;
  std::vector<double> m_crDensity;
};

}  // namespace pome

#endif  // POME_POME_H