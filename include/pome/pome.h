#ifndef POME_POME_H_
#define POME_POME_H_

#include <memory>

#include "modelState.h"

namespace POME {

class Pome {
 public:
  Pome(ModelState state) : m_state(state) {}
  virtual ~Pome() = default;

  void buildSource();
  void dumpSource() const;
  void evolve() const;

 protected:
  ModelState m_state;
  std::unique_ptr<AbstractSource> m_Q;
};

}  // namespace POME

#endif  // POME_POME_H_