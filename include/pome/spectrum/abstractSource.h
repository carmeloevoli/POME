#ifndef POME_ABSTRACTSOURCE_H
#define POME_ABSTRACTSOURCE_H

namespace pome {

class AbstractSource {
 public:
  AbstractSource() {}
  virtual double get(double E, double t = 0) const = 0;
};

}  // namespace pome

#endif  // POME_ABSTRACTSOURCE_H