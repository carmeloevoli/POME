#ifndef POME_UTILITIES_H
#define POME_UTILITIES_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace pome {
namespace utils {

// Axis
std::vector<double> LinAxis(const double &min, const double &max, const size_t &size);
std::vector<double> LogAxis(const double &min, const double &max, const size_t &size);

// Output File
class OutputFile {
  std::string filename;
  std::ofstream out;

 public:
  explicit OutputFile(const std::string &name) : filename(name), out("output/" + name) {}
  ~OutputFile() { std::cout << "created output file " << filename << "\n"; }

  template <typename T>
  OutputFile &operator<<(const T &value) {
    out << value;
    return *this;
  }
};

}  // namespace utils
}  // namespace pome

#endif
