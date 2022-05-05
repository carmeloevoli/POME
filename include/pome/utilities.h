#ifndef POME_UTILITIES_H
#define POME_UTILITIES_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <cmath>
#include <fstream>
#include <functional>
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

// Integration
template <typename T>
T QAGIntegration(std::function<T(T)> f, T start, T stop, int LIMIT, double rel_error = 1e-4) {
  double a = static_cast<double>(start);
  double b = static_cast<double>(stop);
  double abs_error = 0.0;  // disabled
  int key = GSL_INTEG_GAUSS31;
  double result;
  double error;

  gsl_function F;
  F.function = [](double x, void *vf) -> double {
    auto &func = *static_cast<std::function<double(double)> *>(vf);
    return func(x);
  };
  F.params = &f;

  gsl_integration_workspace *workspace_ptr = gsl_integration_workspace_alloc(LIMIT);
  gsl_integration_qag(&F, a, b, abs_error, rel_error, LIMIT, key, workspace_ptr, &result, &error);
  gsl_integration_workspace_free(workspace_ptr);

  return T(result);
}

}  // namespace utils
}  // namespace pome

#endif
