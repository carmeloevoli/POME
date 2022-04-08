#include "pome.h"

using namespace pome;

int main() {
  try {
    Timer timer;
    ModelState state;

    Pome p(state);
    p.buildEnergyAxis();
    p.buildSource();
    p.buildLosses();
    p.evolve(1e-4 * cgs::year, 1.0 * cgs::kyr, 20.1 * cgs::kyr);

  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}