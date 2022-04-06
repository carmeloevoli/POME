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
    p.evolve(1e-3 * cgs::year, 0.1 * cgs::kyr, 10.1 * cgs::kyr);

  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}