#include <memory>

#include "pome.h"

using namespace pome;

int main() {
  try {
    Timer timer;
    ModelState state;

    auto tauEscape = std::unique_ptr<Escape>(new Escape(state));
    auto losses = std::unique_ptr<Losses>(new Losses(state));
    auto energyAxis = utils::LogAxis(1. * cgs::GeV, 1. * cgs::PeV, 600);

    double B = 10. * cgs::muG;
    double rPWN = cgs::pc;
    double tauAdiabatic = 1e3 * cgs::kyr;

    utils::OutputFile out("timescales.txt");
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << tauEscape->get(E, B, rPWN) / (cgs::kyr) << "\t";
      out << E / losses->get(E, B, tauAdiabatic) / (cgs::kyr) << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}