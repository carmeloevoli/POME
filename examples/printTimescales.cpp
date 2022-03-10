#include <memory>

#include "pome.h"

using namespace pome;

int main() {
  try {
    Timer timer;
    ModelState state;

    auto tauEscape = std::unique_ptr<Escape>(new Escape({10. * cgs::muG, 1. * cgs::pc}));
    auto losses = std::unique_ptr<Losses>(new Losses({10. * cgs::muG}));
    auto energyAxis = utils::LogAxis(1. * cgs::GeV, 1. * cgs::PeV, 600);

    utils::OutputFile out("timescales.txt");
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << tauEscape->get(E) / (cgs::kyr) << "\t";
      out << E / losses->get(E) / (cgs::kyr) << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}