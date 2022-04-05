#include <memory>

#include "pome.h"

using namespace pome;

int main() {
  try {
    Timer timer;
    ModelState state;

    auto Q = std::unique_ptr<Spectrum>(new Spectrum(state));
    auto energyAxis = utils::LogAxis(1. * cgs::GeV, 1. * cgs::PeV, 600);
    utils::OutputFile out("source_spectrum.txt");
    out << std::scientific;
    for (auto E : energyAxis) {
      out << E / cgs::GeV << "\t";
      out << E * E * Q->get(E) / (cgs::erg / cgs::sec) << "\t";
      out << "\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}