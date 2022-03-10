#include <memory>

#include "pome.h"

using namespace pome;

int main() {
  try {
    Timer timer;
    ModelState state;

    std::pair<double, double> alpha = std::make_pair(state.alpha_1, state.alpha_2);
    BrokenPowerLawParams p = {alpha, state.E_b, state.V_0, state.tau_0, state.L_B, state.n};
    auto Q = std::unique_ptr<AbstractSource>(new BrokenPowerLaw(p));
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