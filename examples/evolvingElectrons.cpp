#include <memory>

#include "pome.h"

using namespace pome;

void dumpSpectrum(std::vector<double> E, std::vector<double> density, size_t i) {
  std::string filename = "evolved_spectrum_" + std::to_string(i) + ".txt";
  utils::OutputFile out(filename.c_str());
  out << std::scientific;
  auto eSize = E.size();
  for (size_t j = 0; j < eSize; ++j) {
    out << E[j] / cgs::GeV << "\t";
    out << density[j] / (1. / cgs::GeV / cgs::sec) << "\t";
    // out << (Q->get(energyAxis[j]) * tauEscape->get(energyAxis[j])) / (1. / cgs::GeV / cgs::sec) << "\t";
    out << "\n";
  }
}

int main() {
  try {
    Timer timer;
    ModelState state;

    std::pair<double, double> alpha = std::make_pair(state.alpha_1, state.alpha_2);
    BrokenPowerLawParams p = {alpha, state.E_b, state.V_0, state.tau_0, state.L_B, state.n};
    auto Q = std::unique_ptr<AbstractSource>(new BrokenPowerLaw(p));

    auto tauEscape = std::unique_ptr<Escape>(new Escape({10. * cgs::muG, 1. * cgs::pc}));
    auto losses = std::unique_ptr<Losses>(new Losses({10. * cgs::muG}));

    const size_t eSize = 100;
    const double dt = 0.01 * cgs::kyr;

    auto energyAxis = utils::LogAxis(1. * cgs::GeV, 1. * cgs::PeV, eSize);
    auto crDensity = std::vector<double>(eSize);

    for (size_t i = 0; i < 1000; ++i) {
      auto crDensityNew = std::vector<double>(eSize);
      for (size_t j = 0; j < eSize; ++j) {
        auto E = energyAxis[j];
        crDensityNew[j] = dt * Q->get(E) + crDensity[j] * (1. - dt / tauEscape->get(E));
      }
      crDensity.assign(crDensityNew.begin(), crDensityNew.end());

      if (i % 100 == 0) dumpSpectrum(energyAxis, crDensity, i);

      std::cout << i << " " << crDensity.size() << "\n";
    }

    utils::OutputFile out("solution.txt");
    out << std::scientific;
    for (size_t j = 0; j < eSize; ++j) {
      out << energyAxis[j] / cgs::GeV << "\t";
      out << (Q->get(energyAxis[j]) * tauEscape->get(energyAxis[j])) / (1. / cgs::GeV / cgs::sec) << "\t";
      out << "\n";
    }

  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}