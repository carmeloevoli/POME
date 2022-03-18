#include <algorithm>
#include <memory>

#include "pome.h"

using namespace pome;

void dumpSpectrum(std::vector<double> E, std::vector<double> density, size_t i) {
  std::string filename = "evolved_spectrum_" + std::to_string(i) + ".txt";
  utils::OutputFile out(filename.c_str());
  out << std::scientific;
  auto eSize = E.size();
  for (size_t j = 0; j < eSize; ++j) {
    out << E[j] / cgs::TeV << "\t";
    out << (E[j] * E[j] * density[j]) / cgs::TeV << "\t";
    out << "\n";
  }
}

int main() {
  try {
    Timer timer;
    ModelState state;

    auto Q = std::unique_ptr<AbstractSource>(new BrokenPowerLaw(state));
    auto tauEscape = std::unique_ptr<Escape>(new Escape(state));
    auto losses = std::unique_ptr<Losses>(new Losses(state));

    const size_t eSize = 6 * 32;
    const double dt = cgs::year;

    auto energyAxis = utils::LogAxis(1. * cgs::GeV, 1. * cgs::PeV, eSize);
    auto crDensity = std::vector<double>(eSize);

    const double B = 10. * cgs::muG;
    const double rPWN = cgs::pc;
    const double tauAdiabatic = 1e3 * cgs::kyr;

    for (size_t i = 0; i < 20001; ++i) {
      const auto t = (double)i * dt;
      auto crDensityNew = std::vector<double>(eSize);
      for (size_t j = 0; j < eSize - 1; ++j) {
        auto E = energyAxis[j];
        auto Eup = energyAxis[j + 1];
        crDensityNew[j] = crDensity[j] + dt * Q->get(E, t);
        crDensityNew[j] -= crDensity[j] * dt / tauEscape->get(E, B, rPWN);
        crDensityNew[j] +=
            dt / (Eup - E) *
            (losses->get(Eup, B, tauAdiabatic) * crDensity[j + 1] - losses->get(E, B, tauAdiabatic) * crDensity[j]);
      }
      crDensity.assign(crDensityNew.begin(), crDensityNew.end());

      if (i % 1000 == 0) {
        dumpSpectrum(energyAxis, crDensity, i / 1000);
        std::cout << i << " " << crDensity.size() << " " << *std::max_element(crDensity.begin(), crDensity.end())
                  << "\n";
      }
    }

    utils::OutputFile out("solution.txt");
    out << std::scientific;
    for (size_t j = 0; j < eSize; ++j) {
      out << energyAxis[j] / cgs::TeV << "\t";
      auto E = energyAxis[j];
      out << Q->get(E) * tauEscape->get(E, B, rPWN) / (1. / cgs::TeV / cgs::sec) << "\t";
      out << Q->get(E) * E / losses->get(E, B, tauAdiabatic) / (1. / cgs::TeV / cgs::sec) << "\t";
      out << "\n";
    }

  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}