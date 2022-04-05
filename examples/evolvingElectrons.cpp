#include <algorithm>

#include "pome.h"

using namespace pome;

void dumpSpectrum(std::vector<double> E, std::vector<double> density, size_t i) {
  std::string filename = "evolved_spectrum_simple_" + std::to_string(i) + ".txt";
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

    auto Q = Spectrum(state);
    auto tauEscape = Escape(state);
    auto losses = EnergyLosses(state);

    const size_t eSize = 6 * 32;
    auto energyAxis = utils::LogAxis(1. * cgs::GeV, 1. * cgs::PeV, eSize);
    auto crDensity = std::vector<double>(eSize);

    double t = 0;
    const double dt = cgs::year;

    for (size_t i = 0; i < 10001; ++i) {
      t += dt;

      const double r_pwn = cgs::pc;
      const double B_pwn = 10. * cgs::muG;
      const double v_pwn = 1. * cgs::km / cgs::sec;  // TODO check this!

      if (i == 0) {
        std::cout << std::scientific;
        std::cout << t / cgs::year << " year\n";
        std::cout << r_pwn / cgs::pc << " pc\n";
        std::cout << v_pwn / (cgs::km / cgs::sec) << " km/sec\n";
        std::cout << (r_pwn / v_pwn) / cgs::kyr << " kyr\n";
        //   std::cout << E_tot / cgs::erg << "\t";
        //   std::cout << W_B / cgs::erg << "\t";
        std::cout << B_pwn / cgs::muG << " muG\n";
        std::cout << "\n";
      }

      // r_pwn = freeExpansionPwnRadius(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
      // r_pwn = std::max(r_pwn, 1e-5 * cgs::pc);
      // v_pwn = freeExpansionPwnVelocity(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
      // v_pwn = std::max(v_pwn, 10. * cgs::km / cgs::sec);
      // E_tot = Q.getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * E_tot;
      // W_B = state.eta_B * Q.getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * W_B;
      // B_pwn = computeMagneticField(W_B, r_pwn);
      // B_pwn = std::max(B_pwn, cgs::gauss);

      auto b = std::vector<double>(eSize);
      std::transform(energyAxis.begin(), energyAxis.end(), b.begin(),
                     [&](double E) { return losses.get(E, B_pwn, r_pwn, v_pwn); });

      auto tEsc = std::vector<double>(eSize);
      std::transform(energyAxis.begin(), energyAxis.end(), tEsc.begin(),
                     [&](double E) { return tauEscape.get(E, B_pwn, r_pwn); });

      auto crDensityNew = std::vector<double>(eSize);

      // std::cout << computeMinTimescale(energyAxis, b, tEsc) / cgs::year << " yr\n";

      // std::vector<double> centralDiagonal(eSize - 1);
      // std::vector<double> rhs(eSize - 1);
      // std::vector<double> upperDiagonal(eSize - 2);
      // std::vector<double> lowerDiagonal(eSize - 2);

      // const double f = energyAxis[1] / energyAxis[0];
      // const double alpha = dt / 2. / (f - 1);

      // for (size_t j = 0; j < eSize - 1; ++j) {
      //   double Ej = energyAxis[j];
      //   centralDiagonal.at(j) = 1. + alpha / Ej;
      // }

      for (size_t j = 0; j < eSize - 1; ++j) {
        auto E = energyAxis[j];
        auto Eup = energyAxis[j + 1];
        crDensityNew[j] = crDensity[j] + dt * Q.get(E, t);
        crDensityNew[j] -= crDensity[j] * dt / tEsc[j];
        crDensityNew[j] += dt / (Eup - E) * (b[j + 1] * crDensity[j + 1] - b[j] * crDensity[j]);
      }
      crDensity.assign(crDensityNew.begin(), crDensityNew.end());

      if (i % 1000 == 0) {
        dumpSpectrum(energyAxis, crDensity, i / 1000);
        std::cout << i << " " << crDensity.size() << " " << *std::max_element(crDensity.begin(), crDensity.end())
                  << "\n";
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}