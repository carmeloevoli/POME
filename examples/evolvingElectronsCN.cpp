#include <algorithm>
#include <cassert>

#include "pome.h"

using namespace pome;

void dumpSpectrum(std::vector<double> E, std::vector<double> density, size_t i) {
  std::string filename = "evolved_spectrum_cn_" + std::to_string(i) + ".txt";
  utils::OutputFile out(filename.c_str());
  out << std::scientific;
  auto eSize = E.size();
  for (size_t j = 0; j < eSize; ++j) {
    out << E[j] / cgs::TeV << "\t";
    out << (E[j] * E[j] * density[j]) / cgs::TeV << "\t";
    out << "\n";
  }
}

double integrateSpectrum(std::vector<double> E, std::vector<double> density) {
  const auto ln_f = std::log(E[1] / E[0]);
  const auto N = density.size();
  double value = 0;
  for (size_t i = 1; i < N - 1; ++i) {
    value += 2. * pow2(E[i]) * density[i];
  }
  value += pow2(E[0]) * density[0];
  value += pow2(E[N - 1]) * density[N - 1];
  return 0.5 * ln_f * value;
}

std::pair<double, double> computeMinTimescale(std::vector<double> E, std::vector<double> b,
                                              std::vector<double> tEscape) {
  auto eSize = E.size();
  std::vector<double> tLosses;
  for (size_t j = 0; j < eSize; ++j) {
    tLosses.push_back(E[j] / b[j]);
  }
  auto minLosses = *std::min_element(tLosses.begin(), tLosses.end());
  auto minEscape = *std::min_element(tEscape.begin(), tEscape.end());
  return std::make_pair(minLosses, minEscape);
}

int main() {
  try {
    Timer timer;
    ModelState state;

    auto Q = Spectrum(state);
    auto tau_0 = Q.getTau0();
    auto tauEscape = Escape(state);
    auto losses = EnergyLosses(state);

    const size_t eSize = 6 * 32;
    auto energyAxis = utils::LogAxis(1. * cgs::GeV, 1. * cgs::PeV, eSize);
    auto crDensity = std::vector<double>(eSize);

    double t = 0;
    const double dt = 1e-4 * cgs::year;

    double r_pwn = 0.1 * cgs::pc;
    double B_pwn = 100. * cgs::muG;
    double W_B = 0;
    double E_tot = 0;
    double E_tot_exact = 0;
    double v_pwn = 10. * cgs::km / cgs::sec;  // TODO check this!

    utils::OutputFile out("pwn_evolution.txt");
    out << "# t - r - v - E_tot - W_B - B - t_loss - t_esc\n";

    for (size_t i = 0; i < 100001; ++i) {
      t += dt;

      r_pwn = freeExpansionPwnRadius(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
      r_pwn = std::max(r_pwn, 1e-5 * cgs::pc);
      v_pwn = freeExpansionPwnVelocity(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
      v_pwn = std::max(v_pwn, 1. * cgs::km / cgs::sec);
      assert(dt * v_pwn / r_pwn < 1.);
      W_B = state.eta_B * Q.getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * W_B;
      B_pwn = computeMagneticField(W_B, r_pwn);
      B_pwn = std::min(B_pwn, cgs::gauss);

      E_tot = Q.getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * E_tot;
      E_tot_exact = W_B + integrateSpectrum(energyAxis, crDensity);

      auto b = std::vector<double>(eSize);
      std::transform(energyAxis.begin(), energyAxis.end(), b.begin(),
                     [&](double E) { return losses.get(E, B_pwn, r_pwn, v_pwn); });

      auto tEsc = std::vector<double>(eSize);
      std::transform(energyAxis.begin(), energyAxis.end(), tEsc.begin(),
                     [&](double E) { return tauEscape.get(E, B_pwn, r_pwn); });

      auto timescales = computeMinTimescale(energyAxis, b, tEsc);

      if (i % 100 == 0) {
        out << std::scientific;
        out << t / cgs::year << " ";
        out << r_pwn / cgs::pc << " ";
        out << v_pwn / (cgs::km / cgs::sec) << " ";
        out << E_tot / cgs::erg << " ";
        out << W_B / cgs::erg << " ";
        out << E_tot_exact / cgs::erg << " ";
        out << B_pwn / cgs::muG << " ";
        out << timescales.first / cgs::year << " ";
        out << timescales.second / cgs::year << " ";
        out << (r_pwn / v_pwn) / cgs::year << " ";
        out << "\n";
      }

      std::vector<double> centralDiagonal(eSize - 1);
      std::vector<double> rhs(eSize - 1);
      std::vector<double> upperDiagonal(eSize - 2);
      std::vector<double> lowerDiagonal(eSize - 2);

      const double alpha = dt / 2. / (energyAxis[1] / energyAxis[0] - 1);

      for (size_t j = 0; j < eSize - 1; ++j) {
        double Ej = energyAxis[j];
        centralDiagonal.at(j) = 1. + alpha / Ej * b[j];
        if (j != eSize - 2) upperDiagonal.at(j) = -alpha / Ej * b[j + 1];
        rhs.at(j) = dt * Q.get(Ej, t);
        rhs.at(j) += (alpha / Ej * b[j + 1]) * crDensity[j + 1];
        // rhs.at(j) += (1. - dt / tEsc[j] - alpha / Ej * b[j]) * crDensity[j];
        rhs.at(j) += (1. - alpha / Ej * b[j]) * crDensity[j];
      }

      auto crDensityNew = std::vector<double>(eSize - 1);
      GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, crDensityNew);

      for (size_t j = 0; j < eSize - 1; ++j) {
        crDensity[j] = crDensityNew[j];
      }

      for (size_t j = 0; j < eSize - 1; ++j) {
        double Ej = energyAxis[j];
        crDensityNew[j] = crDensity[j] * std::exp(-dt / tEsc[j]) +
                          0.5 * dt * (Q.get(Ej, t) + Q.get(Ej, t + dt) * std::exp(-dt / tEsc[j]));
      }

      if (i % 10000 == 0) {
        dumpSpectrum(energyAxis, crDensity, i / 10000);
        auto maxDensity = *std::max_element(crDensity.begin(), crDensity.end());
        std::cout << i << " " << crDensity.size() << " " << maxDensity << "\n";
      }
    }
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}