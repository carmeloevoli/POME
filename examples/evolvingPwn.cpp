#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <algorithm>
#include <memory>

#include "pome.h"

using namespace pome;

void printTimescales(double B_pwn, double r_pwn, double v_pwn) {
  ModelState state;
  auto tauEscape = std::unique_ptr<Escape>(new Escape(state));
  auto losses = std::unique_ptr<EnergyLosses>(new EnergyLosses(state));
  auto energyAxis = utils::LogAxis(1. * cgs::GeV, 1. * cgs::PeV, 600);

  utils::OutputFile out("timescales.txt");
  out << std::scientific;
  for (auto E : energyAxis) {
    out << E / cgs::GeV << "\t";
    out << tauEscape->get(E, B_pwn, r_pwn) / (cgs::kyr) << "\t";
    out << E / losses->bLosses(E, B_pwn) / (cgs::kyr) << "\t";
    out << E / losses->bAdiabatic(E, r_pwn, v_pwn) / (cgs::kyr) << "\t";
    out << E / losses->get(E, B_pwn, r_pwn, v_pwn) / (cgs::kyr) << "\t";
    out << "\n";
  }
}

int main() {
  try {
    Timer timer;
    ModelState state;
    auto Q = Spectrum(state);

    const double dt = 0.1 * cgs::year;
    const double tau_0 = Q.getTau0();

    double r_pwn = 0;
    double v_pwn = 0;
    double E_tot = 0;
    double W_B = 0;
    double B = 0;

    utils::OutputFile out("pwn_evolution.txt");
    out << std::scientific;
    for (size_t i = 1; i < 1001; ++i) {
      const auto t = (double)i * dt;
      out << t / cgs::year << "\t";
      out << r_pwn / cgs::pc << "\t";
      out << v_pwn / (cgs::km / cgs::sec) << "\t";
      out << E_tot / cgs::erg << "\t";
      out << W_B / cgs::erg << "\t";
      out << B / cgs::muG << "\t";
      out << "\n";

      r_pwn = freeExpansionPwnRadius(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
      r_pwn = std::max(r_pwn, 1e-5 * cgs::pc);
      v_pwn = freeExpansionPwnVelocity(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
      v_pwn = std::max(v_pwn, 10. * cgs::km / cgs::sec);
      E_tot = Q.getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * E_tot;
      W_B = state.eta_B * Q.getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * W_B;
      B = computeMagneticField(W_B, r_pwn);
    }
    std::cout << "r : " << r_pwn / cgs::parsec << " pc\n";
    std::cout << "v : " << v_pwn / (cgs::km / cgs::sec) << " km/s\n";
    std::cout << "t_ad : " << (r_pwn / v_pwn) / cgs::year << " yr\n";
    std::cout << "B : " << B / cgs::muG << " muG\n";

    printTimescales(B, r_pwn, v_pwn);
  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}