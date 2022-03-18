#include <algorithm>
#include <memory>

#include "pome.h"

using namespace pome;

double computeMagneticField(double magneticFieldEnergy, double R) {
  auto V = (4. / 3.) * M_PI * pow3(R);
  auto w_B = magneticFieldEnergy / V;
  return std::sqrt(8. * M_PI * w_B);
}

int main() {
  try {
    Timer timer;
    ModelState state;

    auto Q = std::unique_ptr<AbstractSource>(new BrokenPowerLaw(state));
    const double dt = cgs::year;
    const double tau_0 = 20. * cgs::kyr;

    double r_pwn = 0;
    double v_pwn = 0;
    double E_tot = 0;
    double W_B = 0;
    double B = 0;

    utils::OutputFile out("pwn_evolution.txt");
    out << std::scientific;
    for (size_t i = 0; i < 20001; ++i) {
      const auto t = (double)i * dt;
      out << t / cgs::year << "\t";
      out << r_pwn / cgs::pc << "\t";
      out << v_pwn / (cgs::km / cgs::sec) << "\t";
      out << E_tot / cgs::erg << "\t";
      out << W_B / cgs::erg << "\t";
      out << B / cgs::muG << "\t";
      out << "\n";

      r_pwn = freeExpansionPwnRadius(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
      auto r_pwn_up = freeExpansionPwnRadius(t + dt, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
      v_pwn = (r_pwn_up - r_pwn) / dt;
      auto L_pwn = state.L_0 * luminosityDecay(t, tau_0, state.n);
      E_tot = L_pwn * dt + (1. - dt * v_pwn / r_pwn) * E_tot;
      W_B = state.eta_B * L_pwn * dt + (1. - dt * v_pwn / r_pwn) * W_B;
      B = computeMagneticField(W_B, r_pwn);
    }
  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}