#include "pome.h"

#include <fstream>

namespace POME {

void Pome::buildSource() {
  std::pair<double, double> alpha = std::make_pair(m_state.alpha_1, m_state.alpha_2);
  BrokenPowerLawParams p = {alpha, m_state.E_b, m_state.V_0, m_state.tau_0, m_state.L_B, m_state.n};
  m_Q = std::unique_ptr<AbstractSource>(new BrokenPowerLaw(p));
}

void Pome::dumpSourceSpectrum() const {
  auto energyAxis = utils::LogAxis(1e17 * SI::eV, 1e22 * SI::eV, 500);
  utils::OutputFile out("source_spectrum.txt");
  out << std::scientific;
  for (auto E : energyAxis) {
    out << E / SI::eV << "\t";
    // auto Gamma = E / SI::protonMassC2;
    // out << SI::cLight / adiabatic.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
    // out << SI::cLight / losses.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
    // out << SI::cLight / pair_cmb.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
    // out << SI::cLight / pair_irb.dlnGamma_dt(proton, Gamma) / SI::Mpc << "\t";
    out << "\n";
  }

  // {
  //   std::ofstream outfile("output/source_spectrum.txt");
  //   if (outfile.is_open()) {
  //     const double units = cgs::TeV / cgs::sec / cgs::cm2;
  //     for (double E = 1e-6 * cgs::TeV; E < 1e6 * cgs::TeV; E *= 1.1) {
  //       outfile << E / cgs::TeV << " " << (E * E * m_Q->get(E, 0.)) / units << "\n";
  //     }
  //     outfile.close();
  //   } else {
  //     throw std::runtime_error("output file cannot be opened.");
  //   }
  // }
  // {
  //   std::ofstream outfile("output/source_evolution.txt");
  //   if (outfile.is_open()) {
  //     const double units = 1. / cgs::TeV / cgs::sec / cgs::cm2;
  //     for (double t = cgs::year; t < 1e4 * cgs::year; t *= 1.1) {
  //       outfile << t / cgs::year << " " << m_Q->get(cgs::TeV, t) / units << "\n";
  //     }
  //     outfile.close();
  //   } else {
  //     throw std::runtime_error("output file cannot be opened.");
  //   }
  // }
}

void Pome::evolve() const {
  double rPWN = 10. * cgs::pc;
  double BPWN = 0.;
  for (double t = 0; t < 1e6 * cgs::year; t += 1e2 * cgs::year) {
    std::cout << t / cgs::year << " " << rPWN / cgs::pc << " " << BPWN / cgs::muG << "\n";
  }
}

}  // namespace POME