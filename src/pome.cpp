#include "pome/pome.h"

namespace pome {

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

Pome::Pome(ModelState state) : m_state(state) {}

void Pome::buildEnergyAxis() {
  m_eSize = m_state.eSize;
  m_energyAxis = utils::LogAxis(m_state.minEnergy, m_state.maxEnergy, m_eSize);
  m_crDensity = std::vector<double>(m_eSize);
}

void Pome::buildSource() {
  m_Q = std::make_unique<Spectrum>(m_state);
  m_tau_0 = m_Q->getTau0();
}

void Pome::buildLosses() {
  m_tauEscape = std::make_unique<Escape>(m_state);
  m_losses = std::make_unique<EnergyLosses>(m_state);
}

void Pome::evolve(double dt, double t_dump, double t_max) {
  double r_pwn = 0;  // 0.1 * cgs::pc;
  double B_pwn = 0;  // 100. * cgs::muG;
  double W_B = 0;
  double E_tot = 0;
  double E_cr = 0;
  double v_pwn = 0;  // 10. * cgs::km / cgs::sec;

  utils::OutputFile out("pwn_evolution.txt");
  out << "# t [year] - r [pc] - v [km/s] - E_tot [erg] - W_B [erg] - E_cr [erg] - B [muG]\n";

  double t = 0;
  size_t iDump = 1;
  size_t counter = 0;

  while (t < t_max) {
    t += dt;
    counter++;

    r_pwn = freeExpansionPwnRadius(t, m_state.M_ejecta, m_state.rho_0, m_state.E_SN, m_state.L_0, m_tau_0);
    r_pwn = std::max(r_pwn, 1e-5 * cgs::pc);
    v_pwn = freeExpansionPwnVelocity(t, m_state.M_ejecta, m_state.rho_0, m_state.E_SN, m_state.L_0, m_tau_0);
    v_pwn = std::max(v_pwn, 1. * cgs::km / cgs::sec);
    assert(dt * v_pwn / r_pwn < 1.);
    W_B = m_state.eta_B * m_Q->getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * W_B;
    B_pwn = computeMagneticField(W_B, r_pwn);
    B_pwn = std::min(B_pwn, cgs::gauss);
    E_tot = m_Q->getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * E_tot;
    E_cr = integrateSpectrum(m_energyAxis, m_crDensity);

    auto b = std::vector<double>(m_eSize);
    std::transform(m_energyAxis.begin(), m_energyAxis.end(), b.begin(),
                   [&](double E) { return m_losses->get(E, B_pwn, r_pwn, v_pwn); });

    auto tEsc = std::vector<double>(m_eSize);
    std::transform(m_energyAxis.begin(), m_energyAxis.end(), tEsc.begin(),
                   [&](double E) { return m_tauEscape->get(E, B_pwn, r_pwn); });

    // auto timescales = computeMinTimescale(m_energyAxis, b, tEsc);

    if (counter % 100000 == 0) {
      out << std::scientific;
      out << t / cgs::year << " ";
      out << r_pwn / cgs::pc << " ";
      out << v_pwn / (cgs::km / cgs::sec) << " ";
      out << E_tot / cgs::erg << " ";
      out << W_B / cgs::erg << " ";
      out << E_cr / cgs::erg << " ";
      out << B_pwn / cgs::muG << " ";
      // out << timescales.first / cgs::year << " ";
      // out << timescales.second / cgs::year << " ";
      // out << (r_pwn / v_pwn) / cgs::year << " ";
      out << "\n";
    }

    std::vector<double> centralDiagonal(m_eSize - 1);
    std::vector<double> rhs(m_eSize - 1);
    std::vector<double> upperDiagonal(m_eSize - 2);
    std::vector<double> lowerDiagonal(m_eSize - 2);

    const double alpha = dt / 2. / (m_energyAxis[1] / m_energyAxis[0] - 1);

    for (size_t j = 0; j < m_eSize - 1; ++j) {
      double Ej = m_energyAxis[j];
      centralDiagonal.at(j) = 1. + alpha / Ej * b[j];
      if (j != m_eSize - 2) upperDiagonal.at(j) = -alpha / Ej * b[j + 1];
      rhs.at(j) = dt * m_Q->get(Ej, t);
      rhs.at(j) += (alpha / Ej * b[j + 1]) * m_crDensity[j + 1];
      // rhs.at(j) += (1. - dt / tEsc[j] - alpha / Ej * b[j]) * crDensity[j];
      rhs.at(j) += (1. - alpha / Ej * b[j]) * m_crDensity[j];
    }

    auto crDensityNew = std::vector<double>(m_eSize - 1);
    GSL::gsl_linalg_solve_tridiag(centralDiagonal, upperDiagonal, lowerDiagonal, rhs, crDensityNew);

    for (size_t j = 0; j < m_eSize - 1; ++j) {
      m_crDensity[j] = crDensityNew[j];
    }

    for (size_t j = 0; j < m_eSize - 1; ++j) {
      double Ej = m_energyAxis[j];
      crDensityNew[j] = m_crDensity[j] * std::exp(-dt / tEsc[j]) +
                        0.5 * dt * (m_Q->get(Ej, t) + m_Q->get(Ej, t + dt) * std::exp(-dt / tEsc[j]));
    }

    if (t > (double)iDump * t_dump) {
      dumpSpectrum(m_energyAxis, m_crDensity, iDump);
      auto maxDensity = *std::max_element(m_crDensity.begin(), m_crDensity.end());
      std::cout << t / cgs::year << " " << m_crDensity.size() << " " << maxDensity << "\n";
      iDump++;
    }
  }
}

}  // namespace pome