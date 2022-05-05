
utils::OutputFile out("solution.txt");
out << std::scientific;
for (size_t j = 0; j < eSize; ++j) {
  out << energyAxis[j] / cgs::TeV << "\t";
  auto E = energyAxis[j];
  out << Q->get(E) * tauEscape->get(E, B, rPWN) / (1. / cgs::TeV / cgs::sec) << "\t";
  out << Q->get(E) * E / losses->get(E, B, tauAdiabatic) / (1. / cgs::TeV / cgs::sec) << "\t";
  out << "\n";
}

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

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

double computeMagneticField(double magneticFieldEnergy, double R) {
  auto V = (4. / 3.) * M_PI * pow3(R);
  auto w_B = magneticFieldEnergy / V;
  return std::sqrt(8. * M_PI * w_B);
}

struct odeparams {
  size_t eSize;
  const std::vector<double> &E;
  const ModelState &state;
  const Spectrum &source;
  const Escape &tauEscape;
  const EnergyLosses &losses;
};

int func(double t, const double y[], double f[], void *params) {
  odeparams p = *(odeparams *)params;

  const auto size = p.eSize - 1;
  // const auto W_B = y[eSize - 1];
  // const auto tau0 = p.source.getTau0();
  const auto s = p.state;

  const auto r_pwn = cgs::pc;
  // freeExpansionPwnRadius(t, s.M_ejecta, s.rho_0, s.E_SN, s.L_0, tau0);

  const auto v_pwn = 10. * cgs::km / cgs::sec;
  // freeExpansionPwnVelocity(t, s.M_ejecta, s.rho_0, s.E_SN, s.L_0, tau0);

  const auto B_pwn = 10. * cgs::muG;
  // (W_B, r_pwn);

  for (size_t i = 0; i < size; ++i) {
    auto Q = p.source.get(p.E[i], t);
    // auto tEsc = p.tauEscape.get(p.E[i], B_pwn, r_pwn);
    auto bUp = p.losses.get(p.E[i + 1], B_pwn, r_pwn, v_pwn);
    auto b = p.losses.get(p.E[i], B_pwn, r_pwn, v_pwn);
    auto DeltaE = p.E[i + 1] - p.E[i];
    f[i] = Q + bUp / DeltaE * y[i + 1] - (b / DeltaE /*+ 1 / tEsc */) * y[i];
  }
  // f[eSize - 1] = s.eta_B * p.source.getLuminosity(t) - (v_pwn / r_pwn) * W_B;
  return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
  odeparams p = *(odeparams *)params;

  const auto size = p.eSize - 1;

  const auto r_pwn = cgs::pc;
  // freeExpansionPwnRadius(t, s.M_ejecta, s.rho_0, s.E_SN, s.L_0, tau0);

  const auto v_pwn = 10. * cgs::km / cgs::sec;
  // freeExpansionPwnVelocity(t, s.M_ejecta, s.rho_0, s.E_SN, s.L_0, tau0);

  const auto B_pwn = 10. * cgs::muG;
  // (W_B, r_pwn);

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, size, size);
  gsl_matrix *m = &dfdy_mat.matrix;

  for (size_t i = 0; i < size; i++)
    for (size_t j = 0; j < size; j++) {
      if (i == j) {
        auto b = p.losses.get(p.E[i], B_pwn, r_pwn, v_pwn);
        auto DeltaE = p.E[i + 1] - p.E[i];
        auto J = -b / DeltaE;
        gsl_matrix_set(m, i, j, J);
      } else if (i + 1 == j) {
        auto bUp = p.losses.get(p.E[i + 1], B_pwn, r_pwn, v_pwn);
        auto DeltaE = p.E[i + 1] - p.E[i];
        auto J = bUp / DeltaE;
        gsl_matrix_set(m, i, j, J);
      } else {
        auto J = 0.0;
        gsl_matrix_set(m, i, j, J);
      }
    }

  for (size_t i = 0; i < size; i++) dfdt[i] = 0.0;  // TODO wrong!

  return GSL_SUCCESS;
}

int odeint(void) {
  // double mu = 10;
  // gsl_odeiv2_system sys = {func, jac, 2, &mu};

  // gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);  //
  // int i;
  // double t = 0.0, t1 = 100.0;
  // double y[2] = {1.0, 0.0};

  // for (i = 1; i <= 100; i++) {
  //   double ti = i * t1 / 100.0;
  //   int status = gsl_odeiv2_driver_apply(d, &t, ti, y);

  //   if (status != GSL_SUCCESS) {
  //     printf("error, return value=%d\n", status);
  //     break;
  //   }

  //   printf("%.5e %.5e %.5e\n", t, y[0], y[1]);
  // }

  // gsl_odeiv2_driver_free(d);
  return 0;
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

    odeparams params = {eSize, energyAxis, state, Q, tauEscape, losses};
    gsl_odeiv2_system sys = {func, jac, eSize, &params};

    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    double y[eSize - 1];
    for (size_t i = 0; i < eSize - 1; i++) y[i] = 0.0;

    double t = 0.;
    for (double t_i = 0.0; t_i < 10. * cgs::kyr; t_i += cgs::kyr) {
      std::cout << t_i / cgs::kyr << "\n";
      int status = gsl_odeiv2_driver_apply(d, &t, t_i, y);

      if (status != GSL_SUCCESS) {
        std::cout << "error, return value=" << status << "\n";
        exit(1);
      }
    }
    gsl_odeiv2_driver_free(d);
    // auto crDensity = std::vector<double>(eSize);

    // double dt = 1e-3 * cgs::year;
    // double t = 0;

    // double r_pwn = 0;
    // double v_pwn = 0;
    // double E_tot = 0;
    // double W_B = 0;
    // double B_pwn = 0;

    // for (size_t i = 1; i < 2001; ++i) {
    //   t += dt;

    //   std::cout << std::scientific;
    //   std::cout << t / cgs::year << "\t";
    //   std::cout << r_pwn / cgs::pc << "\t";
    //   std::cout << v_pwn / (cgs::km / cgs::sec) << "\t";
    //   std::cout << E_tot / cgs::erg << "\t";
    //   std::cout << W_B / cgs::erg << "\t";
    //   std::cout << B_pwn / cgs::muG << "\t";
    //   std::cout << "\n";

    // r_pwn = freeExpansionPwnRadius(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
    // r_pwn = std::max(r_pwn, 1e-5 * cgs::pc);
    // v_pwn = freeExpansionPwnVelocity(t, state.M_ejecta, state.rho_0, state.E_SN, state.L_0, tau_0);
    // v_pwn = std::max(v_pwn, 10. * cgs::km / cgs::sec);
    // E_tot = Q.getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * E_tot;
    // W_B = state.eta_B * Q.getLuminosity(t) * dt + (1. - dt * v_pwn / r_pwn) * W_B;
    // B_pwn = computeMagneticField(W_B, r_pwn);
    // B_pwn = std::max(B_pwn, cgs::gauss);

    // auto b = std::vector<double>(eSize);
    // std::transform(energyAxis.begin(), energyAxis.end(), b.begin(),
    //                [&](double E) { return losses.get(E, B_pwn, r_pwn / v_pwn); });

    // auto tEsc = std::vector<double>(eSize);
    // std::transform(energyAxis.begin(), energyAxis.end(), tEsc.begin(),
    //                [&](double E) { return tauEscape.get(E, B_pwn, r_pwn); });

    // auto crDensityNew = std::vector<double>(eSize);

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

    // for (size_t j = 0; j < eSize - 1; ++j) {
    //   auto E = energyAxis[j];
    //   auto Eup = energyAxis[j + 1];
    //   crDensityNew[j] = crDensity[j] + dt * Q->get(E, t);
    //   crDensityNew[j] -= crDensity[j] * dt / tEscape[j];
    //   crDensityNew[j] += dt / (Eup - E) * (bLosses[j + 1] * crDensity[j + 1] - bLosses[j] * crDensity[j]);
    // }
    // crDensity.assign(crDensityNew.begin(), crDensityNew.end());

    // if (i % 1000 == 0) {
    //   dumpSpectrum(energyAxis, crDensity, i / 1000);
    //   std::cout << i << " " << crDensity.size() << " " << *std::max_element(crDensity.begin(), crDensity.end())
    //             << "\n";
    // }

    // year\n";
    //}
  } catch (const std::exception &e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}