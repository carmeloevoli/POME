#include "pome.h"

using namespace pome;

class CoreCollapse {
 protected:
  const double adiabaticIndex = 5. / 3.;
  const double mu = 1.4;
  const double snEnergy = 1e51 * cgs::erg;
  const double massLossRate = 1e-5 * cgs::sunMass / cgs::year;
  const double ejectaMass = 5. * cgs::sunMass;
  const double windSpeed = 1e6 * cgs::cm / cgs::sec;
  const double windLuminosity = 1e36 * cgs::erg / cgs::sec;
  const double windLifetime = 1. * cgs::Myr;
  const double ismDensity = 1. / cgs::cm3;

  double r_w, r_b, rho_b;

  double computeWindRadius() {
    const auto M_5 = massLossRate / (1e-5 * cgs::sunMass / cgs::year);
    const auto u_w = windSpeed / (1e6 * cgs::cm / cgs::sec);
    const auto L_w = windLuminosity / (1e36 * cgs::erg / cgs::sec);
    const auto n = ismDensity / (1. / cgs::cm3);
    const auto t_w = windLifetime / (1. * cgs::Myr);

    auto value = 1.3 * cgs::pc * std::pow(mu / 1.4, -0.5);
    value *= std::pow(M_5 * u_w, 0.5);
    value *= std::pow(L_w, -7. / 35.);
    value *= std::pow(n, -21. / 70.);
    value *= std::pow(t_w, 14. / 35.);
    return value;
  }

  double computeBubbleRadius() {
    const auto L_w = windLuminosity / (1e36 * cgs::erg / cgs::sec);
    const auto n = ismDensity / (1. / cgs::cm3);
    const auto t_w = windLifetime / (1. * cgs::Myr);

    auto value = 28. * cgs::pc;
    value *= std::pow(L_w, 1. / 5.);
    value *= std::pow(n, -1. / 5.);
    value *= std::pow(t_w, 3. / 5.);
    return value;
  }

  double computeBubbleDensity() {
    const auto L_w = windLuminosity / (1e36 * cgs::erg / cgs::sec);
    const auto n = ismDensity / (1. / cgs::cm3);
    const auto t_w = windLifetime / (1. * cgs::Myr);

    auto n_b = 0.01 / cgs::cm3;
    n_b *= std::pow(L_w, 6. / 35.);
    n_b *= std::pow(n, 19. / 35.);
    n_b *= std::pow(t_w, -22. / 35.);
    return cgs::protonMass * mu * n_b;
  }

 public:
  CoreCollapse() {
    r_w = computeWindRadius();
    r_b = computeBubbleRadius();
    rho_b = computeBubbleDensity();
  }
  virtual ~CoreCollapse() = default;

  double Mass(double r) {
    constexpr auto factor = 4. * M_PI / 3.;
    const auto rho_ism = cgs::protonMass * mu * ismDensity;
    auto value = ejectaMass;
    if (r < r_w) {
      value += massLossRate / windSpeed * r;
    } else if (r < r_b) {
      value += massLossRate / windSpeed * r_w;
      value += factor * rho_b * (pow3(r) - pow3(r_w));
    } else {
      value += massLossRate / windSpeed * r_w;
      value += factor * rho_b * (pow3(r_b) - pow3(r_w));
      value += factor * rho_ism * (pow3(r) - pow3(r_b));
    }
    return value;
  }

  double shockVelocity(double R_s) {
    const auto alpha = 6. * (adiabaticIndex - 1.) / (adiabaticIndex + 1.);
    const auto v0 = 0.5 * (adiabaticIndex + 1.) * std::sqrt(2. * alpha * snEnergy / Mass(R_s));
    if (R_s > 1e-3 * cgs::pc) {
      const auto I =
          utils::QAGIntegration<double>([&](double r) { return std::pow(r, alpha - 1.) * Mass(r); }, 0., R_s, 1000);
      return v0 / std::pow(R_s, alpha) / Mass(R_s) * I;
    } else {
      return v0;
    }
  }
};

int main() {
  try {
    Timer timer;
    CoreCollapse cc;

    {
      utils::OutputFile out("cc_mass.txt");
      out << "# r [pc] - Mass [M_odot]\n";
      out << std::scientific;
      for (double r = 0.; r < 100. * cgs::parsec; r += 0.1 * cgs::parsec) {
        out << r / cgs::parsec << " " << cc.Mass(r) / cgs::sunMass << "\n";
      }
    }

    {
      utils::OutputFile out("cc_velocity.txt");
      out << "# t [kyr] - R_s [pc] - u_s [cm/s]\n";
      out << std::scientific;
      double R_s = 0.;
      const double dt = 1. * cgs::year;
      size_t counter = 0;
      for (double t = 0.; t < 5e2 * cgs::kyr; t += dt) {
        auto u_s = cc.shockVelocity(R_s);
        R_s += dt * u_s;
        if (counter % 10 == 0) out << t / cgs::kyr << " " << R_s / cgs::pc << " " << u_s / (cgs::cm / cgs::sec) << "\n";
        counter++;
      }
    }

  } catch (const std::exception& e) {
    std::cerr << "exception caught with message: " << e.what() << "\n";
  }
  return EXIT_SUCCESS;
}