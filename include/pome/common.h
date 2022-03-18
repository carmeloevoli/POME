#ifndef POME_PWN_H
#define POME_PWN_H

namespace pome {

double freeExpansionPwnRadius(const double& t, const double& M_ej, const double& rho_0, const double& E_SN,
                              const double& L_pulsar, const double& tau_0);

double potentialDropEnergy(const double& Omega_0, const double& BSurface, const double& pulsarRadius);

double decayTime(const double& Omega_0, const double& BSurface, const double& pulsarRadius, const double& I);

double luminosityDecay(const double& t, const double& tau_0, const double& n);

}  // namespace pome

#endif  // POME_PWN_H