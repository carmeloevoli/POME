#ifndef INCLUDE_CGS_H_
#define INCLUDE_CGS_H_

#include <cmath>

#define pow2(A) ((A) * (A))
#define pow3(A) ((A) * (A) * (A))
#define pow4(A) ((A) * (A) * (A) * (A))

namespace cgs {

// CGS UNITS
static constexpr double second = 1;
static constexpr double centimeter = 1.;
static constexpr double gram = 1;
static constexpr double erg = 1;
static constexpr double kelvin = 1;
static constexpr double sr = 1;
static constexpr double gauss = 1;
static constexpr double esu = 1;
static constexpr double statvolt = erg / esu;
static constexpr double steradiant = 1;

// TIME UNITS
static constexpr double msec = 1e-3 * second;
static constexpr double year = 3.154e+7 * second;
static constexpr double kiloyear = 1e3 * year;
static constexpr double Megayear = 1e6 * year;
static constexpr double Gigayear = 1e9 * year;

// LENGTH UNITS
static constexpr double fm = 1e-13 * centimeter;
static constexpr double meter = 1e2 * centimeter;
static constexpr double kilometer = 1e3 * meter;
static constexpr double parsec = 3.086e16 * meter;
static constexpr double kiloparsec = 1e3 * parsec;

// MASS UNITS
static constexpr double mgram = 1e-3 * gram;
static constexpr double kilogram = 1e3 * gram;

// ENERGY UNITS
static constexpr double joule = 1e7 * erg;
static constexpr double electronvolt = 1.60217657e-19 * joule;
static constexpr double kiloelectronvolt = 1e3 * electronvolt;
static constexpr double megaelectronvolt = 1e6 * electronvolt;
static constexpr double gigaelectronvolt = 1e9 * electronvolt;
static constexpr double teraelectronvolt = 1e12 * electronvolt;
static constexpr double petaelectronvolt = 1e15 * electronvolt;

// EM UNITS
static constexpr double nanogauss = 1e-9 * gauss;
static constexpr double microgauss = 1e-6 * gauss;
static constexpr double milligauss = 1e-3 * gauss;
static constexpr double volt = statvolt / 299.792458;

// ABBREVIATION
static constexpr double sec = second;
static constexpr double km = kilometer;
static constexpr double kyr = kiloyear;
static constexpr double Myr = Megayear;
static constexpr double kpc = kiloparsec;
static constexpr double eV = electronvolt;
static constexpr double keV = kiloelectronvolt;
static constexpr double MeV = megaelectronvolt;
static constexpr double GeV = gigaelectronvolt;
static constexpr double TeV = teraelectronvolt;
static constexpr double PeV = petaelectronvolt;
static constexpr double cm = centimeter;
static constexpr double cm2 = cm * cm;
static constexpr double cm3 = cm * cm * cm;
static constexpr double m2 = meter * meter;
static constexpr double pc = parsec;
static constexpr double K = kelvin;
static constexpr double muG = microgauss;

// PHYSICAL CONSTANTS
static constexpr double cLight = 2.99792458e10 * centimeter / second;
static constexpr double cSquared = cLight * cLight;
static constexpr double protonMass = 1.67262158e-24 * gram;
static constexpr double protonMassC = protonMass * cLight;
static constexpr double protonMassC2 = protonMass * cSquared;
static constexpr double neutronMass = 1.67492735e-24 * gram;
static constexpr double neutronMassC2 = neutronMass * cSquared;
static constexpr double electronMass = 9.10938291e-28 * gram;
static constexpr double electronMassC2 = electronMass * cSquared;
static constexpr double sunMass = 1.989e33 * gram;
static constexpr double hPlanck = 6.62607015e-34 * joule * second;
static constexpr double kBoltzmann = 1.3806488e-23 * joule / kelvin;
static constexpr double electronRadius = 2.8179403227e-15 * meter;
static constexpr double elementaryCharge = 4.80320427e-10;
static constexpr double IsH = 19 * eV;   // H  eff. ioniz. potential
static constexpr double IsHe = 44 * eV;  // He eff. ioniz. potential
static constexpr double sigmaTh = 6.6524e-25 * cm2;
static constexpr double barn = 1e-24 * cm2;
static constexpr double mbarn = 1e-3 * barn;

}  // namespace cgs

#endif  // INCLUDE_CGS_H_
