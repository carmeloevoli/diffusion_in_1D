#ifndef UNITS_H
#define UNITS_H

// TIME UNITS
static const double sec = 1.;
static const double year = 3.154e+7 * sec;
static const double kyr = 1e3 * year;
static const double Myr = 1e6 * year;
static const double Gyr = 1e9 * year;

// LENGTH UNITS
static const double meter = 1.;
static const double cm = 1e-2 * meter;
static const double km = 1e3 * meter;
static const double parsec = 3.086e16 * meter;
static const double kpc = 1e3 * parsec;

// COMBINED UNITS
static const double cm2 = cm * cm;

#define pow2(A) ((A)*(A))

#endif
