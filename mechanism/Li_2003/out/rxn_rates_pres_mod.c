#include <math.h>
#include "header.h"
#include "rates.h"

void get_rxn_pres_mod (const double T, const double pres, const double * __restrict__ C, double * __restrict__ pres_mod) {
  // third body variable declaration
  double thd;

  // pressure dependence variable declarations
  double k0;
  double kinf;
  double Pr;

  // troe variable declarations
  double logFcent;
  double A;
  double B;

  double logT = log(T);
  double m = pres / (8.31446210e+03 * T);

  // reaction 4;
  pres_mod[0] = m + 1.5 * C[0] + 11.0 * C[2] + 0.8999999999999999 * C[9] + 2.8 * C[10] - 1.0 * C[11];

  // reaction 5;
  pres_mod[1] = m - 1.0 * C[0] - 1.0 * C[1] - 1.0 * C[2] - 1.0 * C[3] - 1.0 * C[4] - 1.0 * C[5] - 1.0 * C[6] - 1.0 * C[7] - 1.0 * C[8] - 1.0 * C[9] - 1.0 * C[10];

  // reaction 6;
  pres_mod[2] = m + 1.5 * C[0] + 11.0 * C[2] + 0.8999999999999999 * C[9] + 2.8 * C[10] - 1.0 * C[11];

  // reaction 7;
  pres_mod[3] = m - 1.0 * C[0] - 1.0 * C[1] - 1.0 * C[2] - 1.0 * C[3] - 1.0 * C[4] - 1.0 * C[5] - 1.0 * C[6] - 1.0 * C[7] - 1.0 * C[8] - 1.0 * C[9] - 1.0 * C[10];

  // reaction 8;
  pres_mod[4] = m + 1.5 * C[0] + 11.0 * C[2] + 0.8999999999999999 * C[9] + 2.8 * C[10] - 0.25 * C[11];

  // reaction 9;
  pres_mod[5] = m + 1.5 * C[0] + 11.0 * C[2] + 0.8999999999999999 * C[9] + 2.8 * C[10] - 0.62 * C[11];

  // reaction 10;
  thd = m + 1.0 * C[0] - 0.21999999999999997 * C[1] + 10.0 * C[2] + 0.8999999999999999 * C[9] + 2.8 * C[10];
  k0 = exp(3.4087162630776540e+01 - 1.72 * logT - (2.6408962763808853e+02 / T));
  kinf = exp(2.1111923826738195e+01 + 0.6 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.00000000e-01 * exp(-T / 1.00000000e-30) + 8.00000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[6] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 17;
  thd = m + 1.5 * C[0] + 11.0 * C[2] + 0.8999999999999999 * C[9] + 2.8 * C[10] - 0.36 * C[11];
  k0 = exp(3.2420178138029655e+01 - (2.2896490201091900e+04 / T));
  kinf = exp(3.3318335397877441e+01 - (2.4370923526129245e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.00000000e-01 * exp(-T / 1.00000000e-30) + 5.00000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[7] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

} // end get_rxn_pres_mod

