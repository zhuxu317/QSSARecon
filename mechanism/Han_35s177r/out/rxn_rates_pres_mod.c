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
  pres_mod[0] = m + 1.5 * C[1] + 11.0 * C[5] + 0.8999999999999999 * C[8] + 2.8 * C[9] - 0.17000000000000004 * C[10];

  // reaction 5;
  pres_mod[1] = m + 1.5 * C[1] + 11.0 * C[5] + 0.8999999999999999 * C[8] + 2.8 * C[9] - 0.17000000000000004 * C[10] - 0.17000000000000004 * C[34];

  // reaction 6;
  pres_mod[2] = m + 1.5 * C[1] + 11.0 * C[5] + 0.5 * C[8] + 1.0 * C[9] - 0.25 * C[10] - 0.25 * C[34];

  // reaction 7;
  pres_mod[3] = m + 1.5 * C[1] + 11.0 * C[5] + 0.8999999999999999 * C[8] + 2.8 * C[9] - 0.56 * C[10] - 0.62 * C[34];

  // reaction 8;
  thd = m + 0.5109999999999999 * C[1] + 10.372 * C[5] + 0.8999999999999999 * C[8] + 2.8 * C[9] - 0.35 * C[10] - 0.526 * C[34];
  k0 = exp(3.1595048163103488e+01 - 1.37367 * logT);
  kinf = exp(2.2260133056545676e+01 + 0.44 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(3.30000000e-01 * exp(-T / 1.00000000e-30) + 6.70000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[4] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 15;
  thd = m + 1.4700000000000002 * C[1] - 0.19999999999999996 * C[3] + 4.0 * C[5] + 4.13 * C[7] + 0.8700000000000001 * C[8] + 0.07000000000000006 * C[9] - 0.5700000000000001 * C[10] - 0.32999999999999996 * C[34];
  k0 = exp(2.8320561800894382e+01 - 1.17797 * logT - (-2.1501062347737443e+03 / T));
  kinf = exp(5.3697073626347098e+00 + 2.3219 * logT - (-1.7121542474768148e+03 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.70000000e-01 * exp(-T / 1.00000000e-30) + 4.30000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[5] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 21;
  thd = m + 1.0 * C[1] + 11.0 * C[5] + 0.75 * C[8] + 2.6 * C[9] - 0.30000000000000004 * C[10] - 0.30000000000000004 * C[34];
  k0 = exp(4.1606096243564160e+01 - 2.79 * logT - (2.1089931963247504e+03 / T));
  kinf = exp(1.6427049858685642e+01 - (1.1996754426242437e+03 / T));
  Pr = k0 * thd / kinf;
  pres_mod[6] =  Pr / (1.0 + Pr);

  // reaction 25;
  thd = m + 1.0 * C[1] + 11.0 * C[5] + 0.5 * C[8] + 1.0 * C[9] - 0.21399999999999997 * C[10] - 0.44999999999999996 * C[34];
  k0 = exp(1.7715987159492048e+01 + 0.95965 * logT - (7.3671344295381414e+03 / T));
  kinf = exp(3.8436700475959327e+01 - 0.93 * logT - (9.9260059168469834e+03 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(1.48000000e-01 * exp(-T / 5.14000000e+01) + 8.52000000e-01 * exp(-T / 3.57000000e+03) + exp(-3.42000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[7] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 33;
  pres_mod[8] = m - 0.6 * C[3] + 5.5 * C[5] - 0.6 * C[33] - 0.65 * C[34];

  // reaction 49;
  pres_mod[9] = m;

  // reaction 69;
  thd = m;
  k0 = exp(6.4942386233079020e+01 - 5.49 * logT - (9.9989727537515614e+02 / T));
  kinf = exp(2.7051202620675607e+01 - 0.414 * logT - (3.3212491280704732e+01 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.90000000e-01 * exp(-T / 1.00000000e-30) + 3.10000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[10] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 84;
  pres_mod[11] = m + 15.25 * C[5];

  // reaction 87;
  pres_mod[12] = m + 15.25 * C[5];

  // reaction 95;
  pres_mod[13] = m - 1.0 * C[0] - 1.0 * C[1] - 1.0 * C[2] - 1.0 * C[4] - 1.0 * C[5] - 1.0 * C[6] - 1.0 * C[7] - 1.0 * C[8] - 1.0 * C[9] - 1.0 * C[10] - 1.0 * C[11] - 1.0 * C[12] - 1.0 * C[13] - 1.0 * C[14] - 1.0 * C[15] - 1.0 * C[16] - 1.0 * C[17] - 1.0 * C[18] - 1.0 * C[19] - 1.0 * C[20] - 1.0 * C[21] - 1.0 * C[22] - 1.0 * C[23] - 1.0 * C[24] - 1.0 * C[25] - 1.0 * C[26] - 1.0 * C[27] - 1.0 * C[28] - 1.0 * C[29] - 1.0 * C[30] - 1.0 * C[31] - 1.0 * C[32] - 1.0 * C[33] - 1.0 * C[34];

  // reaction 122;
  thd = m - 0.19999999999999996 * C[3] + 9.0 * C[5] + 5.0 * C[9] + 0.8 * C[14] + 5.2 * C[15] + 3.4000000000000004 * C[16] - 0.4 * C[34];
  k0 = exp(4.3691487654050235e+01 - 2.87 * logT - (7.8049354509656121e+02 / T));
  kinf = exp(2.7893385380396040e+01 - 0.75 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(3.80000000e-02 * exp(-T / 1.00000000e+01) + 9.62000000e-01 * exp(-T / 7.96200000e+03), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[14] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 125;
  pres_mod[15] = m + 0.25 * C[1] + 3.0999999999999996 * C[5];

  // reaction 128;
  thd = m + 0.3999999999999999 * C[3] + 11.0 * C[5] + 2.0 * C[14] + 2.5 * C[16] + 0.7 * C[33];
  k0 = exp(2.7302517048956513e+01 - (2.8889835218564527e+04 / T));
  kinf = exp(2.5853164551869483e+01 - (2.9012117572825304e+04 / T));
  Pr = k0 * thd / kinf;
  pres_mod[16] =  Pr / (1.0 + Pr);

  // reaction 137;
  thd = m;
  k0 = exp(4.0365366298828434e+01 - 2.5 * logT);
  kinf = exp(2.5423746202738826e+01 - 0.3 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.50000000e-01 * exp(-T / 1.00000000e-30) + 7.50000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[17] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 146;
  thd = m;
  k0 = exp(3.5670178506401783e+01 - (1.5851416293063623e+04 / T));
  kinf = exp(3.3152482033790797e+01 - (1.6253991944950953e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(-1.49000000e-01 * exp(-T / 1.00000000e-30) + 1.14900000e+00 * exp(-T / 3.12500000e+03) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[18] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 152;
  pres_mod[19] = m + 9.0 * C[5];

  // reaction 153;
  pres_mod[20] = m + 9.0 * C[5];

  // reaction 165;
  thd = m;
  k0 = exp(4.4826845844638555e+01 - 3.0 * logT);
  kinf = 30000000000.0;
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.00000000e-01 * exp(-T / 1.00000000e-30) + 4.00000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[21] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 171;
  thd = m;
  k0 = exp(3.3152482033790797e+01 - 1.5 * logT);
  kinf = exp(2.1976028805441779e+01 + 0.24 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.90000000e-01 * exp(-T / 1.00000000e-30) + 7.10000000e-01 * exp(-T / 1.70000000e+03) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[22] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 176;
  pres_mod[23] = m - 1.0 * C[0] - 1.0 * C[1] - 1.0 * C[2] - 1.0 * C[3] - 1.0 * C[4] - 1.0 * C[5] - 1.0 * C[6] - 1.0 * C[7] - 1.0 * C[8] - 1.0 * C[9] - 1.0 * C[10] - 1.0 * C[11] - 1.0 * C[12] - 1.0 * C[13] - 1.0 * C[14] - 1.0 * C[16] - 1.0 * C[17] - 1.0 * C[18] - 1.0 * C[19] - 1.0 * C[20] - 1.0 * C[21] - 1.0 * C[22] - 1.0 * C[23] - 1.0 * C[24] - 1.0 * C[25] - 1.0 * C[26] - 1.0 * C[27] - 1.0 * C[28] - 1.0 * C[29] - 1.0 * C[30] - 1.0 * C[31] - 1.0 * C[32] - 1.0 * C[33] - 1.0 * C[34];

} // end get_rxn_pres_mod

