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

  // reaction 6;
  pres_mod[0] = m + 1.5 * C[2] + 11.0 * C[6] + 0.8999999999999999 * C[10] + 2.8 * C[11] - 1.0 * C[0] - 1.0 * C[3];

  // reaction 9;
  pres_mod[1] = m + 1.5 * C[2] + 11.0 * C[6] - 1.0 * C[0] - 1.0 * C[3] + 0.8999999999999999 * C[10] + 2.8 * C[11];

  // reaction 12;
  pres_mod[2] = m + 1.5 * C[2] + 11.0 * C[6] - 0.25 * C[0] - 0.25 * C[3] + 0.8999999999999999 * C[10] + 2.8 * C[11];

  // reaction 13;
  pres_mod[3] = m + 2.0 * C[2] - 1.0 * C[6] + 0.10000000000000009 * C[3] + 1.0 * C[35] + 0.5 * C[9] + 0.8999999999999999 * C[10] + 2.8 * C[11];

  // reaction 15;
  thd = m + 1.0 * C[2] + 13.0 * C[6] - 0.21999999999999997 * C[9] + 0.8999999999999999 * C[10] + 2.8 * C[11] - 0.32999999999999996 * C[0] - 0.19999999999999996 * C[3];
  k0 = exp(3.4087162630776540e+01 - 1.72 * logT - (2.6408962763808859e+02 / T));
  kinf = exp(2.2260313685392592e+01 + 0.44 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.00000000e-01 * exp(-T / 1.00000000e-30) + 5.00000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[4] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 27;
  thd = m + 6.5 * C[6] + 0.6000000000000001 * C[11] + 0.5 * C[35] + 0.19999999999999996 * C[9] - 0.35 * C[3] + 6.7 * C[8] + 2.7 * C[2] + 1.7999999999999998 * C[10];
  k0 = exp(4.9266569663351575e+01 - 2.3 * logT - (2.4531450567319323e+04 / T));
  kinf = exp(2.8324168296488494e+01 + 0.9 * logT - (2.4531450567319323e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(5.70000000e-01 * exp(-T / 1.00000000e-30) + 4.30000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[5] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 33;
  thd = m + 1.5 * C[2] + 11.0 * C[6] + 0.8999999999999999 * C[10] + 2.8 * C[11] - 0.13 * C[0];
  k0 = exp(4.1884786604823979e+01 - 2.79 * logT - (2.1089931963247509e+03 / T));
  kinf = exp(1.6705882315860439e+01 - (1.1996754426242439e+03 / T));
  Pr = k0 * thd / kinf;
  pres_mod[6] =  Pr / (1.0 + Pr);

  // reaction 38;
  thd = m;
  k0 = exp(5.4751216608091106e+01 - 3.148 * logT - (1.8677497369312685e+04 / T));
  kinf = exp(2.7432570177204710e+01 + 0.413 * logT - (1.7781263324298518e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.10000000e-01 * exp(-T / 1.00000000e-30) + 3.90000000e-01 * exp(-T / 1.00000000e+30), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[7] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 48;
  thd = m + 0.5 * C[35] + 0.5 * C[3] + 0.5 * C[9] + 0.5 * C[10] + 1.0 * C[2] + 2.0 * C[11] + 14.0 * C[6];
  k0 = exp(4.3452057532622490e+01 - 2.36 * logT - (9.7539048256651513e+03 / T));
  kinf = exp(3.8436700475959327e+01 - 0.93 * logT - (9.9255026972821252e+03 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(8.97000000e-01 * exp(-T / 1.39000000e+02) + 1.03000000e-01 * exp(-T / 1.09000000e+04) + exp(-4.55000000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[8] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 55;
  pres_mod[9] = m + 15.25 * C[6] + 0.875 * C[10] + 2.75 * C[11];

  // reaction 56;
  pres_mod[10] = m + 15.25 * C[6] + 0.875 * C[10] + 2.75 * C[11];

  // reaction 77;
  pres_mod[11] = m;

  // reaction 97;
  pres_mod[12] = m + 6.0 * C[6];

  // reaction 109;
  thd = m - 0.41000000000000003 * C[0] - 0.31000000000000005 * C[9] + 3.87 * C[17];
  k0 = exp(6.4942386233079020e+01 - 5.49 * logT - (9.9989727537515637e+02 / T));
  kinf = exp(2.7051202620675607e+01 - 0.414 * logT - (3.3212491280704739e+01 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.90000000e-01 * exp(-T / 1.00000000e-30) + 3.10000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[13] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 121;
  thd = m + 1.0 * C[35];
  k0 = exp(8.6541120807379329e+01 - 6.88 * logT - (2.7425466284824372e+04 / T));
  kinf = exp(2.5575296100866030e+01 + 0.819 * logT - (2.4204861069725728e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(1.58000000e-01 * exp(-T / 8.00000000e+04) + 8.42000000e-01 * exp(-T / 2.80000000e+01) + exp(-7.29800000e+03 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[14] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 151;
  pres_mod[15] = m + 9.0 * C[6];

  // reaction 152;
  pres_mod[16] = m + 9.0 * C[6];

  // reaction 164;
  thd = m;
  k0 = exp(7.9974292115367788e+01 - 5.96 * logT - (3.3605002541294889e+04 / T));
  kinf = exp(4.6388174096502127e+01 - 1.31 * logT - (3.2246309716175150e+04 / T));
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.50000000e-01 * exp(-T / 1.00000000e-30) + 3.50000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[17] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 178;
  pres_mod[18] = m;

  // reaction 206;
  pres_mod[19] = m - 0.30000000000000004 * C[0] + 6.0 * C[6] + 1.0 * C[11];

  // reaction 218;
  thd = m + 5.4 * C[6] + 0.5 * C[11] - 0.55 * C[9] - 0.6 * C[35] - 0.65 * C[0] - 0.65 * C[3] - 0.25 * C[10];
  k0 = exp(3.3440963786413143e+01 - 1.6 * logT);
  kinf = exp(2.6410241193286232e+01 - 0.4 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.00000000e-01 * exp(-T / 1.00000000e-30) + 8.00000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[20] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 219;
  thd = m;
  k0 = exp(4.4826845844638555e+01 - 3.0 * logT);
  kinf = 30000000000.0;
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(6.00000000e-01 * exp(-T / 1.00000000e-30) + 4.00000000e-01 * exp(-T / 1.00000000e+30) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[21] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

  // reaction 223;
  thd = m + 0.7 * C[35] + 0.3999999999999999 * C[9] + 2.0 * C[11] + 11.0 * C[6];
  k0 = exp(2.6714730384054395e+01 - (2.8482227371028610e+04 / T));
  kinf = exp(2.7893385380396040e+01 - (3.1486448173237812e+04 / T));
  Pr = k0 * thd / kinf;
  pres_mod[22] =  Pr / (1.0 + Pr);

  // reaction 231;
  pres_mod[23] = m + 9.0 * C[6];

  // reaction 257;
  thd = m;
  k0 = exp(3.3152482033790797e+01 - 1.5 * logT);
  kinf = exp(2.1976028805441779e+01 + 0.2 * logT);
  Pr = k0 * thd / kinf;
  logFcent = log10( fmax(2.90000000e-01 * exp(-T / 1.00000000e-30) + 7.10000000e-01 * exp(-T / 1.70000000e+03) + exp(-1.00000000e+30 / T), 1.0e-300));
  A = log10(fmax(Pr, 1.0e-300)) - 0.67 * logFcent - 0.4;
  B = 0.806 - 1.1762 * logFcent - 0.14 * log10(fmax(Pr, 1.0e-300));
  pres_mod[24] = pow(10.0, logFcent / (1.0 + A * A / (B * B))) * Pr / (1.0 + Pr);

} // end get_rxn_pres_mod

