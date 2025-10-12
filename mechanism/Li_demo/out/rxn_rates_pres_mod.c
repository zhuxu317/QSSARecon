#include <math.h>
#include "header.h"
#include "rates.h"

void get_rxn_pres_mod (const double T, const double pres, const double * __restrict__ C, double * __restrict__ pres_mod) {
  // third body variable declaration
  double thd;

  double logT = log(T);
  double m = pres / (8.31446210e+03 * T);

  // reaction 3;
  pres_mod[0] = m + 1.5 * C[0] + 11.0 * C[2];

  // reaction 4;
  pres_mod[1] = m + 1.5 * C[0] + 11.0 * C[2];

  // reaction 5;
  pres_mod[2] = m + 1.5 * C[0] + 11.0 * C[2];

} // end get_rxn_pres_mod

