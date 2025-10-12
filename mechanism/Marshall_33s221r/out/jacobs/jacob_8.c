#include "header.h"

void eval_jacob_8 (const double mw_avg, const double rho, const double cp_avg, const double* spec_rates, const double* h, const double* cp, double* jac) {
  double rho_inv = 1.0 / rho;
  double working_temp = (1.0 / cp_avg);
  double j_temp = 1.0 / (rho * cp_avg * cp_avg);
}

