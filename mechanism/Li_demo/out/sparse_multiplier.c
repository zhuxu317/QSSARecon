#include "sparse_multiplier.h"

void sparse_multiplier(const double * A, const double * Vm, double* w) {
  w[0] =  A[0] * Vm[0] +  A[6] * Vm[1] +  A[12] * Vm[2] +  A[18] * Vm[3] +  A[24] * Vm[4] +  A[30] * Vm[5];
  w[1] =  A[1] * Vm[0] +  A[7] * Vm[1] +  A[13] * Vm[2] +  A[19] * Vm[3] +  A[25] * Vm[4] +  A[31] * Vm[5];
  w[2] =  A[2] * Vm[0] +  A[8] * Vm[1] +  A[14] * Vm[2] +  A[20] * Vm[3] +  A[26] * Vm[4] +  A[32] * Vm[5];
  w[3] =  A[3] * Vm[0] +  A[9] * Vm[1] +  A[15] * Vm[2] +  A[21] * Vm[3] +  A[27] * Vm[4] +  A[33] * Vm[5];
  w[4] =  A[4] * Vm[0] +  A[10] * Vm[1] +  A[16] * Vm[2] +  A[22] * Vm[3] +  A[28] * Vm[4] +  A[34] * Vm[5];
  w[5] =  A[5] * Vm[0] +  A[11] * Vm[1] +  A[17] * Vm[2] +  A[23] * Vm[3] +  A[29] * Vm[4] +  A[35] * Vm[5];
}
