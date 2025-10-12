#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 0
  sp_rates[0] = 2.0 * (fwd_rates[0] - rev_rates[0]) * pres_mod[0];
  //sp 1
  sp_rates[1] = -(fwd_rates[0] - rev_rates[0]) * pres_mod[0];

  //rxn 1
  //sp 0
  sp_rates[0] += (fwd_rates[1] - rev_rates[1]);
  //sp 1
  sp_rates[1] -= (fwd_rates[1] - rev_rates[1]);
  //sp 12
  sp_rates[12] = -(fwd_rates[1] - rev_rates[1]);
  //sp 14
  sp_rates[14] = (fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 0
  sp_rates[0] += (fwd_rates[2] - rev_rates[2]);
  //sp 1
  sp_rates[1] -= (fwd_rates[2] - rev_rates[2]);
  //sp 2
  sp_rates[2] = (fwd_rates[2] - rev_rates[2]);
  //sp 14
  sp_rates[14] -= (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 0
  sp_rates[0] -= (fwd_rates[3] - rev_rates[3]);
  //sp 12
  sp_rates[12] += (fwd_rates[3] - rev_rates[3]);
  //sp 13
  sp_rates[13] = -(fwd_rates[3] - rev_rates[3]);
  //sp 14
  sp_rates[14] += (fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 0
  sp_rates[0] -= (fwd_rates[4] - rev_rates[4]) * pres_mod[1];
  //sp 2
  sp_rates[2] += (fwd_rates[4] - rev_rates[4]) * pres_mod[1];
  //sp 14
  sp_rates[14] -= (fwd_rates[4] - rev_rates[4]) * pres_mod[1];

  //rxn 5
  //sp 2
  sp_rates[2] -= (fwd_rates[5] - rev_rates[5]);
  //sp 12
  sp_rates[12] -= (fwd_rates[5] - rev_rates[5]);
  //sp 14
  sp_rates[14] += 2.0 * (fwd_rates[5] - rev_rates[5]);

  //rxn 6
  //sp 0
  sp_rates[0] -= (fwd_rates[6] - rev_rates[6]) * pres_mod[2];
  //sp 12
  sp_rates[12] -= (fwd_rates[6] - rev_rates[6]) * pres_mod[2];
  //sp 14
  sp_rates[14] += (fwd_rates[6] - rev_rates[6]) * pres_mod[2];

  //rxn 7
  //sp 3
  sp_rates[3] = -(fwd_rates[7] - rev_rates[7]) * pres_mod[3];
  //sp 14
  sp_rates[14] += 2.0 * (fwd_rates[7] - rev_rates[7]) * pres_mod[3];

  //rxn 8
  //sp 0
  sp_rates[0] -= (fwd_rates[8] - rev_rates[8]);
  //sp 2
  sp_rates[2] += (fwd_rates[8] - rev_rates[8]);
  //sp 3
  sp_rates[3] -= (fwd_rates[8] - rev_rates[8]);
  //sp 14
  sp_rates[14] += (fwd_rates[8] - rev_rates[8]);

  //rxn 9
  //sp 3
  sp_rates[3] -= (fwd_rates[9] - rev_rates[9]);
  //sp 12
  sp_rates[12] -= (fwd_rates[9] - rev_rates[9]);
  //sp 4
  sp_rates[4] = (fwd_rates[9] - rev_rates[9]);
  //sp 14
  sp_rates[14] += (fwd_rates[9] - rev_rates[9]);

  //rxn 10
  //sp 2
  sp_rates[2] += (fwd_rates[10] - rev_rates[10]);
  //sp 3
  sp_rates[3] -= (fwd_rates[10] - rev_rates[10]);
  //sp 4
  sp_rates[4] += (fwd_rates[10] - rev_rates[10]);
  //sp 14
  sp_rates[14] -= (fwd_rates[10] - rev_rates[10]);

  //rxn 11
  //sp 2
  sp_rates[2] += (fwd_rates[11] - rev_rates[11]);
  //sp 3
  sp_rates[3] -= (fwd_rates[11] - rev_rates[11]);
  //sp 4
  sp_rates[4] += (fwd_rates[11] - rev_rates[11]);
  //sp 14
  sp_rates[14] -= (fwd_rates[11] - rev_rates[11]);

  //rxn 12
  //sp 0
  sp_rates[0] -= (fwd_rates[12] - rev_rates[12]);
  //sp 4
  sp_rates[4] -= (fwd_rates[12] - rev_rates[12]);
  //sp 14
  sp_rates[14] += 2.0 * (fwd_rates[12] - rev_rates[12]);

  //rxn 13
  //sp 0
  sp_rates[0] -= (fwd_rates[13] - rev_rates[13]);
  //sp 1
  sp_rates[1] += (fwd_rates[13] - rev_rates[13]);
  //sp 4
  sp_rates[4] -= (fwd_rates[13] - rev_rates[13]);
  //sp 13
  sp_rates[13] += (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 13
  sp_rates[13] += (fwd_rates[14] - rev_rates[14]);
  //sp 4
  sp_rates[4] -= (fwd_rates[14] - rev_rates[14]);
  //sp 12
  sp_rates[12] -= (fwd_rates[14] - rev_rates[14]);
  //sp 14
  sp_rates[14] += (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 2
  sp_rates[2] += (fwd_rates[15] - rev_rates[15]);
  //sp 4
  sp_rates[4] -= (fwd_rates[15] - rev_rates[15]);
  //sp 13
  sp_rates[13] += (fwd_rates[15] - rev_rates[15]);
  //sp 14
  sp_rates[14] -= (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 2
  sp_rates[2] += (fwd_rates[16] - rev_rates[16]);
  //sp 4
  sp_rates[4] -= (fwd_rates[16] - rev_rates[16]);
  //sp 13
  sp_rates[13] += (fwd_rates[16] - rev_rates[16]);
  //sp 14
  sp_rates[14] -= (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 3
  sp_rates[3] += (fwd_rates[17] - rev_rates[17]);
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[17] - rev_rates[17]);
  //sp 13
  sp_rates[13] += (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 3
  sp_rates[3] += (fwd_rates[18] - rev_rates[18]);
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[18] - rev_rates[18]);
  //sp 13
  sp_rates[13] += (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 0
  sp_rates[0] -= (fwd_rates[19] - rev_rates[19]) * pres_mod[4];
  //sp 4
  sp_rates[4] += (fwd_rates[19] - rev_rates[19]) * pres_mod[4];
  //sp 13
  sp_rates[13] -= (fwd_rates[19] - rev_rates[19]) * pres_mod[4];

  //rxn 20
  //sp 12
  sp_rates[12] -= (fwd_rates[20] - rev_rates[20]) * pres_mod[5];
  //sp 4
  sp_rates[4] += (fwd_rates[20] - rev_rates[20]) * pres_mod[5];
  //sp 14
  sp_rates[14] -= (fwd_rates[20] - rev_rates[20]) * pres_mod[5];

  //rxn 21
  //sp 0
  sp_rates[0] -= (fwd_rates[21] - rev_rates[21]);
  //sp 9
  sp_rates[9] = -(fwd_rates[21] - rev_rates[21]);
  //sp 1
  sp_rates[1] += (fwd_rates[21] - rev_rates[21]);
  //sp 8
  sp_rates[8] = (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 8
  sp_rates[8] += (fwd_rates[22] - rev_rates[22]);
  //sp 9
  sp_rates[9] -= (fwd_rates[22] - rev_rates[22]);
  //sp 12
  sp_rates[12] -= (fwd_rates[22] - rev_rates[22]);
  //sp 14
  sp_rates[14] += (fwd_rates[22] - rev_rates[22]);

  //rxn 23
  //sp 8
  sp_rates[8] += (fwd_rates[23] - rev_rates[23]);
  //sp 9
  sp_rates[9] -= (fwd_rates[23] - rev_rates[23]);
  //sp 2
  sp_rates[2] += (fwd_rates[23] - rev_rates[23]);
  //sp 14
  sp_rates[14] -= (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 0
  sp_rates[0] += 2.0 * (fwd_rates[24] - rev_rates[24]);
  //sp 9
  sp_rates[9] -= (fwd_rates[24] - rev_rates[24]);
  //sp 5
  sp_rates[5] = -(fwd_rates[24] - rev_rates[24]);
  //sp 6
  sp_rates[6] = (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 9
  sp_rates[9] -= (fwd_rates[25] - rev_rates[25]);
  //sp 10
  sp_rates[10] = (fwd_rates[25] - rev_rates[25]);
  //sp 11
  sp_rates[11] = -(fwd_rates[25] - rev_rates[25]);
  //sp 14
  sp_rates[14] += (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 0
  sp_rates[0] -= (fwd_rates[26] - rev_rates[26]);
  //sp 8
  sp_rates[8] -= (fwd_rates[26] - rev_rates[26]);
  //sp 5
  sp_rates[5] += (fwd_rates[26] - rev_rates[26]);
  //sp 1
  sp_rates[1] += (fwd_rates[26] - rev_rates[26]);

  //rxn 27
  //sp 8
  sp_rates[8] -= (fwd_rates[27] - rev_rates[27]);
  //sp 0
  sp_rates[0] += (fwd_rates[27] - rev_rates[27]);
  //sp 11
  sp_rates[11] += (fwd_rates[27] - rev_rates[27]);
  //sp 12
  sp_rates[12] -= (fwd_rates[27] - rev_rates[27]);

  //rxn 28
  //sp 8
  sp_rates[8] -= (fwd_rates[28] - rev_rates[28]);
  //sp 2
  sp_rates[2] += (fwd_rates[28] - rev_rates[28]);
  //sp 5
  sp_rates[5] += (fwd_rates[28] - rev_rates[28]);
  //sp 14
  sp_rates[14] -= (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 8
  sp_rates[8] -= (fwd_rates[29] - rev_rates[29]);
  //sp 1
  sp_rates[1] += (fwd_rates[29] - rev_rates[29]);
  //sp 11
  sp_rates[11] += (fwd_rates[29] - rev_rates[29]);
  //sp 14
  sp_rates[14] -= (fwd_rates[29] - rev_rates[29]);

  //rxn 30
  //sp 8
  sp_rates[8] -= (fwd_rates[30] - rev_rates[30]);
  //sp 11
  sp_rates[11] += (fwd_rates[30] - rev_rates[30]);
  //sp 13
  sp_rates[13] -= (fwd_rates[30] - rev_rates[30]);
  //sp 14
  sp_rates[14] += (fwd_rates[30] - rev_rates[30]);

  //rxn 31
  //sp 8
  sp_rates[8] -= (fwd_rates[31] - rev_rates[31]);
  //sp 11
  sp_rates[11] -= (fwd_rates[31] - rev_rates[31]);
  //sp 6
  sp_rates[6] += (fwd_rates[31] - rev_rates[31]);
  //sp 14
  sp_rates[14] += (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 0
  sp_rates[0] += (fwd_rates[32] - rev_rates[32]);
  //sp 11
  sp_rates[11] += (fwd_rates[32] - rev_rates[32]);
  //sp 5
  sp_rates[5] -= (fwd_rates[32] - rev_rates[32]);
  //sp 14
  sp_rates[14] -= (fwd_rates[32] - rev_rates[32]);

  //rxn 33
  //sp 13
  sp_rates[13] -= (fwd_rates[33] - rev_rates[33]);
  //sp 11
  sp_rates[11] += (fwd_rates[33] - rev_rates[33]);
  //sp 12
  sp_rates[12] += (fwd_rates[33] - rev_rates[33]);
  //sp 5
  sp_rates[5] -= (fwd_rates[33] - rev_rates[33]);

  //rxn 34
  //sp 11
  sp_rates[11] -= (fwd_rates[34] - rev_rates[34]);
  //sp 12
  sp_rates[12] += (fwd_rates[34] - rev_rates[34]);
  //sp 5
  sp_rates[5] -= (fwd_rates[34] - rev_rates[34]);
  //sp 6
  sp_rates[6] += (fwd_rates[34] - rev_rates[34]);

  //rxn 35
  //sp 0
  sp_rates[0] += (fwd_rates[35] - rev_rates[35]);
  //sp 10
  sp_rates[10] -= (fwd_rates[35] - rev_rates[35]);
  //sp 6
  sp_rates[6] += (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 0
  sp_rates[0] -= (fwd_rates[36] - rev_rates[36]);
  //sp 1
  sp_rates[1] += (fwd_rates[36] - rev_rates[36]);
  //sp 10
  sp_rates[10] -= (fwd_rates[36] - rev_rates[36]);
  //sp 6
  sp_rates[6] += (fwd_rates[36] - rev_rates[36]);

  //rxn 37
  //sp 0
  sp_rates[0] += (fwd_rates[37] - rev_rates[37]);
  //sp 10
  sp_rates[10] -= (fwd_rates[37] - rev_rates[37]);
  //sp 12
  sp_rates[12] -= (fwd_rates[37] - rev_rates[37]);
  //sp 7
  sp_rates[7] = (fwd_rates[37] - rev_rates[37]);

  //rxn 38
  //sp 8
  sp_rates[8] += (fwd_rates[38] - rev_rates[38]);
  //sp 10
  sp_rates[10] -= (fwd_rates[38] - rev_rates[38]);
  //sp 11
  sp_rates[11] += (fwd_rates[38] - rev_rates[38]);
  //sp 12
  sp_rates[12] -= (fwd_rates[38] - rev_rates[38]);

  //rxn 39
  //sp 10
  sp_rates[10] -= (fwd_rates[39] - rev_rates[39]);
  //sp 2
  sp_rates[2] += (fwd_rates[39] - rev_rates[39]);
  //sp 14
  sp_rates[14] -= (fwd_rates[39] - rev_rates[39]);
  //sp 6
  sp_rates[6] += (fwd_rates[39] - rev_rates[39]);

  //rxn 40
  //sp 10
  sp_rates[10] -= (fwd_rates[40] - rev_rates[40]);
  //sp 4
  sp_rates[4] += (fwd_rates[40] - rev_rates[40]);
  //sp 13
  sp_rates[13] -= (fwd_rates[40] - rev_rates[40]);
  //sp 6
  sp_rates[6] += (fwd_rates[40] - rev_rates[40]);

  //rxn 41
  //sp 12
  sp_rates[12] += (fwd_rates[41] - rev_rates[41]) * pres_mod[6];
  //sp 6
  sp_rates[6] += (fwd_rates[41] - rev_rates[41]) * pres_mod[6];
  //sp 7
  sp_rates[7] -= (fwd_rates[41] - rev_rates[41]) * pres_mod[6];

  //rxn 42
  //sp 0
  sp_rates[0] -= (fwd_rates[42] - rev_rates[42]);
  //sp 14
  sp_rates[14] += (fwd_rates[42] - rev_rates[42]);
  //sp 6
  sp_rates[6] += (fwd_rates[42] - rev_rates[42]);
  //sp 7
  sp_rates[7] -= (fwd_rates[42] - rev_rates[42]);

  //rxn 43
  //sp 0
  sp_rates[0] -= (fwd_rates[43] - rev_rates[43]);
  //sp 14
  sp_rates[14] += (fwd_rates[43] - rev_rates[43]);
  //sp 6
  sp_rates[6] += (fwd_rates[43] - rev_rates[43]);
  //sp 7
  sp_rates[7] -= (fwd_rates[43] - rev_rates[43]);

  //rxn 44
  //sp 8
  sp_rates[8] -= (fwd_rates[44] - rev_rates[44]);
  //sp 0
  sp_rates[0] += (fwd_rates[44] - rev_rates[44]);
  //sp 11
  sp_rates[11] -= (fwd_rates[44] - rev_rates[44]);
  //sp 7
  sp_rates[7] += (fwd_rates[44] - rev_rates[44]);

  //rxn 45
  //sp 11
  sp_rates[11] += 2.0 * (fwd_rates[45] - rev_rates[45]);
  //sp 12
  sp_rates[12] -= (fwd_rates[45] - rev_rates[45]);
  //sp 7
  sp_rates[7] -= (fwd_rates[45] - rev_rates[45]);

  //rxn 46
  //sp 12
  sp_rates[12] -= (fwd_rates[46] - rev_rates[46]);
  //sp 13
  sp_rates[13] += (fwd_rates[46] - rev_rates[46]);
  //sp 6
  sp_rates[6] += (fwd_rates[46] - rev_rates[46]);
  //sp 7
  sp_rates[7] -= (fwd_rates[46] - rev_rates[46]);

  //sp 15
  (*dy_N) = 0.0;
} // end eval_spec_rates

