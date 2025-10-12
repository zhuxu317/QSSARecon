#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 1
  sp_rates[1] = -2.0 * (fwd_rates[0] - rev_rates[0]) * pres_mod[0];
  //sp 2
  sp_rates[2] = (fwd_rates[0] - rev_rates[0]) * pres_mod[0];

  //rxn 1
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[1] - rev_rates[1]) * pres_mod[1];
  //sp 2
  sp_rates[2] += (fwd_rates[1] - rev_rates[1]) * pres_mod[1];

  //rxn 2
  //sp 1
  sp_rates[1] += (fwd_rates[2] - rev_rates[2]);
  //sp 2
  sp_rates[2] -= (fwd_rates[2] - rev_rates[2]);
  //sp 3
  sp_rates[3] = -(fwd_rates[2] - rev_rates[2]);
  //sp 4
  sp_rates[4] = (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 1
  sp_rates[1] -= (fwd_rates[3] - rev_rates[3]) * pres_mod[2];
  //sp 3
  sp_rates[3] -= (fwd_rates[3] - rev_rates[3]) * pres_mod[2];
  //sp 4
  sp_rates[4] += (fwd_rates[3] - rev_rates[3]) * pres_mod[2];

  //rxn 4
  //sp 3
  sp_rates[3] -= 2.0 * (fwd_rates[4] - rev_rates[4]) * pres_mod[3];
  //sp 5
  sp_rates[5] = (fwd_rates[4] - rev_rates[4]) * pres_mod[3];

  //rxn 5
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[5] - rev_rates[5]) * pres_mod[4];
  //sp 6
  sp_rates[6] = (fwd_rates[5] - rev_rates[5]) * pres_mod[4];

  //rxn 6
  //sp 1
  sp_rates[1] -= (fwd_rates[6] - rev_rates[6]) * pres_mod[5];
  //sp 4
  sp_rates[4] -= (fwd_rates[6] - rev_rates[6]) * pres_mod[5];
  //sp 7
  sp_rates[7] = (fwd_rates[6] - rev_rates[6]) * pres_mod[5];

  //rxn 7
  //sp 3
  sp_rates[3] += (fwd_rates[7] - rev_rates[7]);
  //sp 4
  sp_rates[4] -= 2.0 * (fwd_rates[7] - rev_rates[7]);
  //sp 7
  sp_rates[7] += (fwd_rates[7] - rev_rates[7]);

  //rxn 8
  //sp 1
  sp_rates[1] += (fwd_rates[8] - rev_rates[8]);
  //sp 2
  sp_rates[2] -= (fwd_rates[8] - rev_rates[8]);
  //sp 4
  sp_rates[4] -= (fwd_rates[8] - rev_rates[8]);
  //sp 7
  sp_rates[7] += (fwd_rates[8] - rev_rates[8]);

  //rxn 9
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[9] - rev_rates[9]) * pres_mod[6];
  //sp 2
  sp_rates[2] += (fwd_rates[9] - rev_rates[9]) * pres_mod[6];

  //rxn 10
  //sp 8
  sp_rates[8] = (fwd_rates[10] - rev_rates[10]) * pres_mod[7];
  //sp 1
  sp_rates[1] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[7];
  //sp 5
  sp_rates[5] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[7];

  //rxn 11
  //sp 1
  sp_rates[1] -= (fwd_rates[11] - rev_rates[11]);
  //sp 3
  sp_rates[3] += (fwd_rates[11] - rev_rates[11]);
  //sp 4
  sp_rates[4] += (fwd_rates[11] - rev_rates[11]);
  //sp 5
  sp_rates[5] -= (fwd_rates[11] - rev_rates[11]);

  //rxn 12
  //sp 8
  sp_rates[8] += (fwd_rates[12] - rev_rates[12]);
  //sp 1
  sp_rates[1] += (fwd_rates[12] - rev_rates[12]);
  //sp 2
  sp_rates[2] -= (fwd_rates[12] - rev_rates[12]);
  //sp 5
  sp_rates[5] -= (fwd_rates[12] - rev_rates[12]);

  //rxn 13
  //sp 8
  sp_rates[8] -= (fwd_rates[13] - rev_rates[13]);
  //sp 4
  sp_rates[4] -= (fwd_rates[13] - rev_rates[13]);
  //sp 5
  sp_rates[5] += (fwd_rates[13] - rev_rates[13]);
  //sp 7
  sp_rates[7] += (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 8
  sp_rates[8] -= (fwd_rates[14] - rev_rates[14]);
  //sp 1
  sp_rates[1] -= (fwd_rates[14] - rev_rates[14]);
  //sp 4
  sp_rates[4] += 2.0 * (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 8
  sp_rates[8] -= (fwd_rates[15] - rev_rates[15]);
  //sp 3
  sp_rates[3] -= (fwd_rates[15] - rev_rates[15]);
  //sp 4
  sp_rates[4] += (fwd_rates[15] - rev_rates[15]);
  //sp 5
  sp_rates[5] += (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 8
  sp_rates[8] -= (fwd_rates[16] - rev_rates[16]);
  //sp 1
  sp_rates[1] -= (fwd_rates[16] - rev_rates[16]);
  //sp 3
  sp_rates[3] += (fwd_rates[16] - rev_rates[16]);
  //sp 7
  sp_rates[7] += (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 8
  sp_rates[8] -= (fwd_rates[17] - rev_rates[17]);
  //sp 4
  sp_rates[4] -= (fwd_rates[17] - rev_rates[17]);
  //sp 5
  sp_rates[5] += (fwd_rates[17] - rev_rates[17]);
  //sp 7
  sp_rates[7] += (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 8
  sp_rates[8] += (fwd_rates[18] - rev_rates[18]);
  //sp 3
  sp_rates[3] -= (fwd_rates[18] - rev_rates[18]);
  //sp 4
  sp_rates[4] += (fwd_rates[18] - rev_rates[18]);
  //sp 6
  sp_rates[6] -= (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 1
  sp_rates[1] -= (fwd_rates[19] - rev_rates[19]);
  //sp 4
  sp_rates[4] += (fwd_rates[19] - rev_rates[19]);
  //sp 6
  sp_rates[6] -= (fwd_rates[19] - rev_rates[19]);
  //sp 7
  sp_rates[7] += (fwd_rates[19] - rev_rates[19]);

  //rxn 20
  //sp 8
  sp_rates[8] += (fwd_rates[20] - rev_rates[20]);
  //sp 1
  sp_rates[1] -= (fwd_rates[20] - rev_rates[20]);
  //sp 2
  sp_rates[2] += (fwd_rates[20] - rev_rates[20]);
  //sp 6
  sp_rates[6] -= (fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 8
  sp_rates[8] += (fwd_rates[21] - rev_rates[21]);
  //sp 4
  sp_rates[4] -= (fwd_rates[21] - rev_rates[21]);
  //sp 6
  sp_rates[6] -= (fwd_rates[21] - rev_rates[21]);
  //sp 7
  sp_rates[7] += (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 8
  sp_rates[8] += (fwd_rates[22] - rev_rates[22]);
  //sp 4
  sp_rates[4] -= (fwd_rates[22] - rev_rates[22]);
  //sp 6
  sp_rates[6] -= (fwd_rates[22] - rev_rates[22]);
  //sp 7
  sp_rates[7] += (fwd_rates[22] - rev_rates[22]);

  //rxn 23
  //sp 9
  sp_rates[9] = (fwd_rates[23] - rev_rates[23]);
  //sp 18
  sp_rates[18] = -(fwd_rates[23] - rev_rates[23]);
  //sp 3
  sp_rates[3] += (fwd_rates[23] - rev_rates[23]);
  //sp 5
  sp_rates[5] -= (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 9
  sp_rates[9] += (fwd_rates[24] - rev_rates[24]);
  //sp 18
  sp_rates[18] -= (fwd_rates[24] - rev_rates[24]);
  //sp 4
  sp_rates[4] -= (fwd_rates[24] - rev_rates[24]);
  //sp 1
  sp_rates[1] += (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 1
  sp_rates[1] += (fwd_rates[25] - rev_rates[25]);
  //sp 19
  sp_rates[19] = -(fwd_rates[25] - rev_rates[25]);
  //sp 4
  sp_rates[4] -= (fwd_rates[25] - rev_rates[25]);
  //sp 20
  sp_rates[20] = (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 1
  sp_rates[1] += (fwd_rates[26] - rev_rates[26]);
  //sp 2
  sp_rates[2] -= (fwd_rates[26] - rev_rates[26]);
  //sp 19
  sp_rates[19] -= (fwd_rates[26] - rev_rates[26]);
  //sp 21
  sp_rates[21] = (fwd_rates[26] - rev_rates[26]);

  //rxn 27
  //sp 3
  sp_rates[3] -= (fwd_rates[27] - rev_rates[27]);
  //sp 9
  sp_rates[9] += (fwd_rates[27] - rev_rates[27]);
  //sp 19
  sp_rates[19] -= (fwd_rates[27] - rev_rates[27]);
  //sp 1
  sp_rates[1] += (fwd_rates[27] - rev_rates[27]);

  //rxn 28
  //sp 3
  sp_rates[3] += (fwd_rates[28] - rev_rates[28]);
  //sp 19
  sp_rates[19] -= (fwd_rates[28] - rev_rates[28]);
  //sp 20
  sp_rates[20] += (fwd_rates[28] - rev_rates[28]);
  //sp 5
  sp_rates[5] -= (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 2
  sp_rates[2] += (fwd_rates[29] - rev_rates[29]);
  //sp 1
  sp_rates[1] -= (fwd_rates[29] - rev_rates[29]);
  //sp 18
  sp_rates[18] += (fwd_rates[29] - rev_rates[29]);
  //sp 19
  sp_rates[19] -= (fwd_rates[29] - rev_rates[29]);

  //rxn 30
  //sp 10
  sp_rates[10] = (fwd_rates[30] - rev_rates[30]) * pres_mod[8];
  //sp 2
  sp_rates[2] -= (fwd_rates[30] - rev_rates[30]) * pres_mod[8];
  //sp 19
  sp_rates[19] -= (fwd_rates[30] - rev_rates[30]) * pres_mod[8];

  //rxn 31
  //sp 11
  sp_rates[11] = (fwd_rates[31] - rev_rates[31]);
  //sp 1
  sp_rates[1] += (fwd_rates[31] - rev_rates[31]);
  //sp 19
  sp_rates[19] -= (fwd_rates[31] - rev_rates[31]);
  //sp 7
  sp_rates[7] -= (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 1
  sp_rates[1] -= (fwd_rates[32] - rev_rates[32]) * pres_mod[9];
  //sp 10
  sp_rates[10] += (fwd_rates[32] - rev_rates[32]) * pres_mod[9];
  //sp 21
  sp_rates[21] -= (fwd_rates[32] - rev_rates[32]) * pres_mod[9];

  //rxn 33
  //sp 1
  sp_rates[1] += fwd_rates[33];
  //sp 4
  sp_rates[4] += fwd_rates[33];
  //sp 5
  sp_rates[5] -= fwd_rates[33];
  //sp 9
  sp_rates[9] += fwd_rates[33];
  //sp 21
  sp_rates[21] -= fwd_rates[33];

  //rxn 34
  //sp 3
  sp_rates[3] += (fwd_rates[34] - rev_rates[33]);
  //sp 21
  sp_rates[21] -= (fwd_rates[34] - rev_rates[33]);
  //sp 11
  sp_rates[11] += (fwd_rates[34] - rev_rates[33]);
  //sp 5
  sp_rates[5] -= (fwd_rates[34] - rev_rates[33]);

  //rxn 35
  //sp 1
  sp_rates[1] += (fwd_rates[35] - rev_rates[34]);
  //sp 11
  sp_rates[11] += (fwd_rates[35] - rev_rates[34]);
  //sp 4
  sp_rates[4] -= (fwd_rates[35] - rev_rates[34]);
  //sp 21
  sp_rates[21] -= (fwd_rates[35] - rev_rates[34]);

  //rxn 36
  //sp 8
  sp_rates[8] -= (fwd_rates[36] - rev_rates[35]);
  //sp 11
  sp_rates[11] += (fwd_rates[36] - rev_rates[35]);
  //sp 4
  sp_rates[4] += (fwd_rates[36] - rev_rates[35]);
  //sp 21
  sp_rates[21] -= (fwd_rates[36] - rev_rates[35]);

  //rxn 37
  //sp 1
  sp_rates[1] += 2.0 * fwd_rates[37];
  //sp 21
  sp_rates[21] -= fwd_rates[37];
  //sp 12
  sp_rates[12] = fwd_rates[37];
  //sp 5
  sp_rates[5] -= fwd_rates[37];

  //rxn 38
  //sp 19
  sp_rates[19] += (fwd_rates[38] - rev_rates[36]);
  //sp 4
  sp_rates[4] -= (fwd_rates[38] - rev_rates[36]);
  //sp 21
  sp_rates[21] -= (fwd_rates[38] - rev_rates[36]);
  //sp 7
  sp_rates[7] += (fwd_rates[38] - rev_rates[36]);

  //rxn 39
  //sp 1
  sp_rates[1] += (fwd_rates[39] - rev_rates[37]);
  //sp 3
  sp_rates[3] -= (fwd_rates[39] - rev_rates[37]);
  //sp 20
  sp_rates[20] += (fwd_rates[39] - rev_rates[37]);
  //sp 21
  sp_rates[21] -= (fwd_rates[39] - rev_rates[37]);

  //rxn 40
  //sp 1
  sp_rates[1] += (fwd_rates[40] - rev_rates[38]);
  //sp 2
  sp_rates[2] -= (fwd_rates[40] - rev_rates[38]);
  //sp 10
  sp_rates[10] += (fwd_rates[40] - rev_rates[38]);
  //sp 21
  sp_rates[21] -= (fwd_rates[40] - rev_rates[38]);

  //rxn 41
  //sp 21
  sp_rates[21] += (fwd_rates[41] - rev_rates[39]);
  //sp 22
  sp_rates[22] = -(fwd_rates[41] - rev_rates[39]);

  //rxn 42
  //sp 1
  sp_rates[1] -= (fwd_rates[42] - rev_rates[40]);
  //sp 2
  sp_rates[2] += (fwd_rates[42] - rev_rates[40]);
  //sp 19
  sp_rates[19] += (fwd_rates[42] - rev_rates[40]);
  //sp 22
  sp_rates[22] -= (fwd_rates[42] - rev_rates[40]);

  //rxn 43
  //sp 1
  sp_rates[1] += (fwd_rates[43] - rev_rates[41]);
  //sp 4
  sp_rates[4] += (fwd_rates[43] - rev_rates[41]);
  //sp 5
  sp_rates[5] -= (fwd_rates[43] - rev_rates[41]);
  //sp 9
  sp_rates[9] += (fwd_rates[43] - rev_rates[41]);
  //sp 22
  sp_rates[22] -= (fwd_rates[43] - rev_rates[41]);

  //rxn 44
  //sp 9
  sp_rates[9] += (fwd_rates[44] - rev_rates[42]);
  //sp 2
  sp_rates[2] += (fwd_rates[44] - rev_rates[42]);
  //sp 3
  sp_rates[3] -= (fwd_rates[44] - rev_rates[42]);
  //sp 22
  sp_rates[22] -= (fwd_rates[44] - rev_rates[42]);

  //rxn 45
  //sp 9
  sp_rates[9] += (fwd_rates[45] - rev_rates[43]);
  //sp 5
  sp_rates[5] -= (fwd_rates[45] - rev_rates[43]);
  //sp 22
  sp_rates[22] -= (fwd_rates[45] - rev_rates[43]);
  //sp 7
  sp_rates[7] += (fwd_rates[45] - rev_rates[43]);

  //rxn 46
  //sp 1
  sp_rates[1] += (fwd_rates[46] - rev_rates[44]);
  //sp 2
  sp_rates[2] -= (fwd_rates[46] - rev_rates[44]);
  //sp 10
  sp_rates[10] += (fwd_rates[46] - rev_rates[44]);
  //sp 22
  sp_rates[22] -= (fwd_rates[46] - rev_rates[44]);

  //rxn 47
  //sp 1
  sp_rates[1] += (fwd_rates[47] - rev_rates[45]);
  //sp 3
  sp_rates[3] -= (fwd_rates[47] - rev_rates[45]);
  //sp 20
  sp_rates[20] += (fwd_rates[47] - rev_rates[45]);
  //sp 22
  sp_rates[22] -= (fwd_rates[47] - rev_rates[45]);

  //rxn 48
  //sp 2
  sp_rates[2] += fwd_rates[48];
  //sp 11
  sp_rates[11] += fwd_rates[48];
  //sp 22
  sp_rates[22] -= fwd_rates[48];
  //sp 7
  sp_rates[7] -= fwd_rates[48];

  //rxn 49
  //sp 1
  sp_rates[1] += (fwd_rates[49] - rev_rates[46]);
  //sp 11
  sp_rates[11] += (fwd_rates[49] - rev_rates[46]);
  //sp 4
  sp_rates[4] -= (fwd_rates[49] - rev_rates[46]);
  //sp 22
  sp_rates[22] -= (fwd_rates[49] - rev_rates[46]);

  //rxn 50
  //sp 2
  sp_rates[2] += fwd_rates[50];
  //sp 10
  sp_rates[10] -= fwd_rates[50];
  //sp 11
  sp_rates[11] += fwd_rates[50];
  //sp 4
  sp_rates[4] -= fwd_rates[50];

  //rxn 51
  //sp 8
  sp_rates[8] += (fwd_rates[51] - rev_rates[47]);
  //sp 10
  sp_rates[10] -= (fwd_rates[51] - rev_rates[47]);
  //sp 13
  sp_rates[13] = (fwd_rates[51] - rev_rates[47]);
  //sp 6
  sp_rates[6] -= (fwd_rates[51] - rev_rates[47]);

  //rxn 52
  //sp 10
  sp_rates[10] -= (fwd_rates[52] - rev_rates[48]);
  //sp 11
  sp_rates[11] += (fwd_rates[52] - rev_rates[48]);
  //sp 4
  sp_rates[4] += (fwd_rates[52] - rev_rates[48]);
  //sp 5
  sp_rates[5] -= (fwd_rates[52] - rev_rates[48]);

  //rxn 53
  //sp 1
  sp_rates[1] += (fwd_rates[53] - rev_rates[49]);
  //sp 10
  sp_rates[10] -= (fwd_rates[53] - rev_rates[49]);
  //sp 19
  sp_rates[19] -= (fwd_rates[53] - rev_rates[49]);
  //sp 23
  sp_rates[23] = (fwd_rates[53] - rev_rates[49]);

  //rxn 54
  //sp 11
  sp_rates[11] += (fwd_rates[54] - rev_rates[50]);
  //sp 1
  sp_rates[1] += (fwd_rates[54] - rev_rates[50]);
  //sp 10
  sp_rates[10] -= (fwd_rates[54] - rev_rates[50]);
  //sp 3
  sp_rates[3] -= (fwd_rates[54] - rev_rates[50]);

  //rxn 55
  //sp 1
  sp_rates[1] += (fwd_rates[55] - rev_rates[51]);
  //sp 18
  sp_rates[18] -= (fwd_rates[55] - rev_rates[51]);
  //sp 10
  sp_rates[10] -= (fwd_rates[55] - rev_rates[51]);
  //sp 14
  sp_rates[14] = (fwd_rates[55] - rev_rates[51]);

  //rxn 56
  //sp 1
  sp_rates[1] -= (fwd_rates[56] - rev_rates[52]) * pres_mod[10];
  //sp 10
  sp_rates[10] -= (fwd_rates[56] - rev_rates[52]) * pres_mod[10];
  //sp 13
  sp_rates[13] += (fwd_rates[56] - rev_rates[52]) * pres_mod[10];

  //rxn 57
  //sp 10
  sp_rates[10] -= (fwd_rates[57] - rev_rates[53]);
  //sp 4
  sp_rates[4] -= (fwd_rates[57] - rev_rates[53]);
  //sp 21
  sp_rates[21] += (fwd_rates[57] - rev_rates[53]);
  //sp 7
  sp_rates[7] += (fwd_rates[57] - rev_rates[53]);

  //rxn 58
  //sp 1
  sp_rates[1] += (fwd_rates[58] - rev_rates[54]);
  //sp 10
  sp_rates[10] -= (fwd_rates[58] - rev_rates[54]);
  //sp 22
  sp_rates[22] -= (fwd_rates[58] - rev_rates[54]);
  //sp 15
  sp_rates[15] = (fwd_rates[58] - rev_rates[54]);

  //rxn 59
  //sp 10
  sp_rates[10] -= (fwd_rates[59] - rev_rates[55]);
  //sp 4
  sp_rates[4] -= (fwd_rates[59] - rev_rates[55]);
  //sp 22
  sp_rates[22] += (fwd_rates[59] - rev_rates[55]);
  //sp 7
  sp_rates[7] += (fwd_rates[59] - rev_rates[55]);

  //rxn 60
  //sp 24
  sp_rates[24] = (fwd_rates[60] - rev_rates[56]);
  //sp 1
  sp_rates[1] += (fwd_rates[60] - rev_rates[56]);
  //sp 10
  sp_rates[10] -= 2.0 * (fwd_rates[60] - rev_rates[56]);

  //rxn 61
  //sp 8
  sp_rates[8] -= (fwd_rates[61] - rev_rates[57]);
  //sp 10
  sp_rates[10] -= (fwd_rates[61] - rev_rates[57]);
  //sp 5
  sp_rates[5] += (fwd_rates[61] - rev_rates[57]);
  //sp 13
  sp_rates[13] += (fwd_rates[61] - rev_rates[57]);

  //rxn 62
  //sp 1
  sp_rates[1] += (fwd_rates[62] - rev_rates[58]);
  //sp 10
  sp_rates[10] -= (fwd_rates[62] - rev_rates[58]);
  //sp 21
  sp_rates[21] -= (fwd_rates[62] - rev_rates[58]);
  //sp 15
  sp_rates[15] += (fwd_rates[62] - rev_rates[58]);

  //rxn 63
  //sp 1
  sp_rates[1] += fwd_rates[63];
  //sp 2
  sp_rates[2] += fwd_rates[63];
  //sp 3
  sp_rates[3] -= fwd_rates[63];
  //sp 9
  sp_rates[9] += fwd_rates[63];
  //sp 10
  sp_rates[10] -= fwd_rates[63];

  //rxn 64
  //sp 1
  sp_rates[1] += (fwd_rates[64] - rev_rates[59]);
  //sp 19
  sp_rates[19] -= (fwd_rates[64] - rev_rates[59]);
  //sp 13
  sp_rates[13] -= (fwd_rates[64] - rev_rates[59]);
  //sp 15
  sp_rates[15] += (fwd_rates[64] - rev_rates[59]);

  //rxn 65
  //sp 10
  sp_rates[10] += 2.0 * (fwd_rates[65] - rev_rates[60]);
  //sp 13
  sp_rates[13] -= (fwd_rates[65] - rev_rates[60]);
  //sp 22
  sp_rates[22] -= (fwd_rates[65] - rev_rates[60]);

  //rxn 66
  //sp 10
  sp_rates[10] += (fwd_rates[66] - rev_rates[61]);
  //sp 3
  sp_rates[3] -= (fwd_rates[66] - rev_rates[61]);
  //sp 4
  sp_rates[4] += (fwd_rates[66] - rev_rates[61]);
  //sp 13
  sp_rates[13] -= (fwd_rates[66] - rev_rates[61]);

  //rxn 67
  //sp 10
  sp_rates[10] += (fwd_rates[67] - rev_rates[62]);
  //sp 4
  sp_rates[4] -= (fwd_rates[67] - rev_rates[62]);
  //sp 13
  sp_rates[13] -= (fwd_rates[67] - rev_rates[62]);
  //sp 7
  sp_rates[7] += (fwd_rates[67] - rev_rates[62]);

  //rxn 68
  //sp 21
  sp_rates[21] -= (fwd_rates[68] - rev_rates[63]);
  //sp 10
  sp_rates[10] += 2.0 * (fwd_rates[68] - rev_rates[63]);
  //sp 13
  sp_rates[13] -= (fwd_rates[68] - rev_rates[63]);

  //rxn 69
  //sp 1
  sp_rates[1] -= (fwd_rates[69] - rev_rates[64]);
  //sp 10
  sp_rates[10] += (fwd_rates[69] - rev_rates[64]);
  //sp 2
  sp_rates[2] += (fwd_rates[69] - rev_rates[64]);
  //sp 13
  sp_rates[13] -= (fwd_rates[69] - rev_rates[64]);

  //rxn 70
  //sp 16
  sp_rates[16] = (fwd_rates[70] - rev_rates[65]) * pres_mod[11];
  //sp 9
  sp_rates[9] -= (fwd_rates[70] - rev_rates[65]) * pres_mod[11];
  //sp 21
  sp_rates[21] -= (fwd_rates[70] - rev_rates[65]) * pres_mod[11];

  //rxn 71
  //sp 21
  sp_rates[21] += (fwd_rates[71] - rev_rates[66]);
  //sp 22
  sp_rates[22] -= (fwd_rates[71] - rev_rates[66]);

  //rxn 72
  //sp 9
  sp_rates[9] -= (fwd_rates[72] - rev_rates[67]);
  //sp 3
  sp_rates[3] += (fwd_rates[72] - rev_rates[67]);
  //sp 12
  sp_rates[12] += (fwd_rates[72] - rev_rates[67]);
  //sp 5
  sp_rates[5] -= (fwd_rates[72] - rev_rates[67]);

  //rxn 73
  //sp 9
  sp_rates[9] -= (fwd_rates[73] - rev_rates[68]);
  //sp 4
  sp_rates[4] -= (fwd_rates[73] - rev_rates[68]);
  //sp 12
  sp_rates[12] += (fwd_rates[73] - rev_rates[68]);
  //sp 1
  sp_rates[1] += (fwd_rates[73] - rev_rates[68]);

  //rxn 74
  //sp 9
  sp_rates[9] -= (fwd_rates[74] - rev_rates[69]) * pres_mod[12];
  //sp 2
  sp_rates[2] -= (fwd_rates[74] - rev_rates[69]) * pres_mod[12];
  //sp 11
  sp_rates[11] += (fwd_rates[74] - rev_rates[69]) * pres_mod[12];

  //rxn 75
  //sp 9
  sp_rates[9] -= (fwd_rates[75] - rev_rates[70]) * pres_mod[13];
  //sp 19
  sp_rates[19] -= (fwd_rates[75] - rev_rates[70]) * pres_mod[13];
  //sp 25
  sp_rates[25] = (fwd_rates[75] - rev_rates[70]) * pres_mod[13];

  //rxn 76
  //sp 9
  sp_rates[9] -= (fwd_rates[76] - rev_rates[71]);
  //sp 4
  sp_rates[4] -= (fwd_rates[76] - rev_rates[71]);
  //sp 12
  sp_rates[12] += (fwd_rates[76] - rev_rates[71]);
  //sp 1
  sp_rates[1] += (fwd_rates[76] - rev_rates[71]);

  //rxn 77
  //sp 9
  sp_rates[9] -= (fwd_rates[77] - rev_rates[72]) * pres_mod[14];
  //sp 3
  sp_rates[3] -= (fwd_rates[77] - rev_rates[72]) * pres_mod[14];
  //sp 12
  sp_rates[12] += (fwd_rates[77] - rev_rates[72]) * pres_mod[14];

  //rxn 78
  //sp 8
  sp_rates[8] -= (fwd_rates[78] - rev_rates[73]);
  //sp 9
  sp_rates[9] -= (fwd_rates[78] - rev_rates[73]);
  //sp 12
  sp_rates[12] += (fwd_rates[78] - rev_rates[73]);
  //sp 4
  sp_rates[4] += (fwd_rates[78] - rev_rates[73]);

  //rxn 79
  //sp 1
  sp_rates[1] -= (fwd_rates[79] - rev_rates[74]);
  //sp 2
  sp_rates[2] += (fwd_rates[79] - rev_rates[74]);
  //sp 20
  sp_rates[20] -= (fwd_rates[79] - rev_rates[74]);
  //sp 9
  sp_rates[9] += (fwd_rates[79] - rev_rates[74]);

  //rxn 80
  //sp 1
  sp_rates[1] -= (fwd_rates[80] - rev_rates[75]) * pres_mod[15];
  //sp 11
  sp_rates[11] += (fwd_rates[80] - rev_rates[75]) * pres_mod[15];
  //sp 20
  sp_rates[20] -= (fwd_rates[80] - rev_rates[75]) * pres_mod[15];

  //rxn 81
  //sp 10
  sp_rates[10] -= (fwd_rates[81] - rev_rates[76]);
  //sp 26
  sp_rates[26] = (fwd_rates[81] - rev_rates[76]);
  //sp 20
  sp_rates[20] -= (fwd_rates[81] - rev_rates[76]);

  //rxn 82
  //sp 9
  sp_rates[9] += (fwd_rates[82] - rev_rates[77]) * pres_mod[16];
  //sp 20
  sp_rates[20] -= (fwd_rates[82] - rev_rates[77]) * pres_mod[16];
  //sp 1
  sp_rates[1] += (fwd_rates[82] - rev_rates[77]) * pres_mod[16];

  //rxn 83
  //sp 9
  sp_rates[9] += (fwd_rates[83] - rev_rates[78]) * pres_mod[17];
  //sp 20
  sp_rates[20] -= (fwd_rates[83] - rev_rates[78]) * pres_mod[17];
  //sp 1
  sp_rates[1] += (fwd_rates[83] - rev_rates[78]) * pres_mod[17];

  //rxn 84
  //sp 9
  sp_rates[9] += (fwd_rates[84] - rev_rates[79]);
  //sp 3
  sp_rates[3] -= (fwd_rates[84] - rev_rates[79]);
  //sp 20
  sp_rates[20] -= (fwd_rates[84] - rev_rates[79]);
  //sp 4
  sp_rates[4] += (fwd_rates[84] - rev_rates[79]);

  //rxn 85
  //sp 9
  sp_rates[9] += (fwd_rates[85] - rev_rates[80]);
  //sp 20
  sp_rates[20] -= (fwd_rates[85] - rev_rates[80]);
  //sp 4
  sp_rates[4] -= (fwd_rates[85] - rev_rates[80]);
  //sp 7
  sp_rates[7] += (fwd_rates[85] - rev_rates[80]);

  //rxn 86
  //sp 9
  sp_rates[9] += (fwd_rates[86] - rev_rates[81]);
  //sp 10
  sp_rates[10] -= (fwd_rates[86] - rev_rates[81]);
  //sp 20
  sp_rates[20] -= (fwd_rates[86] - rev_rates[81]);
  //sp 13
  sp_rates[13] += (fwd_rates[86] - rev_rates[81]);

  //rxn 87
  //sp 1
  sp_rates[1] += (fwd_rates[87] - rev_rates[82]);
  //sp 3
  sp_rates[3] -= (fwd_rates[87] - rev_rates[82]);
  //sp 20
  sp_rates[20] -= (fwd_rates[87] - rev_rates[82]);
  //sp 12
  sp_rates[12] += (fwd_rates[87] - rev_rates[82]);

  //rxn 88
  //sp 8
  sp_rates[8] += (fwd_rates[88] - rev_rates[83]);
  //sp 9
  sp_rates[9] += (fwd_rates[88] - rev_rates[83]);
  //sp 20
  sp_rates[20] -= (fwd_rates[88] - rev_rates[83]);
  //sp 5
  sp_rates[5] -= (fwd_rates[88] - rev_rates[83]);

  //rxn 89
  //sp 1
  sp_rates[1] -= (fwd_rates[89] - rev_rates[84]);
  //sp 2
  sp_rates[2] += (fwd_rates[89] - rev_rates[84]);
  //sp 11
  sp_rates[11] -= (fwd_rates[89] - rev_rates[84]);
  //sp 20
  sp_rates[20] += (fwd_rates[89] - rev_rates[84]);

  //rxn 90
  //sp 3
  sp_rates[3] -= (fwd_rates[90] - rev_rates[85]);
  //sp 11
  sp_rates[11] -= (fwd_rates[90] - rev_rates[85]);
  //sp 20
  sp_rates[20] += (fwd_rates[90] - rev_rates[85]);
  //sp 4
  sp_rates[4] += (fwd_rates[90] - rev_rates[85]);

  //rxn 91
  //sp 10
  sp_rates[10] -= (fwd_rates[91] - rev_rates[86]);
  //sp 11
  sp_rates[11] -= (fwd_rates[91] - rev_rates[86]);
  //sp 20
  sp_rates[20] += (fwd_rates[91] - rev_rates[86]);
  //sp 13
  sp_rates[13] += (fwd_rates[91] - rev_rates[86]);

  //rxn 92
  //sp 11
  sp_rates[11] -= (fwd_rates[92] - rev_rates[87]);
  //sp 4
  sp_rates[4] -= (fwd_rates[92] - rev_rates[87]);
  //sp 20
  sp_rates[20] += (fwd_rates[92] - rev_rates[87]);
  //sp 7
  sp_rates[7] += (fwd_rates[92] - rev_rates[87]);

  //rxn 93
  //sp 11
  sp_rates[11] -= (fwd_rates[93] - rev_rates[88]);
  //sp 16
  sp_rates[16] += (fwd_rates[93] - rev_rates[88]);
  //sp 19
  sp_rates[19] -= (fwd_rates[93] - rev_rates[88]);
  //sp 1
  sp_rates[1] += (fwd_rates[93] - rev_rates[88]);

  //rxn 94
  //sp 8
  sp_rates[8] += (fwd_rates[94] - rev_rates[89]);
  //sp 11
  sp_rates[11] -= (fwd_rates[94] - rev_rates[89]);
  //sp 20
  sp_rates[20] += (fwd_rates[94] - rev_rates[89]);
  //sp 5
  sp_rates[5] -= (fwd_rates[94] - rev_rates[89]);

  //rxn 95
  //sp 8
  sp_rates[8] -= (fwd_rates[95] - rev_rates[90]);
  //sp 11
  sp_rates[11] -= (fwd_rates[95] - rev_rates[90]);
  //sp 20
  sp_rates[20] += (fwd_rates[95] - rev_rates[90]);
  //sp 6
  sp_rates[6] += (fwd_rates[95] - rev_rates[90]);

  //rxn 96
  //sp 1
  sp_rates[1] -= 2.0 * (fwd_rates[96] - rev_rates[91]) * pres_mod[18];
  //sp 2
  sp_rates[2] += (fwd_rates[96] - rev_rates[91]) * pres_mod[18];

  //rxn 97
  //sp 21
  sp_rates[21] += (fwd_rates[97] - rev_rates[92]);
  //sp 22
  sp_rates[22] -= (fwd_rates[97] - rev_rates[92]);

  //rxn 98
  //sp 9
  sp_rates[9] += (fwd_rates[98] - rev_rates[93]);
  //sp 11
  sp_rates[11] += (fwd_rates[98] - rev_rates[93]);
  //sp 12
  sp_rates[12] -= (fwd_rates[98] - rev_rates[93]);
  //sp 22
  sp_rates[22] -= (fwd_rates[98] - rev_rates[93]);

  //rxn 99
  //sp 9
  sp_rates[9] += (fwd_rates[99] - rev_rates[94]);
  //sp 19
  sp_rates[19] -= (fwd_rates[99] - rev_rates[94]);
  //sp 12
  sp_rates[12] -= (fwd_rates[99] - rev_rates[94]);
  //sp 20
  sp_rates[20] += (fwd_rates[99] - rev_rates[94]);

  //rxn 100
  //sp 9
  sp_rates[9] += (fwd_rates[100] - rev_rates[95]);
  //sp 3
  sp_rates[3] -= (fwd_rates[100] - rev_rates[95]);
  //sp 21
  sp_rates[21] += (fwd_rates[100] - rev_rates[95]);
  //sp 14
  sp_rates[14] -= (fwd_rates[100] - rev_rates[95]);

  //rxn 101
  //sp 9
  sp_rates[9] += (fwd_rates[101] - rev_rates[96]);
  //sp 10
  sp_rates[10] += (fwd_rates[101] - rev_rates[96]);
  //sp 4
  sp_rates[4] -= (fwd_rates[101] - rev_rates[96]);
  //sp 14
  sp_rates[14] -= (fwd_rates[101] - rev_rates[96]);

  //rxn 102
  //sp 1
  sp_rates[1] -= (fwd_rates[102] - rev_rates[97]) * pres_mod[19];
  //sp 14
  sp_rates[14] -= (fwd_rates[102] - rev_rates[97]) * pres_mod[19];
  //sp 23
  sp_rates[23] += (fwd_rates[102] - rev_rates[97]) * pres_mod[19];

  //rxn 103
  //sp 16
  sp_rates[16] += (fwd_rates[103] - rev_rates[98]);
  //sp 1
  sp_rates[1] += (fwd_rates[103] - rev_rates[98]);
  //sp 4
  sp_rates[4] -= (fwd_rates[103] - rev_rates[98]);
  //sp 14
  sp_rates[14] -= (fwd_rates[103] - rev_rates[98]);

  //rxn 104
  //sp 1
  sp_rates[1] += (fwd_rates[104] - rev_rates[99]);
  //sp 3
  sp_rates[3] -= (fwd_rates[104] - rev_rates[99]);
  //sp 14
  sp_rates[14] -= (fwd_rates[104] - rev_rates[99]);
  //sp 25
  sp_rates[25] += (fwd_rates[104] - rev_rates[99]);

  //rxn 105
  //sp 4
  sp_rates[4] -= (fwd_rates[105] - rev_rates[100]);
  //sp 7
  sp_rates[7] += (fwd_rates[105] - rev_rates[100]);
  //sp 14
  sp_rates[14] += (fwd_rates[105] - rev_rates[100]);
  //sp 23
  sp_rates[23] -= (fwd_rates[105] - rev_rates[100]);

  //rxn 106
  //sp 3
  sp_rates[3] += (fwd_rates[106] - rev_rates[101]);
  //sp 27
  sp_rates[27] = (fwd_rates[106] - rev_rates[101]);
  //sp 5
  sp_rates[5] -= (fwd_rates[106] - rev_rates[101]);
  //sp 23
  sp_rates[23] -= (fwd_rates[106] - rev_rates[101]);

  //rxn 107
  //sp 27
  sp_rates[27] += (fwd_rates[107] - rev_rates[102]);
  //sp 3
  sp_rates[3] -= (fwd_rates[107] - rev_rates[102]);
  //sp 23
  sp_rates[23] -= (fwd_rates[107] - rev_rates[102]);

  //rxn 108
  //sp 1
  sp_rates[1] -= (fwd_rates[108] - rev_rates[103]);
  //sp 2
  sp_rates[2] += (fwd_rates[108] - rev_rates[103]);
  //sp 14
  sp_rates[14] += (fwd_rates[108] - rev_rates[103]);
  //sp 23
  sp_rates[23] -= (fwd_rates[108] - rev_rates[103]);

  //rxn 109
  //sp 10
  sp_rates[10] -= (fwd_rates[109] - rev_rates[104]);
  //sp 13
  sp_rates[13] += (fwd_rates[109] - rev_rates[104]);
  //sp 14
  sp_rates[14] += (fwd_rates[109] - rev_rates[104]);
  //sp 23
  sp_rates[23] -= (fwd_rates[109] - rev_rates[104]);

  //rxn 110
  //sp 11
  sp_rates[11] += (fwd_rates[110] - rev_rates[105]);
  //sp 20
  sp_rates[20] += (fwd_rates[110] - rev_rates[105]);
  //sp 5
  sp_rates[5] -= (fwd_rates[110] - rev_rates[105]);
  //sp 23
  sp_rates[23] -= (fwd_rates[110] - rev_rates[105]);

  //rxn 111
  //sp 1
  sp_rates[1] -= (fwd_rates[111] - rev_rates[106]) * pres_mod[20];
  //sp 15
  sp_rates[15] += (fwd_rates[111] - rev_rates[106]) * pres_mod[20];
  //sp 23
  sp_rates[23] -= (fwd_rates[111] - rev_rates[106]) * pres_mod[20];

  //rxn 112
  //sp 8
  sp_rates[8] += (fwd_rates[112] - rev_rates[107]);
  //sp 15
  sp_rates[15] += (fwd_rates[112] - rev_rates[107]);
  //sp 6
  sp_rates[6] -= (fwd_rates[112] - rev_rates[107]);
  //sp 23
  sp_rates[23] -= (fwd_rates[112] - rev_rates[107]);

  //rxn 113
  //sp 8
  sp_rates[8] += (fwd_rates[113] - rev_rates[108]);
  //sp 5
  sp_rates[5] -= (fwd_rates[113] - rev_rates[108]);
  //sp 14
  sp_rates[14] += (fwd_rates[113] - rev_rates[108]);
  //sp 23
  sp_rates[23] -= (fwd_rates[113] - rev_rates[108]);

  //rxn 114
  //sp 10
  sp_rates[10] -= (fwd_rates[114] - rev_rates[109]);
  //sp 23
  sp_rates[23] += (fwd_rates[114] - rev_rates[109]);
  //sp 13
  sp_rates[13] += (fwd_rates[114] - rev_rates[109]);
  //sp 15
  sp_rates[15] -= (fwd_rates[114] - rev_rates[109]);

  //rxn 115
  //sp 24
  sp_rates[24] += (fwd_rates[115] - rev_rates[110]) * pres_mod[21];
  //sp 1
  sp_rates[1] -= (fwd_rates[115] - rev_rates[110]) * pres_mod[21];
  //sp 15
  sp_rates[15] -= (fwd_rates[115] - rev_rates[110]) * pres_mod[21];

  //rxn 116
  //sp 1
  sp_rates[1] += fwd_rates[116];
  //sp 5
  sp_rates[5] -= fwd_rates[116];
  //sp 10
  sp_rates[10] += fwd_rates[116];
  //sp 12
  sp_rates[12] += fwd_rates[116];
  //sp 15
  sp_rates[15] -= fwd_rates[116];

  //rxn 117
  //sp 7
  sp_rates[7] += (fwd_rates[117] - rev_rates[111]);
  //sp 4
  sp_rates[4] -= (fwd_rates[117] - rev_rates[111]);
  //sp 23
  sp_rates[23] += (fwd_rates[117] - rev_rates[111]);
  //sp 15
  sp_rates[15] -= (fwd_rates[117] - rev_rates[111]);

  //rxn 118
  //sp 4
  sp_rates[4] -= (fwd_rates[118] - rev_rates[112]);
  //sp 28
  (*dy_N) = (fwd_rates[118] - rev_rates[112]);
  //sp 15
  sp_rates[15] -= (fwd_rates[118] - rev_rates[112]);

  //rxn 119
  //sp 27
  sp_rates[27] += (fwd_rates[119] - rev_rates[113]);
  //sp 1
  sp_rates[1] += (fwd_rates[119] - rev_rates[113]);
  //sp 3
  sp_rates[3] -= (fwd_rates[119] - rev_rates[113]);
  //sp 15
  sp_rates[15] -= (fwd_rates[119] - rev_rates[113]);

  //rxn 120
  //sp 10
  sp_rates[10] += (fwd_rates[120] - rev_rates[114]);
  //sp 3
  sp_rates[3] -= (fwd_rates[120] - rev_rates[114]);
  //sp 20
  sp_rates[20] += (fwd_rates[120] - rev_rates[114]);
  //sp 15
  sp_rates[15] -= (fwd_rates[120] - rev_rates[114]);

  //rxn 121
  //sp 8
  sp_rates[8] += (fwd_rates[121] - rev_rates[115]);
  //sp 23
  sp_rates[23] += (fwd_rates[121] - rev_rates[115]);
  //sp 5
  sp_rates[5] -= (fwd_rates[121] - rev_rates[115]);
  //sp 15
  sp_rates[15] -= (fwd_rates[121] - rev_rates[115]);

  //rxn 122
  //sp 1
  sp_rates[1] -= (fwd_rates[122] - rev_rates[116]);
  //sp 2
  sp_rates[2] += (fwd_rates[122] - rev_rates[116]);
  //sp 23
  sp_rates[23] += (fwd_rates[122] - rev_rates[116]);
  //sp 15
  sp_rates[15] -= (fwd_rates[122] - rev_rates[116]);

  //rxn 123
  //sp 11
  sp_rates[11] += (fwd_rates[123] - rev_rates[117]);
  //sp 3
  sp_rates[3] -= (fwd_rates[123] - rev_rates[117]);
  //sp 21
  sp_rates[21] += (fwd_rates[123] - rev_rates[117]);
  //sp 15
  sp_rates[15] -= (fwd_rates[123] - rev_rates[117]);

  //rxn 124
  //sp 24
  sp_rates[24] -= (fwd_rates[124] - rev_rates[118]);
  //sp 8
  sp_rates[8] -= (fwd_rates[124] - rev_rates[118]);
  //sp 6
  sp_rates[6] += (fwd_rates[124] - rev_rates[118]);
  //sp 15
  sp_rates[15] += (fwd_rates[124] - rev_rates[118]);

  //rxn 125
  //sp 24
  sp_rates[24] -= (fwd_rates[125] - rev_rates[119]) * pres_mod[22];
  //sp 1
  sp_rates[1] -= (fwd_rates[125] - rev_rates[119]) * pres_mod[22];
  //sp 17
  sp_rates[17] = (fwd_rates[125] - rev_rates[119]) * pres_mod[22];

  //rxn 126
  //sp 24
  sp_rates[24] -= (fwd_rates[126] - rev_rates[120]);
  //sp 8
  sp_rates[8] -= (fwd_rates[126] - rev_rates[120]);
  //sp 28
  (*dy_N) += (fwd_rates[126] - rev_rates[120]);
  //sp 4
  sp_rates[4] += (fwd_rates[126] - rev_rates[120]);

  //rxn 127
  //sp 24
  sp_rates[24] -= (fwd_rates[127] - rev_rates[121]);
  //sp 3
  sp_rates[3] -= (fwd_rates[127] - rev_rates[121]);
  //sp 28
  (*dy_N) += (fwd_rates[127] - rev_rates[121]);

  //rxn 128
  //sp 24
  sp_rates[24] -= (fwd_rates[128] - rev_rates[122]);
  //sp 1
  sp_rates[1] -= (fwd_rates[128] - rev_rates[122]);
  //sp 2
  sp_rates[2] += (fwd_rates[128] - rev_rates[122]);
  //sp 15
  sp_rates[15] += (fwd_rates[128] - rev_rates[122]);

  //rxn 129
  //sp 24
  sp_rates[24] -= (fwd_rates[129] - rev_rates[123]);
  //sp 8
  sp_rates[8] += (fwd_rates[129] - rev_rates[123]);
  //sp 5
  sp_rates[5] -= (fwd_rates[129] - rev_rates[123]);
  //sp 15
  sp_rates[15] += (fwd_rates[129] - rev_rates[123]);

  //rxn 130
  //sp 24
  sp_rates[24] -= (fwd_rates[130] - rev_rates[124]);
  //sp 8
  sp_rates[8] -= (fwd_rates[130] - rev_rates[124]);
  //sp 5
  sp_rates[5] += (fwd_rates[130] - rev_rates[124]);
  //sp 17
  sp_rates[17] += (fwd_rates[130] - rev_rates[124]);

  //rxn 131
  //sp 24
  sp_rates[24] -= (fwd_rates[131] - rev_rates[125]);
  //sp 10
  sp_rates[10] -= (fwd_rates[131] - rev_rates[125]);
  //sp 13
  sp_rates[13] += (fwd_rates[131] - rev_rates[125]);
  //sp 15
  sp_rates[15] += (fwd_rates[131] - rev_rates[125]);

  //rxn 132
  //sp 24
  sp_rates[24] += (fwd_rates[132] - rev_rates[126]);
  //sp 17
  sp_rates[17] -= (fwd_rates[132] - rev_rates[126]);
  //sp 10
  sp_rates[10] += (fwd_rates[132] - rev_rates[126]);
  //sp 22
  sp_rates[22] -= (fwd_rates[132] - rev_rates[126]);

  //rxn 133
  //sp 24
  sp_rates[24] += (fwd_rates[133] - rev_rates[127]);
  //sp 17
  sp_rates[17] -= (fwd_rates[133] - rev_rates[127]);
  //sp 10
  sp_rates[10] -= (fwd_rates[133] - rev_rates[127]);
  //sp 13
  sp_rates[13] += (fwd_rates[133] - rev_rates[127]);

  //rxn 134
  //sp 24
  sp_rates[24] += (fwd_rates[134] - rev_rates[128]);
  //sp 17
  sp_rates[17] -= (fwd_rates[134] - rev_rates[128]);
  //sp 3
  sp_rates[3] -= (fwd_rates[134] - rev_rates[128]);
  //sp 4
  sp_rates[4] += (fwd_rates[134] - rev_rates[128]);

  //rxn 135
  //sp 17
  sp_rates[17] -= (fwd_rates[135] - rev_rates[129]) * pres_mod[23];
  //sp 10
  sp_rates[10] += 2.0 * (fwd_rates[135] - rev_rates[129]) * pres_mod[23];

  //rxn 136
  //sp 8
  sp_rates[8] -= (fwd_rates[136] - rev_rates[130]);
  //sp 17
  sp_rates[17] -= (fwd_rates[136] - rev_rates[130]);
  //sp 24
  sp_rates[24] += (fwd_rates[136] - rev_rates[130]);
  //sp 6
  sp_rates[6] += (fwd_rates[136] - rev_rates[130]);

  //rxn 137
  //sp 24
  sp_rates[24] += (fwd_rates[137] - rev_rates[131]);
  //sp 17
  sp_rates[17] -= (fwd_rates[137] - rev_rates[131]);
  //sp 2
  sp_rates[2] += (fwd_rates[137] - rev_rates[131]);
  //sp 1
  sp_rates[1] -= (fwd_rates[137] - rev_rates[131]);

  //rxn 138
  //sp 24
  sp_rates[24] += (fwd_rates[138] - rev_rates[132]);
  //sp 17
  sp_rates[17] -= (fwd_rates[138] - rev_rates[132]);
  //sp 4
  sp_rates[4] -= (fwd_rates[138] - rev_rates[132]);
  //sp 7
  sp_rates[7] += (fwd_rates[138] - rev_rates[132]);

  //rxn 139
  //sp 25
  sp_rates[25] -= (fwd_rates[139] - rev_rates[133]);
  //sp 4
  sp_rates[4] += (fwd_rates[139] - rev_rates[133]);
  //sp 5
  sp_rates[5] -= (fwd_rates[139] - rev_rates[133]);
  //sp 9
  sp_rates[9] += 2.0 * (fwd_rates[139] - rev_rates[133]);

  //rxn 140
  //sp 25
  sp_rates[25] -= (fwd_rates[140] - rev_rates[134]);
  //sp 3
  sp_rates[3] -= (fwd_rates[140] - rev_rates[134]);
  //sp 9
  sp_rates[9] += 2.0 * (fwd_rates[140] - rev_rates[134]);
  //sp 1
  sp_rates[1] += (fwd_rates[140] - rev_rates[134]);

  //rxn 141
  //sp 25
  sp_rates[25] -= (fwd_rates[141] - rev_rates[135]);
  //sp 10
  sp_rates[10] -= (fwd_rates[141] - rev_rates[135]);
  //sp 9
  sp_rates[9] += (fwd_rates[141] - rev_rates[135]);
  //sp 15
  sp_rates[15] += (fwd_rates[141] - rev_rates[135]);

  //rxn 142
  //sp 1
  sp_rates[1] -= (fwd_rates[142] - rev_rates[136]);
  //sp 22
  sp_rates[22] += (fwd_rates[142] - rev_rates[136]);
  //sp 25
  sp_rates[25] -= (fwd_rates[142] - rev_rates[136]);
  //sp 9
  sp_rates[9] += (fwd_rates[142] - rev_rates[136]);

  //rxn 143
  //sp 16
  sp_rates[16] -= (fwd_rates[143] - rev_rates[137]);
  //sp 1
  sp_rates[1] -= (fwd_rates[143] - rev_rates[137]);
  //sp 10
  sp_rates[10] += (fwd_rates[143] - rev_rates[137]);
  //sp 9
  sp_rates[9] += (fwd_rates[143] - rev_rates[137]);

  //rxn 144
  //sp 16
  sp_rates[16] -= (fwd_rates[144] - rev_rates[138]);
  //sp 9
  sp_rates[9] += (fwd_rates[144] - rev_rates[138]);
  //sp 21
  sp_rates[21] -= (fwd_rates[144] - rev_rates[138]);
  //sp 15
  sp_rates[15] += (fwd_rates[144] - rev_rates[138]);

  //rxn 145
  //sp 16
  sp_rates[16] -= (fwd_rates[145] - rev_rates[139]);
  //sp 25
  sp_rates[25] += (fwd_rates[145] - rev_rates[139]);
  //sp 3
  sp_rates[3] -= (fwd_rates[145] - rev_rates[139]);
  //sp 4
  sp_rates[4] += (fwd_rates[145] - rev_rates[139]);

  //rxn 146
  //sp 16
  sp_rates[16] -= (fwd_rates[146] - rev_rates[140]);
  //sp 25
  sp_rates[25] += (fwd_rates[146] - rev_rates[140]);
  //sp 10
  sp_rates[10] -= (fwd_rates[146] - rev_rates[140]);
  //sp 13
  sp_rates[13] += (fwd_rates[146] - rev_rates[140]);

  //rxn 147
  //sp 16
  sp_rates[16] -= (fwd_rates[147] - rev_rates[141]);
  //sp 3
  sp_rates[3] -= (fwd_rates[147] - rev_rates[141]);
  //sp 12
  sp_rates[12] += (fwd_rates[147] - rev_rates[141]);
  //sp 21
  sp_rates[21] += (fwd_rates[147] - rev_rates[141]);

  //rxn 148
  //sp 16
  sp_rates[16] -= (fwd_rates[148] - rev_rates[142]);
  //sp 24
  sp_rates[24] += (fwd_rates[148] - rev_rates[142]);
  //sp 10
  sp_rates[10] -= (fwd_rates[148] - rev_rates[142]);
  //sp 9
  sp_rates[9] += (fwd_rates[148] - rev_rates[142]);

  //rxn 149
  //sp 16
  sp_rates[16] -= (fwd_rates[149] - rev_rates[143]);
  //sp 25
  sp_rates[25] += (fwd_rates[149] - rev_rates[143]);
  //sp 4
  sp_rates[4] -= (fwd_rates[149] - rev_rates[143]);
  //sp 7
  sp_rates[7] += (fwd_rates[149] - rev_rates[143]);

  //rxn 150
  //sp 16
  sp_rates[16] -= (fwd_rates[150] - rev_rates[144]);
  //sp 1
  sp_rates[1] -= (fwd_rates[150] - rev_rates[144]);
  //sp 2
  sp_rates[2] += (fwd_rates[150] - rev_rates[144]);
  //sp 25
  sp_rates[25] += (fwd_rates[150] - rev_rates[144]);

  //rxn 151
  //sp 16
  sp_rates[16] -= (fwd_rates[151] - rev_rates[145]);
  //sp 25
  sp_rates[25] += (fwd_rates[151] - rev_rates[145]);
  //sp 10
  sp_rates[10] += (fwd_rates[151] - rev_rates[145]);
  //sp 21
  sp_rates[21] -= (fwd_rates[151] - rev_rates[145]);

  //rxn 152
  //sp 3
  sp_rates[3] -= (fwd_rates[152] - rev_rates[146]);
  //sp 11
  sp_rates[11] += (fwd_rates[152] - rev_rates[146]);
  //sp 27
  sp_rates[27] -= (fwd_rates[152] - rev_rates[146]);
  //sp 20
  sp_rates[20] += (fwd_rates[152] - rev_rates[146]);

  //rxn 153
  //sp 16
  sp_rates[16] += (fwd_rates[153] - rev_rates[147]);
  //sp 1
  sp_rates[1] += (fwd_rates[153] - rev_rates[147]);
  //sp 27
  sp_rates[27] -= (fwd_rates[153] - rev_rates[147]);

  //rxn 154
  //sp 16
  sp_rates[16] += (fwd_rates[154] - rev_rates[148]);
  //sp 27
  sp_rates[27] -= (fwd_rates[154] - rev_rates[148]);
  //sp 4
  sp_rates[4] -= (fwd_rates[154] - rev_rates[148]);
  //sp 7
  sp_rates[7] += (fwd_rates[154] - rev_rates[148]);

  //rxn 155
  //sp 16
  sp_rates[16] += (fwd_rates[155] - rev_rates[149]);
  //sp 1
  sp_rates[1] -= (fwd_rates[155] - rev_rates[149]);
  //sp 2
  sp_rates[2] += (fwd_rates[155] - rev_rates[149]);
  //sp 27
  sp_rates[27] -= (fwd_rates[155] - rev_rates[149]);

  //rxn 156
  //sp 4
  sp_rates[4] += fwd_rates[156];
  //sp 5
  sp_rates[5] -= fwd_rates[156];
  //sp 9
  sp_rates[9] += fwd_rates[156];
  //sp 11
  sp_rates[11] += fwd_rates[156];
  //sp 27
  sp_rates[27] -= fwd_rates[156];

  //rxn 157
  //sp 9
  sp_rates[9] += (fwd_rates[157] - rev_rates[150]);
  //sp 10
  sp_rates[10] += (fwd_rates[157] - rev_rates[150]);
  //sp 27
  sp_rates[27] -= (fwd_rates[157] - rev_rates[150]);

  //rxn 158
  //sp 4
  sp_rates[4] += fwd_rates[158];
  //sp 27
  sp_rates[27] -= fwd_rates[158];
  //sp 20
  sp_rates[20] += 2.0 * fwd_rates[158];
  //sp 5
  sp_rates[5] -= fwd_rates[158];

  //rxn 159
  //sp 1
  sp_rates[1] -= (fwd_rates[159] - rev_rates[151]);
  //sp 10
  sp_rates[10] += (fwd_rates[159] - rev_rates[151]);
  //sp 27
  sp_rates[27] -= (fwd_rates[159] - rev_rates[151]);
  //sp 20
  sp_rates[20] += (fwd_rates[159] - rev_rates[151]);

  //rxn 160
  //sp 3
  sp_rates[3] -= fwd_rates[160];
  //sp 4
  sp_rates[4] += fwd_rates[160];
  //sp 9
  sp_rates[9] += fwd_rates[160];
  //sp 10
  sp_rates[10] += fwd_rates[160];
  //sp 26
  sp_rates[26] -= fwd_rates[160];

  //rxn 161
  //sp 5
  sp_rates[5] -= fwd_rates[161];
  //sp 8
  sp_rates[8] += fwd_rates[161];
  //sp 9
  sp_rates[9] += fwd_rates[161];
  //sp 10
  sp_rates[10] += fwd_rates[161];
  //sp 26
  sp_rates[26] -= fwd_rates[161];

  //rxn 162
  //sp 4
  sp_rates[4] -= fwd_rates[162];
  //sp 7
  sp_rates[7] += fwd_rates[162];
  //sp 9
  sp_rates[9] += fwd_rates[162];
  //sp 10
  sp_rates[10] += fwd_rates[162];
  //sp 26
  sp_rates[26] -= fwd_rates[162];

  //rxn 163
  //sp 2
  sp_rates[2] += (fwd_rates[163] - rev_rates[152]);
  //sp 1
  sp_rates[1] -= (fwd_rates[163] - rev_rates[152]);
  //sp 26
  sp_rates[26] -= (fwd_rates[163] - rev_rates[152]);
  //sp 27
  sp_rates[27] += (fwd_rates[163] - rev_rates[152]);

  //rxn 164
  //sp 1
  sp_rates[1] -= fwd_rates[164];
  //sp 2
  sp_rates[2] += fwd_rates[164];
  //sp 9
  sp_rates[9] += fwd_rates[164];
  //sp 10
  sp_rates[10] += fwd_rates[164];
  //sp 26
  sp_rates[26] -= fwd_rates[164];

  //rxn 165
  //sp 27
  sp_rates[27] += (fwd_rates[165] - rev_rates[153]);
  //sp 26
  sp_rates[26] -= (fwd_rates[165] - rev_rates[153]);
  //sp 3
  sp_rates[3] -= (fwd_rates[165] - rev_rates[153]);
  //sp 4
  sp_rates[4] += (fwd_rates[165] - rev_rates[153]);

  //rxn 166
  //sp 9
  sp_rates[9] += fwd_rates[166] * pres_mod[24];
  //sp 26
  sp_rates[26] -= fwd_rates[166] * pres_mod[24];
  //sp 13
  sp_rates[13] += fwd_rates[166] * pres_mod[24];

  //rxn 167
  //sp 6
  sp_rates[6] += fwd_rates[167];
  //sp 8
  sp_rates[8] -= fwd_rates[167];
  //sp 9
  sp_rates[9] += fwd_rates[167];
  //sp 10
  sp_rates[10] += fwd_rates[167];
  //sp 26
  sp_rates[26] -= fwd_rates[167];

  //rxn 168
  //sp 10
  sp_rates[10] += (fwd_rates[168] - rev_rates[154]);
  //sp 11
  sp_rates[11] += (fwd_rates[168] - rev_rates[154]);
  //sp 28
  (*dy_N) -= (fwd_rates[168] - rev_rates[154]);

  //rxn 169
  //sp 1
  sp_rates[1] += (fwd_rates[169] - rev_rates[155]);
  //sp 26
  sp_rates[26] += (fwd_rates[169] - rev_rates[155]);
  //sp 28
  (*dy_N) -= (fwd_rates[169] - rev_rates[155]);

  //rxn 170
  //sp 8
  sp_rates[8] += (fwd_rates[170] - rev_rates[156]);
  //sp 26
  sp_rates[26] += (fwd_rates[170] - rev_rates[156]);
  //sp 28
  (*dy_N) -= (fwd_rates[170] - rev_rates[156]);
  //sp 5
  sp_rates[5] -= (fwd_rates[170] - rev_rates[156]);

  //rxn 171
  //sp 21
  sp_rates[21] += (fwd_rates[171] - rev_rates[157]);
  //sp 22
  sp_rates[22] -= (fwd_rates[171] - rev_rates[157]);

  //sp 0
  sp_rates[0] = 0.0;
} // end eval_spec_rates

