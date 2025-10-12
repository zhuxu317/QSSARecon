#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 3
  sp_rates[2] = -(fwd_rates[0] - rev_rates[0]);
  //sp 4
  sp_rates[3] = (fwd_rates[0] - rev_rates[0]);
  //sp 5
  sp_rates[4] = (fwd_rates[0] - rev_rates[0]);
  //sp 7
  sp_rates[6] = -(fwd_rates[0] - rev_rates[0]);

  //rxn 1
  //sp 3
  sp_rates[2] += (fwd_rates[1] - rev_rates[1]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1] - rev_rates[1]);
  //sp 5
  sp_rates[4] += (fwd_rates[1] - rev_rates[1]);
  //sp 6
  sp_rates[5] = -(fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 3
  sp_rates[2] += (fwd_rates[2] - rev_rates[2]);
  //sp 4
  sp_rates[3] -= (fwd_rates[2] - rev_rates[2]);
  //sp 5
  sp_rates[4] += (fwd_rates[2] - rev_rates[2]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 10
  sp_rates[9] = (fwd_rates[3] - rev_rates[3]);
  //sp 3
  sp_rates[2] += (fwd_rates[3] - rev_rates[3]);
  //sp 5
  sp_rates[4] -= (fwd_rates[3] - rev_rates[3]);
  //sp 6
  sp_rates[5] -= (fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 10
  sp_rates[9] += (fwd_rates[4] - rev_rates[4]);
  //sp 4
  sp_rates[3] += (fwd_rates[4] - rev_rates[4]);
  //sp 5
  sp_rates[4] -= 2.0 * (fwd_rates[4] - rev_rates[4]);

  //rxn 5
  //sp 3
  sp_rates[2] -= 2.0 * (fwd_rates[5] - rev_rates[5]) * pres_mod[0];
  //sp 6
  sp_rates[5] += (fwd_rates[5] - rev_rates[5]) * pres_mod[0];

  //rxn 6
  //sp 4
  sp_rates[3] -= 2.0 * (fwd_rates[6] - rev_rates[6]) * pres_mod[1];
  //sp 7
  sp_rates[6] += (fwd_rates[6] - rev_rates[6]) * pres_mod[1];

  //rxn 7
  //sp 3
  sp_rates[2] -= (fwd_rates[7] - rev_rates[7]) * pres_mod[2];
  //sp 4
  sp_rates[3] -= (fwd_rates[7] - rev_rates[7]) * pres_mod[2];
  //sp 5
  sp_rates[4] += (fwd_rates[7] - rev_rates[7]) * pres_mod[2];

  //rxn 8
  //sp 10
  sp_rates[9] += (fwd_rates[8] - rev_rates[8]) * pres_mod[3];
  //sp 3
  sp_rates[2] -= (fwd_rates[8] - rev_rates[8]) * pres_mod[3];
  //sp 5
  sp_rates[4] -= (fwd_rates[8] - rev_rates[8]) * pres_mod[3];

  //rxn 9
  //sp 3
  sp_rates[2] -= (fwd_rates[9] - rev_rates[9]) * pres_mod[4];
  //sp 7
  sp_rates[6] -= (fwd_rates[9] - rev_rates[9]) * pres_mod[4];
  //sp 8
  sp_rates[7] = (fwd_rates[9] - rev_rates[9]) * pres_mod[4];

  //rxn 10
  //sp 3
  sp_rates[2] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[5];
  //sp 7
  sp_rates[6] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[5];
  //sp 8
  sp_rates[7] += (fwd_rates[10] - rev_rates[10]) * pres_mod[5];

  //rxn 11
  //sp 3
  sp_rates[2] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[6];
  //sp 7
  sp_rates[6] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[6];
  //sp 8
  sp_rates[7] += (fwd_rates[11] - rev_rates[11]) * pres_mod[6];

  //rxn 12
  //sp 3
  sp_rates[2] -= (fwd_rates[12] - rev_rates[12]) * pres_mod[7];
  //sp 7
  sp_rates[6] -= (fwd_rates[12] - rev_rates[12]) * pres_mod[7];
  //sp 8
  sp_rates[7] += (fwd_rates[12] - rev_rates[12]) * pres_mod[7];

  //rxn 13
  //sp 3
  sp_rates[2] -= (fwd_rates[13] - rev_rates[13]);
  //sp 6
  sp_rates[5] += (fwd_rates[13] - rev_rates[13]);
  //sp 7
  sp_rates[6] += (fwd_rates[13] - rev_rates[13]);
  //sp 8
  sp_rates[7] -= (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 3
  sp_rates[2] -= (fwd_rates[14] - rev_rates[14]);
  //sp 5
  sp_rates[4] += 2.0 * (fwd_rates[14] - rev_rates[14]);
  //sp 8
  sp_rates[7] -= (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 4
  sp_rates[3] -= (fwd_rates[15] - rev_rates[15]);
  //sp 5
  sp_rates[4] += (fwd_rates[15] - rev_rates[15]);
  //sp 7
  sp_rates[6] += (fwd_rates[15] - rev_rates[15]);
  //sp 8
  sp_rates[7] -= (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 4
  sp_rates[3] -= (fwd_rates[16] - rev_rates[16]) * pres_mod[8];
  //sp 5
  sp_rates[4] -= (fwd_rates[16] - rev_rates[16]) * pres_mod[8];
  //sp 8
  sp_rates[7] += (fwd_rates[16] - rev_rates[16]) * pres_mod[8];

  //rxn 17
  //sp 5
  sp_rates[4] += 2.0 * (fwd_rates[17] - rev_rates[17]);
  //sp 6
  sp_rates[5] -= (fwd_rates[17] - rev_rates[17]);
  //sp 7
  sp_rates[6] -= (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 10
  sp_rates[9] += (fwd_rates[18] - rev_rates[18]);
  //sp 5
  sp_rates[4] -= (fwd_rates[18] - rev_rates[18]);
  //sp 7
  sp_rates[6] += (fwd_rates[18] - rev_rates[18]);
  //sp 8
  sp_rates[7] -= (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 10
  sp_rates[9] += (fwd_rates[19] - rev_rates[19]);
  //sp 5
  sp_rates[4] -= (fwd_rates[19] - rev_rates[19]);
  //sp 7
  sp_rates[6] += (fwd_rates[19] - rev_rates[19]);
  //sp 8
  sp_rates[7] -= (fwd_rates[19] - rev_rates[19]);

  //rxn 20
  //sp 9
  sp_rates[8] = (fwd_rates[20] - rev_rates[20]);
  //sp 7
  sp_rates[6] += (fwd_rates[20] - rev_rates[20]);
  //sp 8
  sp_rates[7] -= 2.0 * (fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 9
  sp_rates[8] += (fwd_rates[21] - rev_rates[21]);
  //sp 7
  sp_rates[6] += (fwd_rates[21] - rev_rates[21]);
  //sp 8
  sp_rates[7] -= 2.0 * (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 9
  sp_rates[8] -= (fwd_rates[22] - rev_rates[22]) * pres_mod[9];
  //sp 5
  sp_rates[4] += 2.0 * (fwd_rates[22] - rev_rates[22]) * pres_mod[9];

  //rxn 23
  //sp 9
  sp_rates[8] -= (fwd_rates[23] - rev_rates[23]);
  //sp 10
  sp_rates[9] += (fwd_rates[23] - rev_rates[23]);
  //sp 3
  sp_rates[2] -= (fwd_rates[23] - rev_rates[23]);
  //sp 5
  sp_rates[4] += (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 9
  sp_rates[8] -= (fwd_rates[24] - rev_rates[24]);
  //sp 3
  sp_rates[2] -= (fwd_rates[24] - rev_rates[24]);
  //sp 6
  sp_rates[5] += (fwd_rates[24] - rev_rates[24]);
  //sp 8
  sp_rates[7] += (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 9
  sp_rates[8] -= (fwd_rates[25] - rev_rates[25]);
  //sp 4
  sp_rates[3] -= (fwd_rates[25] - rev_rates[25]);
  //sp 5
  sp_rates[4] += (fwd_rates[25] - rev_rates[25]);
  //sp 8
  sp_rates[7] += (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 9
  sp_rates[8] -= (fwd_rates[26] - rev_rates[26]);
  //sp 10
  sp_rates[9] += (fwd_rates[26] - rev_rates[26]);
  //sp 5
  sp_rates[4] -= (fwd_rates[26] - rev_rates[26]);
  //sp 8
  sp_rates[7] += (fwd_rates[26] - rev_rates[26]);

  //rxn 27
  //sp 9
  sp_rates[8] -= (fwd_rates[27] - rev_rates[27]);
  //sp 10
  sp_rates[9] += (fwd_rates[27] - rev_rates[27]);
  //sp 5
  sp_rates[4] -= (fwd_rates[27] - rev_rates[27]);
  //sp 8
  sp_rates[7] += (fwd_rates[27] - rev_rates[27]);

  //rxn 28
  //sp 9
  sp_rates[8] -= (fwd_rates[28] - rev_rates[28]);
  //sp 10
  sp_rates[9] += (fwd_rates[28] - rev_rates[28]);
  //sp 4
  sp_rates[3] -= (fwd_rates[28] - rev_rates[28]);
  //sp 7
  sp_rates[6] += (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 10
  sp_rates[9] -= (fwd_rates[29] - rev_rates[29]);
  //sp 4
  sp_rates[3] -= (fwd_rates[29] - rev_rates[29]);
  //sp 6
  sp_rates[5] += (fwd_rates[29] - rev_rates[29]);
  //sp 7
  sp_rates[6] += (fwd_rates[29] - rev_rates[29]);

  //rxn 30
  //sp 4
  sp_rates[3] -= (fwd_rates[30] - rev_rates[30]) * pres_mod[10];
  //sp 12
  sp_rates[11] = -(fwd_rates[30] - rev_rates[30]) * pres_mod[10];
  //sp 13
  sp_rates[12] = (fwd_rates[30] - rev_rates[30]) * pres_mod[10];

  //rxn 31
  //sp 4
  sp_rates[3] += (fwd_rates[31] - rev_rates[31]);
  //sp 12
  sp_rates[11] -= (fwd_rates[31] - rev_rates[31]);
  //sp 13
  sp_rates[12] += (fwd_rates[31] - rev_rates[31]);
  //sp 7
  sp_rates[6] -= (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 3
  sp_rates[2] += (fwd_rates[32] - rev_rates[32]);
  //sp 12
  sp_rates[11] -= (fwd_rates[32] - rev_rates[32]);
  //sp 5
  sp_rates[4] -= (fwd_rates[32] - rev_rates[32]);
  //sp 13
  sp_rates[12] += (fwd_rates[32] - rev_rates[32]);

  //rxn 33
  //sp 3
  sp_rates[2] += (fwd_rates[33] - rev_rates[33]);
  //sp 12
  sp_rates[11] -= (fwd_rates[33] - rev_rates[33]);
  //sp 5
  sp_rates[4] -= (fwd_rates[33] - rev_rates[33]);
  //sp 13
  sp_rates[12] += (fwd_rates[33] - rev_rates[33]);

  //rxn 34
  //sp 3
  sp_rates[2] += (fwd_rates[34] - rev_rates[34]);
  //sp 12
  sp_rates[11] -= (fwd_rates[34] - rev_rates[34]);
  //sp 5
  sp_rates[4] -= (fwd_rates[34] - rev_rates[34]);
  //sp 13
  sp_rates[12] += (fwd_rates[34] - rev_rates[34]);

  //rxn 35
  //sp 12
  sp_rates[11] -= (fwd_rates[35] - rev_rates[35]);
  //sp 13
  sp_rates[12] += (fwd_rates[35] - rev_rates[35]);
  //sp 5
  sp_rates[4] += (fwd_rates[35] - rev_rates[35]);
  //sp 8
  sp_rates[7] -= (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 3
  sp_rates[2] += (fwd_rates[36] - rev_rates[36]) * pres_mod[11];
  //sp 12
  sp_rates[11] += (fwd_rates[36] - rev_rates[36]) * pres_mod[11];
  //sp 14
  sp_rates[13] = -(fwd_rates[36] - rev_rates[36]) * pres_mod[11];

  //rxn 37
  //sp 12
  sp_rates[11] += (fwd_rates[37] - rev_rates[37]);
  //sp 14
  sp_rates[13] -= (fwd_rates[37] - rev_rates[37]);
  //sp 7
  sp_rates[6] -= (fwd_rates[37] - rev_rates[37]);
  //sp 8
  sp_rates[7] += (fwd_rates[37] - rev_rates[37]);

  //rxn 38
  //sp 6
  sp_rates[5] += (fwd_rates[38] - rev_rates[38]);
  //sp 3
  sp_rates[2] -= (fwd_rates[38] - rev_rates[38]);
  //sp 12
  sp_rates[11] += (fwd_rates[38] - rev_rates[38]);
  //sp 14
  sp_rates[13] -= (fwd_rates[38] - rev_rates[38]);

  //rxn 39
  //sp 12
  sp_rates[11] += (fwd_rates[39] - rev_rates[39]);
  //sp 4
  sp_rates[3] -= (fwd_rates[39] - rev_rates[39]);
  //sp 5
  sp_rates[4] += (fwd_rates[39] - rev_rates[39]);
  //sp 14
  sp_rates[13] -= (fwd_rates[39] - rev_rates[39]);

  //rxn 40
  //sp 3
  sp_rates[2] += (fwd_rates[40] - rev_rates[40]);
  //sp 4
  sp_rates[3] -= (fwd_rates[40] - rev_rates[40]);
  //sp 13
  sp_rates[12] += (fwd_rates[40] - rev_rates[40]);
  //sp 14
  sp_rates[13] -= (fwd_rates[40] - rev_rates[40]);

  //rxn 41
  //sp 10
  sp_rates[9] += (fwd_rates[41] - rev_rates[41]);
  //sp 12
  sp_rates[11] += (fwd_rates[41] - rev_rates[41]);
  //sp 5
  sp_rates[4] -= (fwd_rates[41] - rev_rates[41]);
  //sp 14
  sp_rates[13] -= (fwd_rates[41] - rev_rates[41]);

  //rxn 42
  //sp 3
  sp_rates[2] += fwd_rates[42];
  //sp 5
  sp_rates[4] += fwd_rates[42];
  //sp 8
  sp_rates[7] -= fwd_rates[42];
  //sp 13
  sp_rates[12] += fwd_rates[42];
  //sp 14
  sp_rates[13] -= fwd_rates[42];

  //rxn 43
  //sp 6
  sp_rates[5] += fwd_rates[43];
  //sp 12
  sp_rates[11] += 2.0 * fwd_rates[43];
  //sp 14
  sp_rates[13] -= 2.0 * fwd_rates[43];

  //rxn 44
  //sp 12
  sp_rates[11] += (fwd_rates[44] - rev_rates[42]);
  //sp 14
  sp_rates[13] -= 2.0 * (fwd_rates[44] - rev_rates[42]);
  //sp 15
  sp_rates[14] = (fwd_rates[44] - rev_rates[42]);

  //rxn 45
  //sp 9
  sp_rates[8] += (fwd_rates[45] - rev_rates[43]);
  //sp 12
  sp_rates[11] += (fwd_rates[45] - rev_rates[43]);
  //sp 14
  sp_rates[13] -= (fwd_rates[45] - rev_rates[43]);
  //sp 8
  sp_rates[7] -= (fwd_rates[45] - rev_rates[43]);

  //rxn 46
  //sp 17
  sp_rates[16] = (fwd_rates[46] - rev_rates[44]);
  //sp 3
  sp_rates[2] -= (fwd_rates[46] - rev_rates[44]);
  //sp 6
  sp_rates[5] += (fwd_rates[46] - rev_rates[44]);
  //sp 16
  sp_rates[15] = -(fwd_rates[46] - rev_rates[44]);

  //rxn 47
  //sp 17
  sp_rates[16] -= (fwd_rates[47] - rev_rates[45]);
  //sp 4
  sp_rates[3] += (fwd_rates[47] - rev_rates[45]);
  //sp 12
  sp_rates[11] += (fwd_rates[47] - rev_rates[45]);
  //sp 7
  sp_rates[6] -= (fwd_rates[47] - rev_rates[45]);

  //rxn 48
  //sp 17
  sp_rates[16] -= (fwd_rates[48] - rev_rates[46]);
  //sp 3
  sp_rates[2] += (fwd_rates[48] - rev_rates[46]);
  //sp 12
  sp_rates[11] += (fwd_rates[48] - rev_rates[46]);
  //sp 5
  sp_rates[4] -= (fwd_rates[48] - rev_rates[46]);

  //rxn 49
  //sp 12
  sp_rates[11] += (fwd_rates[49] - rev_rates[47]);
  //sp 3
  sp_rates[2] += (fwd_rates[49] - rev_rates[47]);
  //sp 4
  sp_rates[3] -= (fwd_rates[49] - rev_rates[47]);
  //sp 16
  sp_rates[15] -= (fwd_rates[49] - rev_rates[47]);

  //rxn 50
  //sp 12
  sp_rates[11] += (fwd_rates[50] - rev_rates[48]);
  //sp 13
  sp_rates[12] -= (fwd_rates[50] - rev_rates[48]);
  //sp 14
  sp_rates[13] += (fwd_rates[50] - rev_rates[48]);
  //sp 16
  sp_rates[15] -= (fwd_rates[50] - rev_rates[48]);

  //rxn 51
  //sp 3
  sp_rates[2] += (fwd_rates[51] - rev_rates[49]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[51] - rev_rates[49]);
  //sp 13
  sp_rates[12] -= (fwd_rates[51] - rev_rates[49]);
  //sp 16
  sp_rates[15] -= (fwd_rates[51] - rev_rates[49]);

  //rxn 52
  //sp 10
  sp_rates[9] -= (fwd_rates[52] - rev_rates[50]);
  //sp 5
  sp_rates[4] += (fwd_rates[52] - rev_rates[50]);
  //sp 18
  sp_rates[17] = (fwd_rates[52] - rev_rates[50]);
  //sp 16
  sp_rates[15] -= (fwd_rates[52] - rev_rates[50]);

  //rxn 53
  //sp 10
  sp_rates[9] -= (fwd_rates[53] - rev_rates[51]);
  //sp 3
  sp_rates[2] += (fwd_rates[53] - rev_rates[51]);
  //sp 15
  sp_rates[14] += (fwd_rates[53] - rev_rates[51]);
  //sp 16
  sp_rates[15] -= (fwd_rates[53] - rev_rates[51]);

  //rxn 54
  //sp 4
  sp_rates[3] += (fwd_rates[54] - rev_rates[52]);
  //sp 14
  sp_rates[13] += (fwd_rates[54] - rev_rates[52]);
  //sp 7
  sp_rates[6] -= (fwd_rates[54] - rev_rates[52]);
  //sp 16
  sp_rates[15] -= (fwd_rates[54] - rev_rates[52]);

  //rxn 55
  //sp 12
  sp_rates[11] += (fwd_rates[55] - rev_rates[53]);
  //sp 5
  sp_rates[4] += (fwd_rates[55] - rev_rates[53]);
  //sp 7
  sp_rates[6] -= (fwd_rates[55] - rev_rates[53]);
  //sp 16
  sp_rates[15] -= (fwd_rates[55] - rev_rates[53]);

  //rxn 56
  //sp 3
  sp_rates[2] += (fwd_rates[56] - rev_rates[54]);
  //sp 13
  sp_rates[12] += (fwd_rates[56] - rev_rates[54]);
  //sp 7
  sp_rates[6] -= (fwd_rates[56] - rev_rates[54]);
  //sp 16
  sp_rates[15] -= (fwd_rates[56] - rev_rates[54]);

  //rxn 57
  //sp 3
  sp_rates[2] += (fwd_rates[57] - rev_rates[55]);
  //sp 4
  sp_rates[3] += (fwd_rates[57] - rev_rates[55]);
  //sp 7
  sp_rates[6] -= (fwd_rates[57] - rev_rates[55]);
  //sp 12
  sp_rates[11] += (fwd_rates[57] - rev_rates[55]);
  //sp 16
  sp_rates[15] -= (fwd_rates[57] - rev_rates[55]);

  //rxn 58
  //sp 3
  sp_rates[2] += (fwd_rates[58] - rev_rates[56]);
  //sp 5
  sp_rates[4] -= (fwd_rates[58] - rev_rates[56]);
  //sp 14
  sp_rates[13] += (fwd_rates[58] - rev_rates[56]);
  //sp 16
  sp_rates[15] -= (fwd_rates[58] - rev_rates[56]);

  //rxn 59
  //sp 18
  sp_rates[17] -= (fwd_rates[59] - rev_rates[57]);
  //sp 3
  sp_rates[2] -= (fwd_rates[59] - rev_rates[57]);
  //sp 6
  sp_rates[5] += (fwd_rates[59] - rev_rates[57]);
  //sp 16
  sp_rates[15] += (fwd_rates[59] - rev_rates[57]);

  //rxn 60
  //sp 3
  sp_rates[2] += fwd_rates[60];
  //sp 5
  sp_rates[4] += fwd_rates[60];
  //sp 7
  sp_rates[6] -= fwd_rates[60];
  //sp 12
  sp_rates[11] += fwd_rates[60];
  //sp 18
  sp_rates[17] -= fwd_rates[60];

  //rxn 61
  //sp 18
  sp_rates[17] -= (fwd_rates[61] - rev_rates[58]);
  //sp 13
  sp_rates[12] += (fwd_rates[61] - rev_rates[58]);
  //sp 6
  sp_rates[5] += (fwd_rates[61] - rev_rates[58]);
  //sp 7
  sp_rates[6] -= (fwd_rates[61] - rev_rates[58]);

  //rxn 62
  //sp 18
  sp_rates[17] -= (fwd_rates[62] - rev_rates[59]);
  //sp 12
  sp_rates[11] += (fwd_rates[62] - rev_rates[59]);
  //sp 7
  sp_rates[6] -= (fwd_rates[62] - rev_rates[59]);
  //sp 10
  sp_rates[9] += (fwd_rates[62] - rev_rates[59]);

  //rxn 63
  //sp 18
  sp_rates[17] -= (fwd_rates[63] - rev_rates[60]);
  //sp 4
  sp_rates[3] += (fwd_rates[63] - rev_rates[60]);
  //sp 7
  sp_rates[6] -= (fwd_rates[63] - rev_rates[60]);
  //sp 15
  sp_rates[14] += (fwd_rates[63] - rev_rates[60]);

  //rxn 64
  //sp 18
  sp_rates[17] -= fwd_rates[64];
  //sp 3
  sp_rates[2] += 2.0 * fwd_rates[64];
  //sp 13
  sp_rates[12] += fwd_rates[64];
  //sp 7
  sp_rates[6] -= fwd_rates[64];

  //rxn 65
  //sp 18
  sp_rates[17] -= 2.0 * (fwd_rates[65] - rev_rates[61]);
  //sp 20
  sp_rates[19] = (fwd_rates[65] - rev_rates[61]);
  //sp 6
  sp_rates[5] += (fwd_rates[65] - rev_rates[61]);

  //rxn 66
  //sp 18
  sp_rates[17] -= 2.0 * (fwd_rates[66] - rev_rates[62]);
  //sp 3
  sp_rates[2] += 2.0 * (fwd_rates[66] - rev_rates[62]);
  //sp 20
  sp_rates[19] += (fwd_rates[66] - rev_rates[62]);

  //rxn 67
  //sp 18
  sp_rates[17] -= (fwd_rates[67] - rev_rates[63]);
  //sp 5
  sp_rates[4] += (fwd_rates[67] - rev_rates[63]);
  //sp 14
  sp_rates[13] += (fwd_rates[67] - rev_rates[63]);
  //sp 7
  sp_rates[6] -= (fwd_rates[67] - rev_rates[63]);

  //rxn 68
  //sp 18
  sp_rates[17] -= (fwd_rates[68] - rev_rates[64]);
  //sp 3
  sp_rates[2] += (fwd_rates[68] - rev_rates[64]);
  //sp 5
  sp_rates[4] -= (fwd_rates[68] - rev_rates[64]);
  //sp 15
  sp_rates[14] += (fwd_rates[68] - rev_rates[64]);

  //rxn 69
  //sp 12
  sp_rates[11] += fwd_rates[69];
  //sp 18
  sp_rates[17] -= fwd_rates[69];
  //sp 3
  sp_rates[2] += 2.0 * fwd_rates[69];
  //sp 4
  sp_rates[3] -= fwd_rates[69];

  //rxn 70
  //sp 12
  sp_rates[11] += (fwd_rates[70] - rev_rates[65]);
  //sp 18
  sp_rates[17] -= (fwd_rates[70] - rev_rates[65]);
  //sp 4
  sp_rates[3] -= (fwd_rates[70] - rev_rates[65]);
  //sp 6
  sp_rates[5] += (fwd_rates[70] - rev_rates[65]);

  //rxn 71
  //sp 17
  sp_rates[16] += (fwd_rates[71] - rev_rates[66]) * pres_mod[12];
  //sp 18
  sp_rates[17] -= (fwd_rates[71] - rev_rates[66]) * pres_mod[12];
  //sp 6
  sp_rates[5] += (fwd_rates[71] - rev_rates[66]) * pres_mod[12];

  //rxn 72
  //sp 18
  sp_rates[17] -= (fwd_rates[72] - rev_rates[67]) * pres_mod[13];
  //sp 3
  sp_rates[2] += (fwd_rates[72] - rev_rates[67]) * pres_mod[13];
  //sp 16
  sp_rates[15] += (fwd_rates[72] - rev_rates[67]) * pres_mod[13];

  //rxn 73
  //sp 18
  sp_rates[17] -= (fwd_rates[73] - rev_rates[68]);
  //sp 3
  sp_rates[2] += (fwd_rates[73] - rev_rates[68]);
  //sp 20
  sp_rates[19] += (fwd_rates[73] - rev_rates[68]);
  //sp 16
  sp_rates[15] -= (fwd_rates[73] - rev_rates[68]);

  //rxn 74
  //sp 18
  sp_rates[17] -= (fwd_rates[74] - rev_rates[69]);
  //sp 12
  sp_rates[11] += (fwd_rates[74] - rev_rates[69]);
  //sp 13
  sp_rates[12] -= (fwd_rates[74] - rev_rates[69]);
  //sp 15
  sp_rates[14] += (fwd_rates[74] - rev_rates[69]);

  //rxn 75
  //sp 18
  sp_rates[17] -= (fwd_rates[75] - rev_rates[70]);
  //sp 5
  sp_rates[4] += (fwd_rates[75] - rev_rates[70]);
  //sp 15
  sp_rates[14] += (fwd_rates[75] - rev_rates[70]);
  //sp 8
  sp_rates[7] -= (fwd_rates[75] - rev_rates[70]);

  //rxn 76
  //sp 18
  sp_rates[17] += (fwd_rates[76] - rev_rates[71]) * pres_mod[14];
  //sp 19
  sp_rates[18] = -(fwd_rates[76] - rev_rates[71]) * pres_mod[14];

  //rxn 77
  //sp 19
  sp_rates[18] -= (fwd_rates[77] - rev_rates[72]);
  //sp 3
  sp_rates[2] += (fwd_rates[77] - rev_rates[72]);
  //sp 21
  sp_rates[20] = (fwd_rates[77] - rev_rates[72]);
  //sp 6
  sp_rates[5] -= (fwd_rates[77] - rev_rates[72]);

  //rxn 78
  //sp 18
  sp_rates[17] += (fwd_rates[78] - rev_rates[73]);
  //sp 19
  sp_rates[18] -= (fwd_rates[78] - rev_rates[73]);

  //rxn 79
  //sp 3
  sp_rates[2] += (fwd_rates[79] - rev_rates[74]);
  //sp 5
  sp_rates[4] += (fwd_rates[79] - rev_rates[74]);
  //sp 7
  sp_rates[6] -= (fwd_rates[79] - rev_rates[74]);
  //sp 12
  sp_rates[11] += (fwd_rates[79] - rev_rates[74]);
  //sp 19
  sp_rates[18] -= (fwd_rates[79] - rev_rates[74]);

  //rxn 80
  //sp 18
  sp_rates[17] += (fwd_rates[80] - rev_rates[75]);
  //sp 19
  sp_rates[18] -= (fwd_rates[80] - rev_rates[75]);

  //rxn 81
  //sp 19
  sp_rates[18] -= (fwd_rates[81] - rev_rates[76]);
  //sp 3
  sp_rates[2] += (fwd_rates[81] - rev_rates[76]);
  //sp 5
  sp_rates[4] -= (fwd_rates[81] - rev_rates[76]);
  //sp 15
  sp_rates[14] += (fwd_rates[81] - rev_rates[76]);

  //rxn 82
  //sp 19
  sp_rates[18] -= (fwd_rates[82] - rev_rates[77]);
  //sp 12
  sp_rates[11] += (fwd_rates[82] - rev_rates[77]);
  //sp 13
  sp_rates[12] -= (fwd_rates[82] - rev_rates[77]);
  //sp 15
  sp_rates[14] += (fwd_rates[82] - rev_rates[77]);

  //rxn 83
  //sp 3
  sp_rates[2] -= (fwd_rates[83] - rev_rates[78]);
  //sp 14
  sp_rates[13] += (fwd_rates[83] - rev_rates[78]);
  //sp 6
  sp_rates[5] += (fwd_rates[83] - rev_rates[78]);
  //sp 15
  sp_rates[14] -= (fwd_rates[83] - rev_rates[78]);

  //rxn 84
  //sp 10
  sp_rates[9] += (fwd_rates[84] - rev_rates[79]);
  //sp 5
  sp_rates[4] -= (fwd_rates[84] - rev_rates[79]);
  //sp 14
  sp_rates[13] += (fwd_rates[84] - rev_rates[79]);
  //sp 15
  sp_rates[14] -= (fwd_rates[84] - rev_rates[79]);

  //rxn 85
  //sp 4
  sp_rates[3] -= (fwd_rates[85] - rev_rates[80]);
  //sp 5
  sp_rates[4] += (fwd_rates[85] - rev_rates[80]);
  //sp 14
  sp_rates[13] += (fwd_rates[85] - rev_rates[80]);
  //sp 15
  sp_rates[14] -= (fwd_rates[85] - rev_rates[80]);

  //rxn 86
  //sp 8
  sp_rates[7] += (fwd_rates[86] - rev_rates[81]);
  //sp 14
  sp_rates[13] += (fwd_rates[86] - rev_rates[81]);
  //sp 15
  sp_rates[14] -= (fwd_rates[86] - rev_rates[81]);
  //sp 7
  sp_rates[6] -= (fwd_rates[86] - rev_rates[81]);

  //rxn 87
  //sp 9
  sp_rates[8] += (fwd_rates[87] - rev_rates[82]);
  //sp 14
  sp_rates[13] += (fwd_rates[87] - rev_rates[82]);
  //sp 15
  sp_rates[14] -= (fwd_rates[87] - rev_rates[82]);
  //sp 8
  sp_rates[7] -= (fwd_rates[87] - rev_rates[82]);

  //rxn 88
  //sp 14
  sp_rates[13] += (fwd_rates[88] - rev_rates[83]);
  //sp 21
  sp_rates[20] -= (fwd_rates[88] - rev_rates[83]);
  //sp 22
  sp_rates[21] = (fwd_rates[88] - rev_rates[83]);
  //sp 15
  sp_rates[14] -= (fwd_rates[88] - rev_rates[83]);

  //rxn 89
  //sp 3
  sp_rates[2] += (fwd_rates[89] - rev_rates[84]) * pres_mod[15];
  //sp 14
  sp_rates[13] += (fwd_rates[89] - rev_rates[84]) * pres_mod[15];
  //sp 15
  sp_rates[14] -= (fwd_rates[89] - rev_rates[84]) * pres_mod[15];

  //rxn 90
  //sp 12
  sp_rates[11] += (fwd_rates[90] - rev_rates[85]) * pres_mod[16];
  //sp 6
  sp_rates[5] += (fwd_rates[90] - rev_rates[85]) * pres_mod[16];
  //sp 15
  sp_rates[14] -= (fwd_rates[90] - rev_rates[85]) * pres_mod[16];

  //rxn 91
  //sp 10
  sp_rates[9] += (fwd_rates[91] - rev_rates[86]) * pres_mod[17];
  //sp 12
  sp_rates[11] += (fwd_rates[91] - rev_rates[86]) * pres_mod[17];
  //sp 23
  sp_rates[22] = -(fwd_rates[91] - rev_rates[86]) * pres_mod[17];

  //rxn 92
  //sp 13
  sp_rates[12] += (fwd_rates[92] - rev_rates[87]) * pres_mod[18];
  //sp 6
  sp_rates[5] += (fwd_rates[92] - rev_rates[87]) * pres_mod[18];
  //sp 23
  sp_rates[22] -= (fwd_rates[92] - rev_rates[87]) * pres_mod[18];

  //rxn 93
  //sp 10
  sp_rates[9] += (fwd_rates[93] - rev_rates[88]);
  //sp 8
  sp_rates[7] -= (fwd_rates[93] - rev_rates[88]);
  //sp 23
  sp_rates[22] += (fwd_rates[93] - rev_rates[88]);
  //sp 24
  sp_rates[23] = -(fwd_rates[93] - rev_rates[88]);

  //rxn 94
  //sp 3
  sp_rates[2] += (fwd_rates[94] - rev_rates[89]);
  //sp 5
  sp_rates[4] -= (fwd_rates[94] - rev_rates[89]);
  //sp 15
  sp_rates[14] -= (fwd_rates[94] - rev_rates[89]);
  //sp 23
  sp_rates[22] += (fwd_rates[94] - rev_rates[89]);

  //rxn 95
  //sp 10
  sp_rates[9] += (fwd_rates[95] - rev_rates[90]) * pres_mod[19];
  //sp 12
  sp_rates[11] += (fwd_rates[95] - rev_rates[90]) * pres_mod[19];
  //sp 23
  sp_rates[22] -= (fwd_rates[95] - rev_rates[90]) * pres_mod[19];

  //rxn 96
  //sp 13
  sp_rates[12] += (fwd_rates[96] - rev_rates[91]) * pres_mod[20];
  //sp 6
  sp_rates[5] += (fwd_rates[96] - rev_rates[91]) * pres_mod[20];
  //sp 23
  sp_rates[22] -= (fwd_rates[96] - rev_rates[91]) * pres_mod[20];

  //rxn 97
  //sp 3
  sp_rates[2] -= (fwd_rates[97] - rev_rates[92]);
  //sp 5
  sp_rates[4] += (fwd_rates[97] - rev_rates[92]);
  //sp 6
  sp_rates[5] += (fwd_rates[97] - rev_rates[92]);
  //sp 12
  sp_rates[11] += (fwd_rates[97] - rev_rates[92]);
  //sp 23
  sp_rates[22] -= (fwd_rates[97] - rev_rates[92]);

  //rxn 98
  //sp 5
  sp_rates[4] += (fwd_rates[98] - rev_rates[93]);
  //sp 12
  sp_rates[11] += (fwd_rates[98] - rev_rates[93]);
  //sp 21
  sp_rates[20] -= (fwd_rates[98] - rev_rates[93]);
  //sp 22
  sp_rates[21] += (fwd_rates[98] - rev_rates[93]);
  //sp 23
  sp_rates[22] -= (fwd_rates[98] - rev_rates[93]);

  //rxn 99
  //sp 12
  sp_rates[11] += (fwd_rates[99] - rev_rates[94]);
  //sp 4
  sp_rates[3] -= (fwd_rates[99] - rev_rates[94]);
  //sp 5
  sp_rates[4] += 2.0 * (fwd_rates[99] - rev_rates[94]);
  //sp 23
  sp_rates[22] -= (fwd_rates[99] - rev_rates[94]);

  //rxn 100
  //sp 5
  sp_rates[4] += fwd_rates[100];
  //sp 8
  sp_rates[7] -= fwd_rates[100];
  //sp 9
  sp_rates[8] += fwd_rates[100];
  //sp 12
  sp_rates[11] += fwd_rates[100];
  //sp 23
  sp_rates[22] -= fwd_rates[100];

  //rxn 101
  //sp 25
  sp_rates[24] = -(fwd_rates[101] - rev_rates[95]);
  //sp 13
  sp_rates[12] += (fwd_rates[101] - rev_rates[95]);
  //sp 7
  sp_rates[6] -= (fwd_rates[101] - rev_rates[95]);
  //sp 8
  sp_rates[7] += (fwd_rates[101] - rev_rates[95]);

  //rxn 102
  //sp 18
  sp_rates[17] += (fwd_rates[102] - rev_rates[96]) * pres_mod[21];
  //sp 3
  sp_rates[2] += (fwd_rates[102] - rev_rates[96]) * pres_mod[21];
  //sp 21
  sp_rates[20] -= (fwd_rates[102] - rev_rates[96]) * pres_mod[21];

  //rxn 103
  //sp 21
  sp_rates[20] -= (fwd_rates[103] - rev_rates[97]) * pres_mod[22];
  //sp 6
  sp_rates[5] += (fwd_rates[103] - rev_rates[97]) * pres_mod[22];
  //sp 16
  sp_rates[15] += (fwd_rates[103] - rev_rates[97]) * pres_mod[22];

  //rxn 104
  //sp 3
  sp_rates[2] += (fwd_rates[104] - rev_rates[98]);
  //sp 4
  sp_rates[3] -= (fwd_rates[104] - rev_rates[98]);
  //sp 21
  sp_rates[20] -= (fwd_rates[104] - rev_rates[98]);
  //sp 15
  sp_rates[14] += (fwd_rates[104] - rev_rates[98]);

  //rxn 105
  //sp 3
  sp_rates[2] += fwd_rates[105];
  //sp 4
  sp_rates[3] -= fwd_rates[105];
  //sp 6
  sp_rates[5] += fwd_rates[105];
  //sp 12
  sp_rates[11] += fwd_rates[105];
  //sp 21
  sp_rates[20] -= fwd_rates[105];

  //rxn 106
  //sp 21
  sp_rates[20] -= (fwd_rates[106] - rev_rates[99]);
  //sp 5
  sp_rates[4] += (fwd_rates[106] - rev_rates[99]);
  //sp 7
  sp_rates[6] -= (fwd_rates[106] - rev_rates[99]);
  //sp 15
  sp_rates[14] += (fwd_rates[106] - rev_rates[99]);

  //rxn 107
  //sp 26
  sp_rates[25] = (fwd_rates[107] - rev_rates[100]);
  //sp 4
  sp_rates[3] += (fwd_rates[107] - rev_rates[100]);
  //sp 21
  sp_rates[20] -= (fwd_rates[107] - rev_rates[100]);
  //sp 7
  sp_rates[6] -= (fwd_rates[107] - rev_rates[100]);

  //rxn 108
  //sp 26
  sp_rates[25] += (fwd_rates[108] - rev_rates[101]);
  //sp 21
  sp_rates[20] -= (fwd_rates[108] - rev_rates[101]);
  //sp 5
  sp_rates[4] += (fwd_rates[108] - rev_rates[101]);
  //sp 8
  sp_rates[7] -= (fwd_rates[108] - rev_rates[101]);

  //rxn 109
  //sp 18
  sp_rates[17] -= (fwd_rates[109] - rev_rates[102]);
  //sp 3
  sp_rates[2] += (fwd_rates[109] - rev_rates[102]);
  //sp 21
  sp_rates[20] += (fwd_rates[109] - rev_rates[102]);
  //sp 6
  sp_rates[5] -= (fwd_rates[109] - rev_rates[102]);

  //rxn 110
  //sp 18
  sp_rates[17] -= (fwd_rates[110] - rev_rates[103]);
  //sp 27
  sp_rates[26] = (fwd_rates[110] - rev_rates[103]);
  //sp 3
  sp_rates[2] += (fwd_rates[110] - rev_rates[103]);
  //sp 21
  sp_rates[20] -= (fwd_rates[110] - rev_rates[103]);

  //rxn 111
  //sp 3
  sp_rates[2] -= (fwd_rates[111] - rev_rates[104]) * pres_mod[23];
  //sp 21
  sp_rates[20] -= (fwd_rates[111] - rev_rates[104]) * pres_mod[23];
  //sp 22
  sp_rates[21] += (fwd_rates[111] - rev_rates[104]) * pres_mod[23];

  //rxn 112
  //sp 27
  sp_rates[26] += (fwd_rates[112] - rev_rates[105]);
  //sp 21
  sp_rates[20] -= 2.0 * (fwd_rates[112] - rev_rates[105]);
  //sp 6
  sp_rates[5] += (fwd_rates[112] - rev_rates[105]);

  //rxn 113
  //sp 22
  sp_rates[21] += (fwd_rates[113] - rev_rates[106]);
  //sp 12
  sp_rates[11] += (fwd_rates[113] - rev_rates[106]);
  //sp 21
  sp_rates[20] -= (fwd_rates[113] - rev_rates[106]);
  //sp 14
  sp_rates[13] -= (fwd_rates[113] - rev_rates[106]);

  //rxn 114
  //sp 6
  sp_rates[5] += (fwd_rates[114] - rev_rates[107]);
  //sp 28
  sp_rates[27] = (fwd_rates[114] - rev_rates[107]);
  //sp 21
  sp_rates[20] -= (fwd_rates[114] - rev_rates[107]);
  //sp 14
  sp_rates[13] -= (fwd_rates[114] - rev_rates[107]);

  //rxn 115
  //sp 21
  sp_rates[20] -= 2.0 * (fwd_rates[115] - rev_rates[108]) * pres_mod[24];
  //sp 29
  sp_rates[28] = (fwd_rates[115] - rev_rates[108]) * pres_mod[24];

  //rxn 116
  //sp 10
  sp_rates[9] += (fwd_rates[116] - rev_rates[109]);
  //sp 19
  sp_rates[18] += (fwd_rates[116] - rev_rates[109]);
  //sp 21
  sp_rates[20] -= (fwd_rates[116] - rev_rates[109]);
  //sp 5
  sp_rates[4] -= (fwd_rates[116] - rev_rates[109]);

  //rxn 117
  //sp 6
  sp_rates[5] += (fwd_rates[117] - rev_rates[110]);
  //sp 21
  sp_rates[20] -= (fwd_rates[117] - rev_rates[110]);
  //sp 5
  sp_rates[4] -= (fwd_rates[117] - rev_rates[110]);
  //sp 15
  sp_rates[14] += (fwd_rates[117] - rev_rates[110]);

  //rxn 118
  //sp 3
  sp_rates[2] += (fwd_rates[118] - rev_rates[111]);
  //sp 21
  sp_rates[20] -= (fwd_rates[118] - rev_rates[111]);
  //sp 5
  sp_rates[4] -= (fwd_rates[118] - rev_rates[111]);
  //sp 24
  sp_rates[23] += (fwd_rates[118] - rev_rates[111]);

  //rxn 119
  //sp 26
  sp_rates[25] += (fwd_rates[119] - rev_rates[112]);
  //sp 3
  sp_rates[2] += (fwd_rates[119] - rev_rates[112]);
  //sp 21
  sp_rates[20] -= (fwd_rates[119] - rev_rates[112]);
  //sp 5
  sp_rates[4] -= (fwd_rates[119] - rev_rates[112]);

  //rxn 120
  //sp 18
  sp_rates[17] += (fwd_rates[120] - rev_rates[113]);
  //sp 21
  sp_rates[20] -= (fwd_rates[120] - rev_rates[113]);
  //sp 5
  sp_rates[4] -= (fwd_rates[120] - rev_rates[113]);
  //sp 10
  sp_rates[9] += (fwd_rates[120] - rev_rates[113]);

  //rxn 121
  //sp 26
  sp_rates[25] -= (fwd_rates[121] - rev_rates[114]) * pres_mod[25];
  //sp 24
  sp_rates[23] += (fwd_rates[121] - rev_rates[114]) * pres_mod[25];

  //rxn 122
  //sp 26
  sp_rates[25] -= (fwd_rates[122] - rev_rates[115]) * pres_mod[26];
  //sp 3
  sp_rates[2] += (fwd_rates[122] - rev_rates[115]) * pres_mod[26];
  //sp 15
  sp_rates[14] += (fwd_rates[122] - rev_rates[115]) * pres_mod[26];

  //rxn 123
  //sp 26
  sp_rates[25] -= (fwd_rates[123] - rev_rates[116]);
  //sp 3
  sp_rates[2] -= (fwd_rates[123] - rev_rates[116]);
  //sp 6
  sp_rates[5] += (fwd_rates[123] - rev_rates[116]);
  //sp 15
  sp_rates[14] += (fwd_rates[123] - rev_rates[116]);

  //rxn 124
  //sp 26
  sp_rates[25] -= (fwd_rates[124] - rev_rates[117]);
  //sp 4
  sp_rates[3] -= (fwd_rates[124] - rev_rates[117]);
  //sp 5
  sp_rates[4] += (fwd_rates[124] - rev_rates[117]);
  //sp 15
  sp_rates[14] += (fwd_rates[124] - rev_rates[117]);

  //rxn 125
  //sp 26
  sp_rates[25] -= (fwd_rates[125] - rev_rates[118]);
  //sp 5
  sp_rates[4] -= (fwd_rates[125] - rev_rates[118]);
  //sp 15
  sp_rates[14] += (fwd_rates[125] - rev_rates[118]);
  //sp 10
  sp_rates[9] += (fwd_rates[125] - rev_rates[118]);

  //rxn 126
  //sp 26
  sp_rates[25] -= (fwd_rates[126] - rev_rates[119]);
  //sp 8
  sp_rates[7] += (fwd_rates[126] - rev_rates[119]);
  //sp 7
  sp_rates[6] -= (fwd_rates[126] - rev_rates[119]);
  //sp 15
  sp_rates[14] += (fwd_rates[126] - rev_rates[119]);

  //rxn 127
  //sp 9
  sp_rates[8] += (fwd_rates[127] - rev_rates[120]);
  //sp 26
  sp_rates[25] -= (fwd_rates[127] - rev_rates[120]);
  //sp 15
  sp_rates[14] += (fwd_rates[127] - rev_rates[120]);
  //sp 8
  sp_rates[7] -= (fwd_rates[127] - rev_rates[120]);

  //rxn 128
  //sp 26
  sp_rates[25] -= (fwd_rates[128] - rev_rates[121]);
  //sp 12
  sp_rates[11] -= (fwd_rates[128] - rev_rates[121]);
  //sp 21
  sp_rates[20] += (fwd_rates[128] - rev_rates[121]);
  //sp 13
  sp_rates[12] += (fwd_rates[128] - rev_rates[121]);

  //rxn 129
  //sp 26
  sp_rates[25] -= (fwd_rates[129] - rev_rates[122]);
  //sp 21
  sp_rates[20] -= (fwd_rates[129] - rev_rates[122]);
  //sp 22
  sp_rates[21] += (fwd_rates[129] - rev_rates[122]);
  //sp 15
  sp_rates[14] += (fwd_rates[129] - rev_rates[122]);

  //rxn 130
  //sp 26
  sp_rates[25] -= 2.0 * (fwd_rates[130] - rev_rates[123]);
  //sp 30
  sp_rates[29] = (fwd_rates[130] - rev_rates[123]);
  //sp 15
  sp_rates[14] += (fwd_rates[130] - rev_rates[123]);

  //rxn 131
  //sp 26
  sp_rates[25] -= (fwd_rates[131] - rev_rates[124]);
  //sp 14
  sp_rates[13] += (fwd_rates[131] - rev_rates[124]);
  //sp 30
  sp_rates[29] += (fwd_rates[131] - rev_rates[124]);
  //sp 15
  sp_rates[14] -= (fwd_rates[131] - rev_rates[124]);

  //rxn 132
  //sp 3
  sp_rates[2] += (fwd_rates[132] - rev_rates[125]) * pres_mod[27];
  //sp 15
  sp_rates[14] += (fwd_rates[132] - rev_rates[125]) * pres_mod[27];
  //sp 24
  sp_rates[23] -= (fwd_rates[132] - rev_rates[125]) * pres_mod[27];

  //rxn 133
  //sp 3
  sp_rates[2] -= (fwd_rates[133] - rev_rates[126]);
  //sp 6
  sp_rates[5] += (fwd_rates[133] - rev_rates[126]);
  //sp 15
  sp_rates[14] += (fwd_rates[133] - rev_rates[126]);
  //sp 24
  sp_rates[23] -= (fwd_rates[133] - rev_rates[126]);

  //rxn 134
  //sp 15
  sp_rates[14] += (fwd_rates[134] - rev_rates[127]);
  //sp 8
  sp_rates[7] += (fwd_rates[134] - rev_rates[127]);
  //sp 7
  sp_rates[6] -= (fwd_rates[134] - rev_rates[127]);
  //sp 24
  sp_rates[23] -= (fwd_rates[134] - rev_rates[127]);

  //rxn 135
  //sp 15
  sp_rates[14] += (fwd_rates[135] - rev_rates[128]);
  //sp 8
  sp_rates[7] += (fwd_rates[135] - rev_rates[128]);
  //sp 7
  sp_rates[6] -= (fwd_rates[135] - rev_rates[128]);
  //sp 24
  sp_rates[23] -= (fwd_rates[135] - rev_rates[128]);

  //rxn 136
  //sp 4
  sp_rates[3] -= (fwd_rates[136] - rev_rates[129]);
  //sp 5
  sp_rates[4] += (fwd_rates[136] - rev_rates[129]);
  //sp 15
  sp_rates[14] += (fwd_rates[136] - rev_rates[129]);
  //sp 24
  sp_rates[23] -= (fwd_rates[136] - rev_rates[129]);

  //rxn 137
  //sp 10
  sp_rates[9] += (fwd_rates[137] - rev_rates[130]);
  //sp 5
  sp_rates[4] -= (fwd_rates[137] - rev_rates[130]);
  //sp 15
  sp_rates[14] += (fwd_rates[137] - rev_rates[130]);
  //sp 24
  sp_rates[23] -= (fwd_rates[137] - rev_rates[130]);

  //rxn 138
  //sp 9
  sp_rates[8] += (fwd_rates[138] - rev_rates[131]);
  //sp 8
  sp_rates[7] -= (fwd_rates[138] - rev_rates[131]);
  //sp 15
  sp_rates[14] += (fwd_rates[138] - rev_rates[131]);
  //sp 24
  sp_rates[23] -= (fwd_rates[138] - rev_rates[131]);

  //rxn 139
  //sp 30
  sp_rates[29] += (fwd_rates[139] - rev_rates[132]);
  //sp 12
  sp_rates[11] += (fwd_rates[139] - rev_rates[132]);
  //sp 14
  sp_rates[13] -= (fwd_rates[139] - rev_rates[132]);
  //sp 24
  sp_rates[23] -= (fwd_rates[139] - rev_rates[132]);

  //rxn 140
  //sp 14
  sp_rates[13] -= (fwd_rates[140] - rev_rates[133]);
  //sp 15
  sp_rates[14] += 2.0 * (fwd_rates[140] - rev_rates[133]);
  //sp 24
  sp_rates[23] -= (fwd_rates[140] - rev_rates[133]);

  //rxn 141
  //sp 30
  sp_rates[29] += (fwd_rates[141] - rev_rates[134]);
  //sp 15
  sp_rates[14] += (fwd_rates[141] - rev_rates[134]);
  //sp 24
  sp_rates[23] -= 2.0 * (fwd_rates[141] - rev_rates[134]);

  //rxn 142
  //sp 26
  sp_rates[25] -= (fwd_rates[142] - rev_rates[135]);
  //sp 30
  sp_rates[29] += (fwd_rates[142] - rev_rates[135]);
  //sp 15
  sp_rates[14] += (fwd_rates[142] - rev_rates[135]);
  //sp 24
  sp_rates[23] -= (fwd_rates[142] - rev_rates[135]);

  //rxn 143
  //sp 21
  sp_rates[20] -= (fwd_rates[143] - rev_rates[136]) * pres_mod[28];
  //sp 7
  sp_rates[6] -= (fwd_rates[143] - rev_rates[136]) * pres_mod[28];
  //sp 31
  sp_rates[30] = (fwd_rates[143] - rev_rates[136]) * pres_mod[28];

  //rxn 144
  //sp 26
  sp_rates[25] += (fwd_rates[144] - rev_rates[137]);
  //sp 3
  sp_rates[2] -= (fwd_rates[144] - rev_rates[137]);
  //sp 5
  sp_rates[4] += (fwd_rates[144] - rev_rates[137]);
  //sp 31
  sp_rates[30] -= (fwd_rates[144] - rev_rates[137]);

  //rxn 145
  //sp 26
  sp_rates[25] += (fwd_rates[145] - rev_rates[138]);
  //sp 4
  sp_rates[3] -= (fwd_rates[145] - rev_rates[138]);
  //sp 31
  sp_rates[30] -= (fwd_rates[145] - rev_rates[138]);
  //sp 7
  sp_rates[6] += (fwd_rates[145] - rev_rates[138]);

  //rxn 146
  //sp 26
  sp_rates[25] += (fwd_rates[146] - rev_rates[139]);
  //sp 5
  sp_rates[4] -= (fwd_rates[146] - rev_rates[139]);
  //sp 31
  sp_rates[30] -= (fwd_rates[146] - rev_rates[139]);
  //sp 8
  sp_rates[7] += (fwd_rates[146] - rev_rates[139]);

  //rxn 147
  //sp 32
  sp_rates[31] = (fwd_rates[147] - rev_rates[140]);
  //sp 7
  sp_rates[6] += (fwd_rates[147] - rev_rates[140]);
  //sp 31
  sp_rates[30] -= (fwd_rates[147] - rev_rates[140]);
  //sp 8
  sp_rates[7] -= (fwd_rates[147] - rev_rates[140]);

  //rxn 148
  //sp 26
  sp_rates[25] += 2.0 * (fwd_rates[148] - rev_rates[141]);
  //sp 21
  sp_rates[20] -= (fwd_rates[148] - rev_rates[141]);
  //sp 31
  sp_rates[30] -= (fwd_rates[148] - rev_rates[141]);

  //rxn 149
  //sp 21
  sp_rates[20] += (fwd_rates[149] - rev_rates[142]);
  //sp 22
  sp_rates[21] -= (fwd_rates[149] - rev_rates[142]);
  //sp 31
  sp_rates[30] -= (fwd_rates[149] - rev_rates[142]);
  //sp 32
  sp_rates[31] += (fwd_rates[149] - rev_rates[142]);

  //rxn 150
  //sp 26
  sp_rates[25] += (fwd_rates[150] - rev_rates[143]);
  //sp 12
  sp_rates[11] -= (fwd_rates[150] - rev_rates[143]);
  //sp 13
  sp_rates[12] += (fwd_rates[150] - rev_rates[143]);
  //sp 31
  sp_rates[30] -= (fwd_rates[150] - rev_rates[143]);

  //rxn 151
  //sp 3
  sp_rates[2] += (fwd_rates[151] - rev_rates[144]);
  //sp 13
  sp_rates[12] += (fwd_rates[151] - rev_rates[144]);
  //sp 14
  sp_rates[13] -= (fwd_rates[151] - rev_rates[144]);
  //sp 26
  sp_rates[25] += (fwd_rates[151] - rev_rates[144]);
  //sp 31
  sp_rates[30] -= (fwd_rates[151] - rev_rates[144]);

  //rxn 152
  //sp 32
  sp_rates[31] += (fwd_rates[152] - rev_rates[145]);
  //sp 14
  sp_rates[13] += (fwd_rates[152] - rev_rates[145]);
  //sp 15
  sp_rates[14] -= (fwd_rates[152] - rev_rates[145]);
  //sp 31
  sp_rates[30] -= (fwd_rates[152] - rev_rates[145]);

  //rxn 153
  //sp 3
  sp_rates[2] += (fwd_rates[153] - rev_rates[146]);
  //sp 12
  sp_rates[11] += (fwd_rates[153] - rev_rates[146]);
  //sp 15
  sp_rates[14] -= (fwd_rates[153] - rev_rates[146]);
  //sp 31
  sp_rates[30] -= (fwd_rates[153] - rev_rates[146]);
  //sp 32
  sp_rates[31] += (fwd_rates[153] - rev_rates[146]);

  //rxn 154
  //sp 15
  sp_rates[14] += (fwd_rates[154] - rev_rates[147]);
  //sp 32
  sp_rates[31] += (fwd_rates[154] - rev_rates[147]);
  //sp 31
  sp_rates[30] -= (fwd_rates[154] - rev_rates[147]);
  //sp 24
  sp_rates[23] -= (fwd_rates[154] - rev_rates[147]);

  //rxn 155
  //sp 26
  sp_rates[25] -= (fwd_rates[155] - rev_rates[148]);
  //sp 32
  sp_rates[31] += (fwd_rates[155] - rev_rates[148]);
  //sp 31
  sp_rates[30] -= (fwd_rates[155] - rev_rates[148]);
  //sp 15
  sp_rates[14] += (fwd_rates[155] - rev_rates[148]);

  //rxn 156
  //sp 32
  sp_rates[31] += (fwd_rates[156] - rev_rates[149]);
  //sp 30
  sp_rates[29] -= (fwd_rates[156] - rev_rates[149]);
  //sp 31
  sp_rates[30] -= (fwd_rates[156] - rev_rates[149]);
  //sp 24
  sp_rates[23] += (fwd_rates[156] - rev_rates[149]);

  //rxn 157
  //sp 26
  sp_rates[25] += 2.0 * (fwd_rates[157] - rev_rates[150]);
  //sp 31
  sp_rates[30] -= 2.0 * (fwd_rates[157] - rev_rates[150]);
  //sp 7
  sp_rates[6] += (fwd_rates[157] - rev_rates[150]);

  //rxn 158
  //sp 26
  sp_rates[25] += 2.0 * (fwd_rates[158] - rev_rates[151]);
  //sp 31
  sp_rates[30] -= 2.0 * (fwd_rates[158] - rev_rates[151]);
  //sp 7
  sp_rates[6] += (fwd_rates[158] - rev_rates[151]);

  //rxn 159
  //sp 7
  sp_rates[6] += (fwd_rates[159] - rev_rates[152]);
  //sp 30
  sp_rates[29] += (fwd_rates[159] - rev_rates[152]);
  //sp 31
  sp_rates[30] -= 2.0 * (fwd_rates[159] - rev_rates[152]);
  //sp 15
  sp_rates[14] += (fwd_rates[159] - rev_rates[152]);

  //rxn 160
  //sp 33
  sp_rates[32] = (fwd_rates[160] - rev_rates[153]);
  //sp 27
  sp_rates[26] -= (fwd_rates[160] - rev_rates[153]);
  //sp 31
  sp_rates[30] -= (fwd_rates[160] - rev_rates[153]);
  //sp 32
  sp_rates[31] += (fwd_rates[160] - rev_rates[153]);

  //rxn 161
  //sp 34
  sp_rates[33] = -(fwd_rates[161] - rev_rates[154]);
  //sp 36
  sp_rates[35] = (fwd_rates[161] - rev_rates[154]);
  //sp 31
  sp_rates[30] -= (fwd_rates[161] - rev_rates[154]);
  //sp 26
  sp_rates[25] += (fwd_rates[161] - rev_rates[154]);

  //rxn 162
  //sp 34
  sp_rates[33] += (fwd_rates[162] - rev_rates[155]);
  //sp 29
  sp_rates[28] -= (fwd_rates[162] - rev_rates[155]);
  //sp 31
  sp_rates[30] -= (fwd_rates[162] - rev_rates[155]);
  //sp 32
  sp_rates[31] += (fwd_rates[162] - rev_rates[155]);

  //rxn 163
  //sp 37
  sp_rates[36] = -(fwd_rates[163] - rev_rates[156]);
  //sp 38
  sp_rates[37] = (fwd_rates[163] - rev_rates[156]);
  //sp 31
  sp_rates[30] -= (fwd_rates[163] - rev_rates[156]);
  //sp 32
  sp_rates[31] += (fwd_rates[163] - rev_rates[156]);

  //rxn 164
  //sp 3
  sp_rates[2] -= (fwd_rates[164] - rev_rates[157]);
  //sp 6
  sp_rates[5] += (fwd_rates[164] - rev_rates[157]);
  //sp 21
  sp_rates[20] += (fwd_rates[164] - rev_rates[157]);
  //sp 22
  sp_rates[21] -= (fwd_rates[164] - rev_rates[157]);

  //rxn 165
  //sp 5
  sp_rates[4] += (fwd_rates[165] - rev_rates[158]);
  //sp 4
  sp_rates[3] -= (fwd_rates[165] - rev_rates[158]);
  //sp 21
  sp_rates[20] += (fwd_rates[165] - rev_rates[158]);
  //sp 22
  sp_rates[21] -= (fwd_rates[165] - rev_rates[158]);

  //rxn 166
  //sp 10
  sp_rates[9] += (fwd_rates[166] - rev_rates[159]);
  //sp 21
  sp_rates[20] += (fwd_rates[166] - rev_rates[159]);
  //sp 5
  sp_rates[4] -= (fwd_rates[166] - rev_rates[159]);
  //sp 22
  sp_rates[21] -= (fwd_rates[166] - rev_rates[159]);

  //rxn 167
  //sp 21
  sp_rates[20] += (fwd_rates[167] - rev_rates[160]);
  //sp 22
  sp_rates[21] -= (fwd_rates[167] - rev_rates[160]);
  //sp 7
  sp_rates[6] -= (fwd_rates[167] - rev_rates[160]);
  //sp 8
  sp_rates[7] += (fwd_rates[167] - rev_rates[160]);

  //rxn 168
  //sp 9
  sp_rates[8] += (fwd_rates[168] - rev_rates[161]);
  //sp 21
  sp_rates[20] += (fwd_rates[168] - rev_rates[161]);
  //sp 22
  sp_rates[21] -= (fwd_rates[168] - rev_rates[161]);
  //sp 8
  sp_rates[7] -= (fwd_rates[168] - rev_rates[161]);

  //rxn 169
  //sp 27
  sp_rates[26] += (fwd_rates[169] - rev_rates[162]);
  //sp 3
  sp_rates[2] += (fwd_rates[169] - rev_rates[162]);
  //sp 22
  sp_rates[21] -= (fwd_rates[169] - rev_rates[162]);
  //sp 16
  sp_rates[15] -= (fwd_rates[169] - rev_rates[162]);

  //rxn 170
  //sp 18
  sp_rates[17] -= (fwd_rates[170] - rev_rates[163]);
  //sp 21
  sp_rates[20] += 2.0 * (fwd_rates[170] - rev_rates[163]);
  //sp 22
  sp_rates[21] -= (fwd_rates[170] - rev_rates[163]);

  //rxn 171
  //sp 19
  sp_rates[18] -= (fwd_rates[171] - rev_rates[164]);
  //sp 21
  sp_rates[20] += 2.0 * (fwd_rates[171] - rev_rates[164]);
  //sp 22
  sp_rates[21] -= (fwd_rates[171] - rev_rates[164]);

  //rxn 172
  //sp 18
  sp_rates[17] += (fwd_rates[172] - rev_rates[165]);
  //sp 19
  sp_rates[18] -= (fwd_rates[172] - rev_rates[165]);

  //rxn 173
  //sp 26
  sp_rates[25] -= (fwd_rates[173] - rev_rates[166]) * pres_mod[29];
  //sp 3
  sp_rates[2] -= (fwd_rates[173] - rev_rates[166]) * pres_mod[29];
  //sp 30
  sp_rates[29] += (fwd_rates[173] - rev_rates[166]) * pres_mod[29];

  //rxn 174
  //sp 5
  sp_rates[4] += (fwd_rates[174] - rev_rates[167]) * pres_mod[30];
  //sp 21
  sp_rates[20] += (fwd_rates[174] - rev_rates[167]) * pres_mod[30];
  //sp 30
  sp_rates[29] -= (fwd_rates[174] - rev_rates[167]) * pres_mod[30];

  //rxn 175
  //sp 10
  sp_rates[9] += (fwd_rates[175] - rev_rates[168]) * pres_mod[31];
  //sp 19
  sp_rates[18] += (fwd_rates[175] - rev_rates[168]) * pres_mod[31];
  //sp 30
  sp_rates[29] -= (fwd_rates[175] - rev_rates[168]) * pres_mod[31];

  //rxn 176
  //sp 3
  sp_rates[2] += (fwd_rates[176] - rev_rates[169]) * pres_mod[32];
  //sp 30
  sp_rates[29] -= (fwd_rates[176] - rev_rates[169]) * pres_mod[32];
  //sp 24
  sp_rates[23] += (fwd_rates[176] - rev_rates[169]) * pres_mod[32];

  //rxn 177
  //sp 3
  sp_rates[2] -= (fwd_rates[177] - rev_rates[170]);
  //sp 6
  sp_rates[5] += (fwd_rates[177] - rev_rates[170]);
  //sp 30
  sp_rates[29] -= (fwd_rates[177] - rev_rates[170]);
  //sp 24
  sp_rates[23] += (fwd_rates[177] - rev_rates[170]);

  //rxn 178
  //sp 26
  sp_rates[25] += (fwd_rates[178] - rev_rates[171]);
  //sp 3
  sp_rates[2] -= (fwd_rates[178] - rev_rates[171]);
  //sp 6
  sp_rates[5] += (fwd_rates[178] - rev_rates[171]);
  //sp 30
  sp_rates[29] -= (fwd_rates[178] - rev_rates[171]);

  //rxn 179
  //sp 4
  sp_rates[3] -= (fwd_rates[179] - rev_rates[172]);
  //sp 5
  sp_rates[4] += (fwd_rates[179] - rev_rates[172]);
  //sp 30
  sp_rates[29] -= (fwd_rates[179] - rev_rates[172]);
  //sp 24
  sp_rates[23] += (fwd_rates[179] - rev_rates[172]);

  //rxn 180
  //sp 26
  sp_rates[25] += (fwd_rates[180] - rev_rates[173]);
  //sp 4
  sp_rates[3] -= (fwd_rates[180] - rev_rates[173]);
  //sp 5
  sp_rates[4] += (fwd_rates[180] - rev_rates[173]);
  //sp 30
  sp_rates[29] -= (fwd_rates[180] - rev_rates[173]);

  //rxn 181
  //sp 10
  sp_rates[9] += (fwd_rates[181] - rev_rates[174]);
  //sp 5
  sp_rates[4] -= (fwd_rates[181] - rev_rates[174]);
  //sp 30
  sp_rates[29] -= (fwd_rates[181] - rev_rates[174]);
  //sp 24
  sp_rates[23] += (fwd_rates[181] - rev_rates[174]);

  //rxn 182
  //sp 26
  sp_rates[25] += (fwd_rates[182] - rev_rates[175]);
  //sp 5
  sp_rates[4] -= (fwd_rates[182] - rev_rates[175]);
  //sp 30
  sp_rates[29] -= (fwd_rates[182] - rev_rates[175]);
  //sp 10
  sp_rates[9] += (fwd_rates[182] - rev_rates[175]);

  //rxn 183
  //sp 8
  sp_rates[7] += (fwd_rates[183] - rev_rates[176]);
  //sp 30
  sp_rates[29] -= (fwd_rates[183] - rev_rates[176]);
  //sp 7
  sp_rates[6] -= (fwd_rates[183] - rev_rates[176]);
  //sp 24
  sp_rates[23] += (fwd_rates[183] - rev_rates[176]);

  //rxn 184
  //sp 26
  sp_rates[25] += (fwd_rates[184] - rev_rates[177]);
  //sp 30
  sp_rates[29] -= (fwd_rates[184] - rev_rates[177]);
  //sp 7
  sp_rates[6] -= (fwd_rates[184] - rev_rates[177]);
  //sp 8
  sp_rates[7] += (fwd_rates[184] - rev_rates[177]);

  //rxn 185
  //sp 14
  sp_rates[13] -= (fwd_rates[185] - rev_rates[178]);
  //sp 30
  sp_rates[29] -= (fwd_rates[185] - rev_rates[178]);
  //sp 15
  sp_rates[14] += (fwd_rates[185] - rev_rates[178]);
  //sp 24
  sp_rates[23] += (fwd_rates[185] - rev_rates[178]);

  //rxn 186
  //sp 9
  sp_rates[8] += (fwd_rates[186] - rev_rates[179]);
  //sp 24
  sp_rates[23] += (fwd_rates[186] - rev_rates[179]);
  //sp 30
  sp_rates[29] -= (fwd_rates[186] - rev_rates[179]);
  //sp 8
  sp_rates[7] -= (fwd_rates[186] - rev_rates[179]);

  //rxn 187
  //sp 9
  sp_rates[8] += (fwd_rates[187] - rev_rates[180]);
  //sp 26
  sp_rates[25] += (fwd_rates[187] - rev_rates[180]);
  //sp 30
  sp_rates[29] -= (fwd_rates[187] - rev_rates[180]);
  //sp 8
  sp_rates[7] -= (fwd_rates[187] - rev_rates[180]);

  //rxn 188
  //sp 22
  sp_rates[21] += (fwd_rates[188] - rev_rates[181]);
  //sp 21
  sp_rates[20] -= (fwd_rates[188] - rev_rates[181]);
  //sp 30
  sp_rates[29] -= (fwd_rates[188] - rev_rates[181]);
  //sp 24
  sp_rates[23] += (fwd_rates[188] - rev_rates[181]);

  //rxn 189
  //sp 26
  sp_rates[25] += (fwd_rates[189] - rev_rates[182]);
  //sp 22
  sp_rates[21] += (fwd_rates[189] - rev_rates[182]);
  //sp 21
  sp_rates[20] -= (fwd_rates[189] - rev_rates[182]);
  //sp 30
  sp_rates[29] -= (fwd_rates[189] - rev_rates[182]);

  //rxn 190
  //sp 26
  sp_rates[25] -= (fwd_rates[190] - rev_rates[183]);
  //sp 24
  sp_rates[23] += (fwd_rates[190] - rev_rates[183]);

  //rxn 191
  //sp 26
  sp_rates[25] += (fwd_rates[191] - rev_rates[184]) * pres_mod[33];
  //sp 5
  sp_rates[4] += (fwd_rates[191] - rev_rates[184]) * pres_mod[33];
  //sp 32
  sp_rates[31] -= (fwd_rates[191] - rev_rates[184]) * pres_mod[33];

  //rxn 192
  //sp 3
  sp_rates[2] -= (fwd_rates[192] - rev_rates[185]);
  //sp 6
  sp_rates[5] += (fwd_rates[192] - rev_rates[185]);
  //sp 31
  sp_rates[30] += (fwd_rates[192] - rev_rates[185]);
  //sp 32
  sp_rates[31] -= (fwd_rates[192] - rev_rates[185]);

  //rxn 193
  //sp 3
  sp_rates[2] -= (fwd_rates[193] - rev_rates[186]);
  //sp 6
  sp_rates[5] += (fwd_rates[193] - rev_rates[186]);
  //sp 39
  sp_rates[38] = (fwd_rates[193] - rev_rates[186]);
  //sp 32
  sp_rates[31] -= (fwd_rates[193] - rev_rates[186]);

  //rxn 194
  //sp 26
  sp_rates[25] += (fwd_rates[194] - rev_rates[187]);
  //sp 3
  sp_rates[2] -= (fwd_rates[194] - rev_rates[187]);
  //sp 10
  sp_rates[9] += (fwd_rates[194] - rev_rates[187]);
  //sp 32
  sp_rates[31] -= (fwd_rates[194] - rev_rates[187]);

  //rxn 195
  //sp 10
  sp_rates[9] += (fwd_rates[195] - rev_rates[188]);
  //sp 5
  sp_rates[4] -= (fwd_rates[195] - rev_rates[188]);
  //sp 31
  sp_rates[30] += (fwd_rates[195] - rev_rates[188]);
  //sp 32
  sp_rates[31] -= (fwd_rates[195] - rev_rates[188]);

  //rxn 196
  //sp 10
  sp_rates[9] += (fwd_rates[196] - rev_rates[189]);
  //sp 5
  sp_rates[4] -= (fwd_rates[196] - rev_rates[189]);
  //sp 39
  sp_rates[38] += (fwd_rates[196] - rev_rates[189]);
  //sp 32
  sp_rates[31] -= (fwd_rates[196] - rev_rates[189]);

  //rxn 197
  //sp 4
  sp_rates[3] -= (fwd_rates[197] - rev_rates[190]);
  //sp 5
  sp_rates[4] += (fwd_rates[197] - rev_rates[190]);
  //sp 39
  sp_rates[38] += (fwd_rates[197] - rev_rates[190]);
  //sp 32
  sp_rates[31] -= (fwd_rates[197] - rev_rates[190]);

  //rxn 198
  //sp 4
  sp_rates[3] -= (fwd_rates[198] - rev_rates[191]);
  //sp 5
  sp_rates[4] += (fwd_rates[198] - rev_rates[191]);
  //sp 31
  sp_rates[30] += (fwd_rates[198] - rev_rates[191]);
  //sp 32
  sp_rates[31] -= (fwd_rates[198] - rev_rates[191]);

  //rxn 199
  //sp 9
  sp_rates[8] += (fwd_rates[199] - rev_rates[192]);
  //sp 8
  sp_rates[7] -= (fwd_rates[199] - rev_rates[192]);
  //sp 31
  sp_rates[30] += (fwd_rates[199] - rev_rates[192]);
  //sp 32
  sp_rates[31] -= (fwd_rates[199] - rev_rates[192]);

  //rxn 200
  //sp 5
  sp_rates[4] += (fwd_rates[200] - rev_rates[193]);
  //sp 39
  sp_rates[38] -= (fwd_rates[200] - rev_rates[193]);
  //sp 15
  sp_rates[14] += (fwd_rates[200] - rev_rates[193]);

  //rxn 201
  //sp 41
  sp_rates[40] = (fwd_rates[201] - rev_rates[194]);
  //sp 4
  sp_rates[3] -= (fwd_rates[201] - rev_rates[194]);
  //sp 5
  sp_rates[4] += (fwd_rates[201] - rev_rates[194]);
  //sp 40
  sp_rates[39] = -(fwd_rates[201] - rev_rates[194]);

  //rxn 202
  //sp 41
  sp_rates[40] += (fwd_rates[202] - rev_rates[195]);
  //sp 3
  sp_rates[2] -= (fwd_rates[202] - rev_rates[195]);
  //sp 6
  sp_rates[5] += (fwd_rates[202] - rev_rates[195]);
  //sp 40
  sp_rates[39] -= (fwd_rates[202] - rev_rates[195]);

  //rxn 203
  //sp 41
  sp_rates[40] -= (fwd_rates[203] - rev_rates[196]);
  //sp 42
  sp_rates[41] = (fwd_rates[203] - rev_rates[196]);
  //sp 3
  sp_rates[2] += (fwd_rates[203] - rev_rates[196]);
  //sp 5
  sp_rates[4] -= (fwd_rates[203] - rev_rates[196]);

  //rxn 204
  //sp 41
  sp_rates[40] -= (fwd_rates[204] - rev_rates[197]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[204] - rev_rates[197]);
  //sp 7
  sp_rates[6] -= (fwd_rates[204] - rev_rates[197]);

  //rxn 205
  //sp 41
  sp_rates[40] -= (fwd_rates[205] - rev_rates[198]);
  //sp 12
  sp_rates[11] += (fwd_rates[205] - rev_rates[198]);
  //sp 5
  sp_rates[4] -= (fwd_rates[205] - rev_rates[198]);
  //sp 16
  sp_rates[15] += (fwd_rates[205] - rev_rates[198]);

  //rxn 206
  //sp 41
  sp_rates[40] -= (fwd_rates[206] - rev_rates[199]);
  //sp 17
  sp_rates[16] += (fwd_rates[206] - rev_rates[199]);
  //sp 4
  sp_rates[3] -= (fwd_rates[206] - rev_rates[199]);
  //sp 12
  sp_rates[11] += (fwd_rates[206] - rev_rates[199]);

  //rxn 207
  //sp 41
  sp_rates[40] -= (fwd_rates[207] - rev_rates[200]);
  //sp 21
  sp_rates[20] += (fwd_rates[207] - rev_rates[200]);
  //sp 22
  sp_rates[21] -= (fwd_rates[207] - rev_rates[200]);
  //sp 40
  sp_rates[39] += (fwd_rates[207] - rev_rates[200]);

  //rxn 208
  //sp 41
  sp_rates[40] -= (fwd_rates[208] - rev_rates[201]);
  //sp 20
  sp_rates[19] -= (fwd_rates[208] - rev_rates[201]);
  //sp 40
  sp_rates[39] += 2.0 * (fwd_rates[208] - rev_rates[201]);

  //rxn 209
  //sp 41
  sp_rates[40] -= (fwd_rates[209] - rev_rates[202]);
  //sp 27
  sp_rates[26] -= (fwd_rates[209] - rev_rates[202]);
  //sp 33
  sp_rates[32] += (fwd_rates[209] - rev_rates[202]);
  //sp 40
  sp_rates[39] += (fwd_rates[209] - rev_rates[202]);

  //rxn 210
  //sp 41
  sp_rates[40] -= (fwd_rates[210] - rev_rates[203]);
  //sp 34
  sp_rates[33] += (fwd_rates[210] - rev_rates[203]);
  //sp 29
  sp_rates[28] -= (fwd_rates[210] - rev_rates[203]);
  //sp 40
  sp_rates[39] += (fwd_rates[210] - rev_rates[203]);

  //rxn 211
  //sp 42
  sp_rates[41] += (fwd_rates[211] - rev_rates[204]);
  //sp 43
  sp_rates[42] = -(fwd_rates[211] - rev_rates[204]);
  //sp 5
  sp_rates[4] -= (fwd_rates[211] - rev_rates[204]);
  //sp 10
  sp_rates[9] += (fwd_rates[211] - rev_rates[204]);

  //rxn 212
  //sp 42
  sp_rates[41] -= (fwd_rates[212] - rev_rates[205]);
  //sp 3
  sp_rates[2] -= (fwd_rates[212] - rev_rates[205]);
  //sp 12
  sp_rates[11] += (fwd_rates[212] - rev_rates[205]);
  //sp 16
  sp_rates[15] += (fwd_rates[212] - rev_rates[205]);

  //rxn 213
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[213] - rev_rates[206]);
  //sp 42
  sp_rates[41] -= (fwd_rates[213] - rev_rates[206]);
  //sp 4
  sp_rates[3] -= (fwd_rates[213] - rev_rates[206]);

  //rxn 214
  //sp 42
  sp_rates[41] -= (fwd_rates[214] - rev_rates[207]);
  //sp 3
  sp_rates[2] += (fwd_rates[214] - rev_rates[207]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[214] - rev_rates[207]);
  //sp 5
  sp_rates[4] -= (fwd_rates[214] - rev_rates[207]);

  //rxn 215
  //sp 4
  sp_rates[3] += (fwd_rates[215] - rev_rates[208]);
  //sp 42
  sp_rates[41] -= (fwd_rates[215] - rev_rates[208]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[215] - rev_rates[208]);
  //sp 7
  sp_rates[6] -= (fwd_rates[215] - rev_rates[208]);

  //rxn 216
  //sp 42
  sp_rates[41] -= (fwd_rates[216] - rev_rates[209]);
  //sp 12
  sp_rates[11] += (fwd_rates[216] - rev_rates[209]);
  //sp 13
  sp_rates[12] += (fwd_rates[216] - rev_rates[209]);
  //sp 7
  sp_rates[6] -= (fwd_rates[216] - rev_rates[209]);

  //rxn 217
  //sp 17
  sp_rates[16] -= (fwd_rates[217] - rev_rates[210]);
  //sp 42
  sp_rates[41] -= (fwd_rates[217] - rev_rates[210]);
  //sp 41
  sp_rates[40] += (fwd_rates[217] - rev_rates[210]);
  //sp 12
  sp_rates[11] += (fwd_rates[217] - rev_rates[210]);

  //rxn 218
  //sp 12
  sp_rates[11] += (fwd_rates[218] - rev_rates[211]);
  //sp 4
  sp_rates[3] -= (fwd_rates[218] - rev_rates[211]);
  //sp 16
  sp_rates[15] += (fwd_rates[218] - rev_rates[211]);
  //sp 40
  sp_rates[39] -= (fwd_rates[218] - rev_rates[211]);

  //rxn 219
  //sp 43
  sp_rates[42] += (fwd_rates[219] - rev_rates[212]);
  //sp 4
  sp_rates[3] += (fwd_rates[219] - rev_rates[212]);
  //sp 7
  sp_rates[6] -= (fwd_rates[219] - rev_rates[212]);
  //sp 40
  sp_rates[39] -= (fwd_rates[219] - rev_rates[212]);

  //rxn 220
  //sp 3
  sp_rates[2] += (fwd_rates[220] - rev_rates[213]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[220] - rev_rates[213]);
  //sp 7
  sp_rates[6] -= (fwd_rates[220] - rev_rates[213]);
  //sp 40
  sp_rates[39] -= (fwd_rates[220] - rev_rates[213]);

  //rxn 221
  //sp 12
  sp_rates[11] += (fwd_rates[221] - rev_rates[214]);
  //sp 14
  sp_rates[13] += (fwd_rates[221] - rev_rates[214]);
  //sp 7
  sp_rates[6] -= (fwd_rates[221] - rev_rates[214]);
  //sp 40
  sp_rates[39] -= (fwd_rates[221] - rev_rates[214]);

  //rxn 222
  //sp 3
  sp_rates[2] += (fwd_rates[222] - rev_rates[215]);
  //sp 21
  sp_rates[20] -= (fwd_rates[222] - rev_rates[215]);
  //sp 45
  sp_rates[44] = (fwd_rates[222] - rev_rates[215]);
  //sp 40
  sp_rates[39] -= (fwd_rates[222] - rev_rates[215]);

  //rxn 223
  //sp 3
  sp_rates[2] += (fwd_rates[223] - rev_rates[216]);
  //sp 5
  sp_rates[4] -= (fwd_rates[223] - rev_rates[216]);
  //sp 43
  sp_rates[42] += (fwd_rates[223] - rev_rates[216]);
  //sp 40
  sp_rates[39] -= (fwd_rates[223] - rev_rates[216]);

  //rxn 224
  //sp 43
  sp_rates[42] += (fwd_rates[224] - rev_rates[217]) * pres_mod[34];
  //sp 12
  sp_rates[11] -= (fwd_rates[224] - rev_rates[217]) * pres_mod[34];
  //sp 16
  sp_rates[15] -= (fwd_rates[224] - rev_rates[217]) * pres_mod[34];

  //rxn 225
  //sp 18
  sp_rates[17] += (fwd_rates[225] - rev_rates[218]);
  //sp 3
  sp_rates[2] -= (fwd_rates[225] - rev_rates[218]);
  //sp 12
  sp_rates[11] += (fwd_rates[225] - rev_rates[218]);
  //sp 43
  sp_rates[42] -= (fwd_rates[225] - rev_rates[218]);

  //rxn 226
  //sp 12
  sp_rates[11] += (fwd_rates[226] - rev_rates[219]);
  //sp 3
  sp_rates[2] -= (fwd_rates[226] - rev_rates[219]);
  //sp 19
  sp_rates[18] += (fwd_rates[226] - rev_rates[219]);
  //sp 43
  sp_rates[42] -= (fwd_rates[226] - rev_rates[219]);

  //rxn 227
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[227] - rev_rates[220]);
  //sp 3
  sp_rates[2] += (fwd_rates[227] - rev_rates[220]);
  //sp 43
  sp_rates[42] -= (fwd_rates[227] - rev_rates[220]);
  //sp 4
  sp_rates[3] -= (fwd_rates[227] - rev_rates[220]);

  //rxn 228
  //sp 3
  sp_rates[2] += (fwd_rates[228] - rev_rates[221]);
  //sp 5
  sp_rates[4] -= (fwd_rates[228] - rev_rates[221]);
  //sp 43
  sp_rates[42] -= (fwd_rates[228] - rev_rates[221]);
  //sp 12
  sp_rates[11] += (fwd_rates[228] - rev_rates[221]);
  //sp 14
  sp_rates[13] += (fwd_rates[228] - rev_rates[221]);

  //rxn 229
  //sp 12
  sp_rates[11] += (fwd_rates[229] - rev_rates[222]);
  //sp 43
  sp_rates[42] -= (fwd_rates[229] - rev_rates[222]);
  //sp 20
  sp_rates[19] += (fwd_rates[229] - rev_rates[222]);
  //sp 16
  sp_rates[15] -= (fwd_rates[229] - rev_rates[222]);

  //rxn 230
  //sp 33
  sp_rates[32] += (fwd_rates[230] - rev_rates[223]);
  //sp 19
  sp_rates[18] -= (fwd_rates[230] - rev_rates[223]);
  //sp 12
  sp_rates[11] += (fwd_rates[230] - rev_rates[223]);
  //sp 43
  sp_rates[42] -= (fwd_rates[230] - rev_rates[223]);

  //rxn 231
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[231] - rev_rates[224]);
  //sp 43
  sp_rates[42] -= 2.0 * (fwd_rates[231] - rev_rates[224]);
  //sp 20
  sp_rates[19] += (fwd_rates[231] - rev_rates[224]);

  //rxn 232
  //sp 43
  sp_rates[42] -= (fwd_rates[232] - rev_rates[225]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[232] - rev_rates[225]);
  //sp 5
  sp_rates[4] += (fwd_rates[232] - rev_rates[225]);
  //sp 7
  sp_rates[6] -= (fwd_rates[232] - rev_rates[225]);

  //rxn 233
  //sp 3
  sp_rates[2] += (fwd_rates[233] - rev_rates[226]);
  //sp 7
  sp_rates[6] -= (fwd_rates[233] - rev_rates[226]);
  //sp 43
  sp_rates[42] -= (fwd_rates[233] - rev_rates[226]);
  //sp 12
  sp_rates[11] += (fwd_rates[233] - rev_rates[226]);
  //sp 13
  sp_rates[12] += (fwd_rates[233] - rev_rates[226]);

  //rxn 234
  //sp 4
  sp_rates[3] += (fwd_rates[234] - rev_rates[227]);
  //sp 7
  sp_rates[6] -= (fwd_rates[234] - rev_rates[227]);
  //sp 43
  sp_rates[42] -= (fwd_rates[234] - rev_rates[227]);
  //sp 12
  sp_rates[11] += (fwd_rates[234] - rev_rates[227]);
  //sp 14
  sp_rates[13] += (fwd_rates[234] - rev_rates[227]);

  //rxn 235
  //sp 43
  sp_rates[42] -= (fwd_rates[235] - rev_rates[228]);
  //sp 4
  sp_rates[3] -= (fwd_rates[235] - rev_rates[228]);
  //sp 13
  sp_rates[12] += (fwd_rates[235] - rev_rates[228]);
  //sp 16
  sp_rates[15] += (fwd_rates[235] - rev_rates[228]);

  //rxn 236
  //sp 12
  sp_rates[11] += (fwd_rates[236] - rev_rates[229]);
  //sp 43
  sp_rates[42] -= (fwd_rates[236] - rev_rates[229]);
  //sp 27
  sp_rates[26] += (fwd_rates[236] - rev_rates[229]);
  //sp 21
  sp_rates[20] -= (fwd_rates[236] - rev_rates[229]);

  //rxn 237
  //sp 12
  sp_rates[11] += (fwd_rates[237] - rev_rates[230]);
  //sp 43
  sp_rates[42] -= (fwd_rates[237] - rev_rates[230]);
  //sp 20
  sp_rates[19] -= (fwd_rates[237] - rev_rates[230]);
  //sp 45
  sp_rates[44] += (fwd_rates[237] - rev_rates[230]);

  //rxn 238
  //sp 3
  sp_rates[2] -= (fwd_rates[238] - rev_rates[231]) * pres_mod[35];
  //sp 20
  sp_rates[19] += (fwd_rates[238] - rev_rates[231]) * pres_mod[35];
  //sp 40
  sp_rates[39] -= (fwd_rates[238] - rev_rates[231]) * pres_mod[35];

  //rxn 239
  //sp 3
  sp_rates[2] -= (fwd_rates[239] - rev_rates[232]);
  //sp 20
  sp_rates[19] -= (fwd_rates[239] - rev_rates[232]);
  //sp 6
  sp_rates[5] += (fwd_rates[239] - rev_rates[232]);
  //sp 40
  sp_rates[39] += (fwd_rates[239] - rev_rates[232]);

  //rxn 240
  //sp 4
  sp_rates[3] -= (fwd_rates[240] - rev_rates[233]);
  //sp 3
  sp_rates[2] += (fwd_rates[240] - rev_rates[233]);
  //sp 20
  sp_rates[19] -= (fwd_rates[240] - rev_rates[233]);
  //sp 43
  sp_rates[42] += (fwd_rates[240] - rev_rates[233]);

  //rxn 241
  //sp 4
  sp_rates[3] -= (fwd_rates[241] - rev_rates[234]);
  //sp 18
  sp_rates[17] += (fwd_rates[241] - rev_rates[234]);
  //sp 20
  sp_rates[19] -= (fwd_rates[241] - rev_rates[234]);
  //sp 12
  sp_rates[11] += (fwd_rates[241] - rev_rates[234]);

  //rxn 242
  //sp 4
  sp_rates[3] -= (fwd_rates[242] - rev_rates[235]);
  //sp 20
  sp_rates[19] -= (fwd_rates[242] - rev_rates[235]);
  //sp 5
  sp_rates[4] += (fwd_rates[242] - rev_rates[235]);
  //sp 40
  sp_rates[39] += (fwd_rates[242] - rev_rates[235]);

  //rxn 243
  //sp 44
  sp_rates[43] = (fwd_rates[243] - rev_rates[236]) * pres_mod[36];
  //sp 20
  sp_rates[19] -= (fwd_rates[243] - rev_rates[236]) * pres_mod[36];

  //rxn 244
  //sp 3
  sp_rates[2] += (fwd_rates[244] - rev_rates[237]);
  //sp 20
  sp_rates[19] -= (fwd_rates[244] - rev_rates[237]);
  //sp 46
  sp_rates[45] = (fwd_rates[244] - rev_rates[237]);
  //sp 40
  sp_rates[39] -= (fwd_rates[244] - rev_rates[237]);

  //rxn 245
  //sp 20
  sp_rates[19] -= (fwd_rates[245] - rev_rates[238]) * pres_mod[37];
  //sp 21
  sp_rates[20] -= (fwd_rates[245] - rev_rates[238]) * pres_mod[37];
  //sp 47
  sp_rates[46] = (fwd_rates[245] - rev_rates[238]) * pres_mod[37];

  //rxn 246
  //sp 6
  sp_rates[5] += (fwd_rates[246] - rev_rates[239]);
  //sp 20
  sp_rates[19] -= 2.0 * (fwd_rates[246] - rev_rates[239]);
  //sp 46
  sp_rates[45] += (fwd_rates[246] - rev_rates[239]);

  //rxn 247
  //sp 33
  sp_rates[32] += (fwd_rates[247] - rev_rates[240]);
  //sp 12
  sp_rates[11] += (fwd_rates[247] - rev_rates[240]);
  //sp 20
  sp_rates[19] -= (fwd_rates[247] - rev_rates[240]);
  //sp 14
  sp_rates[13] -= (fwd_rates[247] - rev_rates[240]);

  //rxn 248
  //sp 4
  sp_rates[3] += (fwd_rates[248] - rev_rates[241]);
  //sp 20
  sp_rates[19] -= (fwd_rates[248] - rev_rates[241]);
  //sp 48
  sp_rates[47] = (fwd_rates[248] - rev_rates[241]);
  //sp 8
  sp_rates[7] -= (fwd_rates[248] - rev_rates[241]);

  //rxn 249
  //sp 20
  sp_rates[19] -= (fwd_rates[249] - rev_rates[242]);
  //sp 14
  sp_rates[13] += (fwd_rates[249] - rev_rates[242]);
  //sp 15
  sp_rates[14] += (fwd_rates[249] - rev_rates[242]);
  //sp 8
  sp_rates[7] -= (fwd_rates[249] - rev_rates[242]);

  //rxn 250
  //sp 49
  sp_rates[48] = (fwd_rates[250] - rev_rates[243]);
  //sp 20
  sp_rates[19] -= (fwd_rates[250] - rev_rates[243]);
  //sp 8
  sp_rates[7] -= (fwd_rates[250] - rev_rates[243]);

  //rxn 251
  //sp 20
  sp_rates[19] += (fwd_rates[251] - rev_rates[244]);
  //sp 21
  sp_rates[20] += (fwd_rates[251] - rev_rates[244]);
  //sp 22
  sp_rates[21] -= (fwd_rates[251] - rev_rates[244]);
  //sp 40
  sp_rates[39] -= (fwd_rates[251] - rev_rates[244]);

  //rxn 252
  //sp 3
  sp_rates[2] += (fwd_rates[252] - rev_rates[245]);
  //sp 19
  sp_rates[18] -= (fwd_rates[252] - rev_rates[245]);
  //sp 20
  sp_rates[19] -= (fwd_rates[252] - rev_rates[245]);
  //sp 45
  sp_rates[44] += (fwd_rates[252] - rev_rates[245]);

  //rxn 253
  //sp 18
  sp_rates[17] += (fwd_rates[253] - rev_rates[246]);
  //sp 19
  sp_rates[18] -= (fwd_rates[253] - rev_rates[246]);

  //rxn 254
  //sp 18
  sp_rates[17] -= (fwd_rates[254] - rev_rates[247]);
  //sp 3
  sp_rates[2] += (fwd_rates[254] - rev_rates[247]);
  //sp 20
  sp_rates[19] -= (fwd_rates[254] - rev_rates[247]);
  //sp 45
  sp_rates[44] += (fwd_rates[254] - rev_rates[247]);

  //rxn 255
  //sp 50
  sp_rates[49] = -(fwd_rates[255] - rev_rates[248]);
  //sp 3
  sp_rates[2] -= (fwd_rates[255] - rev_rates[248]);
  //sp 43
  sp_rates[42] += (fwd_rates[255] - rev_rates[248]);
  //sp 6
  sp_rates[5] += (fwd_rates[255] - rev_rates[248]);

  //rxn 256
  //sp 50
  sp_rates[49] -= (fwd_rates[256] - rev_rates[249]);
  //sp 43
  sp_rates[42] += (fwd_rates[256] - rev_rates[249]);
  //sp 4
  sp_rates[3] -= (fwd_rates[256] - rev_rates[249]);
  //sp 5
  sp_rates[4] += (fwd_rates[256] - rev_rates[249]);

  //rxn 257
  //sp 50
  sp_rates[49] -= (fwd_rates[257] - rev_rates[250]);
  //sp 43
  sp_rates[42] += (fwd_rates[257] - rev_rates[250]);
  //sp 5
  sp_rates[4] -= (fwd_rates[257] - rev_rates[250]);
  //sp 10
  sp_rates[9] += (fwd_rates[257] - rev_rates[250]);

  //rxn 258
  //sp 20
  sp_rates[19] += (fwd_rates[258] - rev_rates[251]);
  //sp 44
  sp_rates[43] -= (fwd_rates[258] - rev_rates[251]);

  //rxn 259
  //sp 28
  sp_rates[27] += (fwd_rates[259] - rev_rates[252]);
  //sp 3
  sp_rates[2] += (fwd_rates[259] - rev_rates[252]);
  //sp 44
  sp_rates[43] -= (fwd_rates[259] - rev_rates[252]);
  //sp 5
  sp_rates[4] -= (fwd_rates[259] - rev_rates[252]);

  //rxn 260
  //sp 44
  sp_rates[43] -= (fwd_rates[260] - rev_rates[253]);
  //sp 14
  sp_rates[13] += 2.0 * (fwd_rates[260] - rev_rates[253]);
  //sp 7
  sp_rates[6] -= (fwd_rates[260] - rev_rates[253]);

  //rxn 261
  //sp 44
  sp_rates[43] -= (fwd_rates[261] - rev_rates[254]) * pres_mod[38];
  //sp 51
  sp_rates[50] = (fwd_rates[261] - rev_rates[254]) * pres_mod[38];
  //sp 20
  sp_rates[19] -= (fwd_rates[261] - rev_rates[254]) * pres_mod[38];

  //rxn 262
  //sp 4
  sp_rates[3] -= (fwd_rates[262] - rev_rates[255]);
  //sp 18
  sp_rates[17] += (fwd_rates[262] - rev_rates[255]);
  //sp 44
  sp_rates[43] -= (fwd_rates[262] - rev_rates[255]);
  //sp 12
  sp_rates[11] += (fwd_rates[262] - rev_rates[255]);

  //rxn 263
  //sp 10
  sp_rates[9] += (fwd_rates[263] - rev_rates[256]);
  //sp 20
  sp_rates[19] -= (fwd_rates[263] - rev_rates[256]);
  //sp 5
  sp_rates[4] -= (fwd_rates[263] - rev_rates[256]);
  //sp 40
  sp_rates[39] += (fwd_rates[263] - rev_rates[256]);

  //rxn 264
  //sp 50
  sp_rates[49] += (fwd_rates[264] - rev_rates[257]);
  //sp 3
  sp_rates[2] += (fwd_rates[264] - rev_rates[257]);
  //sp 20
  sp_rates[19] -= (fwd_rates[264] - rev_rates[257]);
  //sp 5
  sp_rates[4] -= (fwd_rates[264] - rev_rates[257]);

  //rxn 265
  //sp 28
  sp_rates[27] += (fwd_rates[265] - rev_rates[258]);
  //sp 3
  sp_rates[2] += (fwd_rates[265] - rev_rates[258]);
  //sp 20
  sp_rates[19] -= (fwd_rates[265] - rev_rates[258]);
  //sp 5
  sp_rates[4] -= (fwd_rates[265] - rev_rates[258]);

  //rxn 266
  //sp 12
  sp_rates[11] += (fwd_rates[266] - rev_rates[259]);
  //sp 20
  sp_rates[19] -= (fwd_rates[266] - rev_rates[259]);
  //sp 5
  sp_rates[4] -= (fwd_rates[266] - rev_rates[259]);
  //sp 21
  sp_rates[20] += (fwd_rates[266] - rev_rates[259]);

  //rxn 267
  //sp 52
  sp_rates[51] = (fwd_rates[267] - rev_rates[260]);
  //sp 20
  sp_rates[19] -= (fwd_rates[267] - rev_rates[260]);
  //sp 5
  sp_rates[4] -= (fwd_rates[267] - rev_rates[260]);

  //rxn 268
  //sp 52
  sp_rates[51] += (fwd_rates[268] - rev_rates[261]) * pres_mod[39];
  //sp 20
  sp_rates[19] -= (fwd_rates[268] - rev_rates[261]) * pres_mod[39];
  //sp 5
  sp_rates[4] -= (fwd_rates[268] - rev_rates[261]) * pres_mod[39];

  //rxn 269
  //sp 28
  sp_rates[27] += (fwd_rates[269] - rev_rates[262]);
  //sp 3
  sp_rates[2] += (fwd_rates[269] - rev_rates[262]);
  //sp 52
  sp_rates[51] -= (fwd_rates[269] - rev_rates[262]);

  //rxn 270
  //sp 28
  sp_rates[27] += (fwd_rates[270] - rev_rates[263]);
  //sp 3
  sp_rates[2] -= (fwd_rates[270] - rev_rates[263]);
  //sp 52
  sp_rates[51] -= (fwd_rates[270] - rev_rates[263]);
  //sp 6
  sp_rates[5] += (fwd_rates[270] - rev_rates[263]);

  //rxn 271
  //sp 4
  sp_rates[3] -= (fwd_rates[271] - rev_rates[264]);
  //sp 28
  sp_rates[27] += (fwd_rates[271] - rev_rates[264]);
  //sp 52
  sp_rates[51] -= (fwd_rates[271] - rev_rates[264]);
  //sp 5
  sp_rates[4] += (fwd_rates[271] - rev_rates[264]);

  //rxn 272
  //sp 28
  sp_rates[27] += (fwd_rates[272] - rev_rates[265]);
  //sp 10
  sp_rates[9] += (fwd_rates[272] - rev_rates[265]);
  //sp 52
  sp_rates[51] -= (fwd_rates[272] - rev_rates[265]);
  //sp 5
  sp_rates[4] -= (fwd_rates[272] - rev_rates[265]);

  //rxn 273
  //sp 6
  sp_rates[5] += fwd_rates[273];
  //sp 7
  sp_rates[6] -= fwd_rates[273];
  //sp 13
  sp_rates[12] += fwd_rates[273];
  //sp 14
  sp_rates[13] += fwd_rates[273];
  //sp 52
  sp_rates[51] -= fwd_rates[273];

  //rxn 274
  //sp 50
  sp_rates[49] += (fwd_rates[274] - rev_rates[266]);
  //sp 3
  sp_rates[2] += (fwd_rates[274] - rev_rates[266]);
  //sp 52
  sp_rates[51] -= (fwd_rates[274] - rev_rates[266]);

  //rxn 275
  //sp 52
  sp_rates[51] -= (fwd_rates[275] - rev_rates[267]);
  //sp 48
  sp_rates[47] += (fwd_rates[275] - rev_rates[267]);

  //rxn 276
  //sp 52
  sp_rates[51] -= (fwd_rates[276] - rev_rates[268]);
  //sp 14
  sp_rates[13] += (fwd_rates[276] - rev_rates[268]);
  //sp 7
  sp_rates[6] -= (fwd_rates[276] - rev_rates[268]);
  //sp 23
  sp_rates[22] += (fwd_rates[276] - rev_rates[268]);

  //rxn 277
  //sp 50
  sp_rates[49] += (fwd_rates[277] - rev_rates[269]);
  //sp 52
  sp_rates[51] -= (fwd_rates[277] - rev_rates[269]);
  //sp 21
  sp_rates[20] -= (fwd_rates[277] - rev_rates[269]);
  //sp 22
  sp_rates[21] += (fwd_rates[277] - rev_rates[269]);

  //rxn 278
  //sp 28
  sp_rates[27] += (fwd_rates[278] - rev_rates[270]) * pres_mod[40];
  //sp 18
  sp_rates[17] -= (fwd_rates[278] - rev_rates[270]) * pres_mod[40];
  //sp 12
  sp_rates[11] -= (fwd_rates[278] - rev_rates[270]) * pres_mod[40];

  //rxn 279
  //sp 12
  sp_rates[11] += (fwd_rates[279] - rev_rates[271]);
  //sp 3
  sp_rates[2] -= (fwd_rates[279] - rev_rates[271]);
  //sp 28
  sp_rates[27] -= (fwd_rates[279] - rev_rates[271]);
  //sp 21
  sp_rates[20] += (fwd_rates[279] - rev_rates[271]);

  //rxn 280
  //sp 3
  sp_rates[2] -= (fwd_rates[280] - rev_rates[272]);
  //sp 28
  sp_rates[27] -= (fwd_rates[280] - rev_rates[272]);
  //sp 43
  sp_rates[42] += (fwd_rates[280] - rev_rates[272]);
  //sp 6
  sp_rates[5] += (fwd_rates[280] - rev_rates[272]);

  //rxn 281
  //sp 4
  sp_rates[3] -= (fwd_rates[281] - rev_rates[273]);
  //sp 18
  sp_rates[17] += (fwd_rates[281] - rev_rates[273]);
  //sp 28
  sp_rates[27] -= (fwd_rates[281] - rev_rates[273]);
  //sp 13
  sp_rates[12] += (fwd_rates[281] - rev_rates[273]);

  //rxn 282
  //sp 4
  sp_rates[3] -= (fwd_rates[282] - rev_rates[274]);
  //sp 28
  sp_rates[27] -= (fwd_rates[282] - rev_rates[274]);
  //sp 14
  sp_rates[13] += 2.0 * (fwd_rates[282] - rev_rates[274]);

  //rxn 283
  //sp 4
  sp_rates[3] -= (fwd_rates[283] - rev_rates[275]);
  //sp 43
  sp_rates[42] += (fwd_rates[283] - rev_rates[275]);
  //sp 28
  sp_rates[27] -= (fwd_rates[283] - rev_rates[275]);
  //sp 5
  sp_rates[4] += (fwd_rates[283] - rev_rates[275]);

  //rxn 284
  //sp 12
  sp_rates[11] += (fwd_rates[284] - rev_rates[276]);
  //sp 28
  sp_rates[27] -= (fwd_rates[284] - rev_rates[276]);
  //sp 5
  sp_rates[4] -= (fwd_rates[284] - rev_rates[276]);
  //sp 24
  sp_rates[23] += (fwd_rates[284] - rev_rates[276]);

  //rxn 285
  //sp 13
  sp_rates[12] += (fwd_rates[285] - rev_rates[277]);
  //sp 28
  sp_rates[27] -= (fwd_rates[285] - rev_rates[277]);
  //sp 5
  sp_rates[4] -= (fwd_rates[285] - rev_rates[277]);
  //sp 21
  sp_rates[20] += (fwd_rates[285] - rev_rates[277]);

  //rxn 286
  //sp 10
  sp_rates[9] += (fwd_rates[286] - rev_rates[278]);
  //sp 43
  sp_rates[42] += (fwd_rates[286] - rev_rates[278]);
  //sp 28
  sp_rates[27] -= (fwd_rates[286] - rev_rates[278]);
  //sp 5
  sp_rates[4] -= (fwd_rates[286] - rev_rates[278]);

  //rxn 287
  //sp 12
  sp_rates[11] += (fwd_rates[287] - rev_rates[279]);
  //sp 34
  sp_rates[33] += (fwd_rates[287] - rev_rates[279]);
  //sp 28
  sp_rates[27] -= (fwd_rates[287] - rev_rates[279]);
  //sp 21
  sp_rates[20] -= (fwd_rates[287] - rev_rates[279]);

  //rxn 288
  //sp 3
  sp_rates[2] += (fwd_rates[288] - rev_rates[280]);
  //sp 28
  sp_rates[27] += (fwd_rates[288] - rev_rates[280]);
  //sp 15
  sp_rates[14] -= (fwd_rates[288] - rev_rates[280]);
  //sp 16
  sp_rates[15] -= (fwd_rates[288] - rev_rates[280]);

  //rxn 289
  //sp 27
  sp_rates[26] += (fwd_rates[289] - rev_rates[281]);
  //sp 12
  sp_rates[11] += (fwd_rates[289] - rev_rates[281]);
  //sp 19
  sp_rates[18] -= (fwd_rates[289] - rev_rates[281]);
  //sp 28
  sp_rates[27] -= (fwd_rates[289] - rev_rates[281]);

  //rxn 290
  //sp 33
  sp_rates[32] += (fwd_rates[290] - rev_rates[282]) * pres_mod[41];
  //sp 3
  sp_rates[2] -= (fwd_rates[290] - rev_rates[282]) * pres_mod[41];
  //sp 20
  sp_rates[19] -= (fwd_rates[290] - rev_rates[282]) * pres_mod[41];

  //rxn 291
  //sp 33
  sp_rates[32] -= (fwd_rates[291] - rev_rates[283]);
  //sp 3
  sp_rates[2] -= (fwd_rates[291] - rev_rates[283]);
  //sp 20
  sp_rates[19] += (fwd_rates[291] - rev_rates[283]);
  //sp 6
  sp_rates[5] += (fwd_rates[291] - rev_rates[283]);

  //rxn 292
  //sp 33
  sp_rates[32] -= (fwd_rates[292] - rev_rates[284]);
  //sp 3
  sp_rates[2] -= (fwd_rates[292] - rev_rates[284]);
  //sp 44
  sp_rates[43] += (fwd_rates[292] - rev_rates[284]);
  //sp 6
  sp_rates[5] += (fwd_rates[292] - rev_rates[284]);

  //rxn 293
  //sp 33
  sp_rates[32] -= (fwd_rates[293] - rev_rates[285]);
  //sp 10
  sp_rates[9] += (fwd_rates[293] - rev_rates[285]);
  //sp 20
  sp_rates[19] += (fwd_rates[293] - rev_rates[285]);
  //sp 5
  sp_rates[4] -= (fwd_rates[293] - rev_rates[285]);

  //rxn 294
  //sp 33
  sp_rates[32] -= (fwd_rates[294] - rev_rates[286]);
  //sp 20
  sp_rates[19] += (fwd_rates[294] - rev_rates[286]);
  //sp 4
  sp_rates[3] -= (fwd_rates[294] - rev_rates[286]);
  //sp 5
  sp_rates[4] += (fwd_rates[294] - rev_rates[286]);

  //rxn 295
  //sp 33
  sp_rates[32] -= (fwd_rates[295] - rev_rates[287]);
  //sp 12
  sp_rates[11] += (fwd_rates[295] - rev_rates[287]);
  //sp 4
  sp_rates[3] -= (fwd_rates[295] - rev_rates[287]);
  //sp 21
  sp_rates[20] += (fwd_rates[295] - rev_rates[287]);

  //rxn 296
  //sp 33
  sp_rates[32] -= (fwd_rates[296] - rev_rates[288]);
  //sp 18
  sp_rates[17] += (fwd_rates[296] - rev_rates[288]);
  //sp 4
  sp_rates[3] -= (fwd_rates[296] - rev_rates[288]);
  //sp 14
  sp_rates[13] += (fwd_rates[296] - rev_rates[288]);

  //rxn 297
  //sp 33
  sp_rates[32] -= (fwd_rates[297] - rev_rates[289]);
  //sp 20
  sp_rates[19] += (fwd_rates[297] - rev_rates[289]);
  //sp 21
  sp_rates[20] -= (fwd_rates[297] - rev_rates[289]);
  //sp 22
  sp_rates[21] += (fwd_rates[297] - rev_rates[289]);

  //rxn 298
  //sp 33
  sp_rates[32] -= (fwd_rates[298] - rev_rates[290]);
  //sp 20
  sp_rates[19] += 2.0 * (fwd_rates[298] - rev_rates[290]);
  //sp 40
  sp_rates[39] -= (fwd_rates[298] - rev_rates[290]);

  //rxn 299
  //sp 33
  sp_rates[32] -= (fwd_rates[299] - rev_rates[291]);
  //sp 18
  sp_rates[17] += (fwd_rates[299] - rev_rates[291]);
  //sp 20
  sp_rates[19] += (fwd_rates[299] - rev_rates[291]);
  //sp 16
  sp_rates[15] -= (fwd_rates[299] - rev_rates[291]);

  //rxn 300
  //sp 33
  sp_rates[32] -= (fwd_rates[300] - rev_rates[292]);
  //sp 5
  sp_rates[4] += (fwd_rates[300] - rev_rates[292]);
  //sp 48
  sp_rates[47] += (fwd_rates[300] - rev_rates[292]);
  //sp 8
  sp_rates[7] -= (fwd_rates[300] - rev_rates[292]);

  //rxn 301
  //sp 33
  sp_rates[32] -= (fwd_rates[301] - rev_rates[293]);
  //sp 34
  sp_rates[33] += (fwd_rates[301] - rev_rates[293]);
  //sp 27
  sp_rates[26] += (fwd_rates[301] - rev_rates[293]);
  //sp 29
  sp_rates[28] -= (fwd_rates[301] - rev_rates[293]);

  //rxn 302
  //sp 33
  sp_rates[32] -= (fwd_rates[302] - rev_rates[294]);
  //sp 27
  sp_rates[26] += (fwd_rates[302] - rev_rates[294]);
  //sp 12
  sp_rates[11] += (fwd_rates[302] - rev_rates[294]);
  //sp 14
  sp_rates[13] -= (fwd_rates[302] - rev_rates[294]);

  //rxn 303
  //sp 33
  sp_rates[32] -= (fwd_rates[303] - rev_rates[295]);
  //sp 27
  sp_rates[26] += (fwd_rates[303] - rev_rates[295]);
  //sp 14
  sp_rates[13] += (fwd_rates[303] - rev_rates[295]);
  //sp 15
  sp_rates[14] -= (fwd_rates[303] - rev_rates[295]);

  //rxn 304
  //sp 33
  sp_rates[32] -= (fwd_rates[304] - rev_rates[296]);
  //sp 3
  sp_rates[2] += (fwd_rates[304] - rev_rates[296]);
  //sp 53
  sp_rates[52] = (fwd_rates[304] - rev_rates[296]);
  //sp 15
  sp_rates[14] -= (fwd_rates[304] - rev_rates[296]);

  //rxn 305
  //sp 33
  sp_rates[32] -= (fwd_rates[305] - rev_rates[297]);
  //sp 54
  sp_rates[53] = (fwd_rates[305] - rev_rates[297]);
  //sp 14
  sp_rates[13] -= (fwd_rates[305] - rev_rates[297]);

  //rxn 306
  //sp 33
  sp_rates[32] -= (fwd_rates[306] - rev_rates[298]);
  //sp 55
  sp_rates[54] = (fwd_rates[306] - rev_rates[298]);
  //sp 5
  sp_rates[4] += (fwd_rates[306] - rev_rates[298]);
  //sp 7
  sp_rates[6] -= (fwd_rates[306] - rev_rates[298]);

  //rxn 307
  //sp 33
  sp_rates[32] -= (fwd_rates[307] - rev_rates[299]);
  //sp 55
  sp_rates[54] += (fwd_rates[307] - rev_rates[299]);
  //sp 5
  sp_rates[4] += (fwd_rates[307] - rev_rates[299]);
  //sp 7
  sp_rates[6] -= (fwd_rates[307] - rev_rates[299]);

  //rxn 308
  //sp 33
  sp_rates[32] -= (fwd_rates[308] - rev_rates[300]);
  //sp 20
  sp_rates[19] += (fwd_rates[308] - rev_rates[300]);
  //sp 7
  sp_rates[6] -= (fwd_rates[308] - rev_rates[300]);
  //sp 8
  sp_rates[7] += (fwd_rates[308] - rev_rates[300]);

  //rxn 309
  //sp 33
  sp_rates[32] -= (fwd_rates[309] - rev_rates[301]);
  //sp 20
  sp_rates[19] += (fwd_rates[309] - rev_rates[301]);
  //sp 7
  sp_rates[6] -= (fwd_rates[309] - rev_rates[301]);
  //sp 8
  sp_rates[7] += (fwd_rates[309] - rev_rates[301]);

  //rxn 310
  //sp 33
  sp_rates[32] -= (fwd_rates[310] - rev_rates[302]);
  //sp 3
  sp_rates[2] += (fwd_rates[310] - rev_rates[302]);
  //sp 7
  sp_rates[6] -= (fwd_rates[310] - rev_rates[302]);
  //sp 56
  sp_rates[55] = (fwd_rates[310] - rev_rates[302]);

  //rxn 311
  //sp 33
  sp_rates[32] -= (fwd_rates[311] - rev_rates[303]);
  //sp 3
  sp_rates[2] += (fwd_rates[311] - rev_rates[303]);
  //sp 7
  sp_rates[6] -= (fwd_rates[311] - rev_rates[303]);
  //sp 56
  sp_rates[55] += (fwd_rates[311] - rev_rates[303]);

  //rxn 312
  //sp 33
  sp_rates[32] -= (fwd_rates[312] - rev_rates[304]);
  //sp 28
  sp_rates[27] += (fwd_rates[312] - rev_rates[304]);
  //sp 5
  sp_rates[4] += (fwd_rates[312] - rev_rates[304]);
  //sp 7
  sp_rates[6] -= (fwd_rates[312] - rev_rates[304]);

  //rxn 313
  //sp 33
  sp_rates[32] -= (fwd_rates[313] - rev_rates[305]);
  //sp 28
  sp_rates[27] += (fwd_rates[313] - rev_rates[305]);
  //sp 5
  sp_rates[4] += (fwd_rates[313] - rev_rates[305]);
  //sp 7
  sp_rates[6] -= (fwd_rates[313] - rev_rates[305]);

  //rxn 314
  //sp 33
  sp_rates[32] -= (fwd_rates[314] - rev_rates[306]);
  //sp 26
  sp_rates[25] += (fwd_rates[314] - rev_rates[306]);
  //sp 12
  sp_rates[11] += (fwd_rates[314] - rev_rates[306]);
  //sp 7
  sp_rates[6] -= (fwd_rates[314] - rev_rates[306]);

  //rxn 315
  //sp 33
  sp_rates[32] -= (fwd_rates[315] - rev_rates[307]);
  //sp 26
  sp_rates[25] += (fwd_rates[315] - rev_rates[307]);
  //sp 12
  sp_rates[11] += (fwd_rates[315] - rev_rates[307]);
  //sp 7
  sp_rates[6] -= (fwd_rates[315] - rev_rates[307]);

  //rxn 316
  //sp 33
  sp_rates[32] -= (fwd_rates[316] - rev_rates[308]);
  //sp 21
  sp_rates[20] += (fwd_rates[316] - rev_rates[308]);
  //sp 13
  sp_rates[12] += (fwd_rates[316] - rev_rates[308]);
  //sp 7
  sp_rates[6] -= (fwd_rates[316] - rev_rates[308]);

  //rxn 317
  //sp 33
  sp_rates[32] -= (fwd_rates[317] - rev_rates[309]);
  //sp 21
  sp_rates[20] += (fwd_rates[317] - rev_rates[309]);
  //sp 13
  sp_rates[12] += (fwd_rates[317] - rev_rates[309]);
  //sp 7
  sp_rates[6] -= (fwd_rates[317] - rev_rates[309]);

  //rxn 318
  //sp 49
  sp_rates[48] -= (fwd_rates[318] - rev_rates[310]);
  //sp 5
  sp_rates[4] += (fwd_rates[318] - rev_rates[310]);
  //sp 55
  sp_rates[54] += (fwd_rates[318] - rev_rates[310]);

  //rxn 319
  //sp 49
  sp_rates[48] -= (fwd_rates[319] - rev_rates[311]);
  //sp 5
  sp_rates[4] += (fwd_rates[319] - rev_rates[311]);
  //sp 55
  sp_rates[54] += (fwd_rates[319] - rev_rates[311]);

  //rxn 320
  //sp 49
  sp_rates[48] -= (fwd_rates[320] - rev_rates[312]);
  //sp 4
  sp_rates[3] += (fwd_rates[320] - rev_rates[312]);
  //sp 48
  sp_rates[47] += (fwd_rates[320] - rev_rates[312]);

  //rxn 321
  //sp 49
  sp_rates[48] -= (fwd_rates[321] - rev_rates[313]);
  //sp 4
  sp_rates[3] += (fwd_rates[321] - rev_rates[313]);
  //sp 48
  sp_rates[47] += (fwd_rates[321] - rev_rates[313]);

  //rxn 322
  //sp 49
  sp_rates[48] -= (fwd_rates[322] - rev_rates[314]);
  //sp 3
  sp_rates[2] += (fwd_rates[322] - rev_rates[314]);
  //sp 56
  sp_rates[55] += (fwd_rates[322] - rev_rates[314]);

  //rxn 323
  //sp 49
  sp_rates[48] -= (fwd_rates[323] - rev_rates[315]);
  //sp 3
  sp_rates[2] += (fwd_rates[323] - rev_rates[315]);
  //sp 56
  sp_rates[55] += (fwd_rates[323] - rev_rates[315]);

  //rxn 324
  //sp 49
  sp_rates[48] -= (fwd_rates[324] - rev_rates[316]);
  //sp 28
  sp_rates[27] += (fwd_rates[324] - rev_rates[316]);
  //sp 5
  sp_rates[4] += (fwd_rates[324] - rev_rates[316]);

  //rxn 325
  //sp 49
  sp_rates[48] -= (fwd_rates[325] - rev_rates[317]);
  //sp 28
  sp_rates[27] += (fwd_rates[325] - rev_rates[317]);
  //sp 5
  sp_rates[4] += (fwd_rates[325] - rev_rates[317]);

  //rxn 326
  //sp 49
  sp_rates[48] -= (fwd_rates[326] - rev_rates[318]);
  //sp 14
  sp_rates[13] += (fwd_rates[326] - rev_rates[318]);
  //sp 15
  sp_rates[14] += (fwd_rates[326] - rev_rates[318]);

  //rxn 327
  //sp 49
  sp_rates[48] -= (fwd_rates[327] - rev_rates[319]);
  //sp 14
  sp_rates[13] += (fwd_rates[327] - rev_rates[319]);
  //sp 15
  sp_rates[14] += (fwd_rates[327] - rev_rates[319]);

  //rxn 328
  //sp 49
  sp_rates[48] -= fwd_rates[328];
  //sp 3
  sp_rates[2] += fwd_rates[328];
  //sp 12
  sp_rates[11] += fwd_rates[328];
  //sp 15
  sp_rates[14] += fwd_rates[328];

  //rxn 329
  //sp 49
  sp_rates[48] -= fwd_rates[329];
  //sp 3
  sp_rates[2] += fwd_rates[329];
  //sp 12
  sp_rates[11] += fwd_rates[329];
  //sp 15
  sp_rates[14] += fwd_rates[329];

  //rxn 330
  //sp 49
  sp_rates[48] -= (fwd_rates[330] - rev_rates[320]);
  //sp 26
  sp_rates[25] += (fwd_rates[330] - rev_rates[320]);
  //sp 12
  sp_rates[11] += (fwd_rates[330] - rev_rates[320]);

  //rxn 331
  //sp 49
  sp_rates[48] -= (fwd_rates[331] - rev_rates[321]);
  //sp 26
  sp_rates[25] += (fwd_rates[331] - rev_rates[321]);
  //sp 12
  sp_rates[11] += (fwd_rates[331] - rev_rates[321]);

  //rxn 332
  //sp 33
  sp_rates[32] -= (fwd_rates[332] - rev_rates[322]);
  //sp 49
  sp_rates[48] += (fwd_rates[332] - rev_rates[322]);
  //sp 7
  sp_rates[6] -= (fwd_rates[332] - rev_rates[322]);

  //rxn 333
  //sp 33
  sp_rates[32] -= (fwd_rates[333] - rev_rates[323]);
  //sp 4
  sp_rates[3] += (fwd_rates[333] - rev_rates[323]);
  //sp 7
  sp_rates[6] -= (fwd_rates[333] - rev_rates[323]);
  //sp 48
  sp_rates[47] += (fwd_rates[333] - rev_rates[323]);

  //rxn 334
  //sp 33
  sp_rates[32] -= (fwd_rates[334] - rev_rates[324]);
  //sp 14
  sp_rates[13] += (fwd_rates[334] - rev_rates[324]);
  //sp 7
  sp_rates[6] -= (fwd_rates[334] - rev_rates[324]);
  //sp 15
  sp_rates[14] += (fwd_rates[334] - rev_rates[324]);

  //rxn 335
  //sp 33
  sp_rates[32] -= (fwd_rates[335] - rev_rates[325]);
  //sp 3
  sp_rates[2] += (fwd_rates[335] - rev_rates[325]);
  //sp 7
  sp_rates[6] -= (fwd_rates[335] - rev_rates[325]);
  //sp 12
  sp_rates[11] += (fwd_rates[335] - rev_rates[325]);
  //sp 15
  sp_rates[14] += (fwd_rates[335] - rev_rates[325]);

  //rxn 336
  //sp 5
  sp_rates[4] -= fwd_rates[336];
  //sp 10
  sp_rates[9] += fwd_rates[336];
  //sp 12
  sp_rates[11] += fwd_rates[336];
  //sp 14
  sp_rates[13] += fwd_rates[336];
  //sp 56
  sp_rates[55] -= fwd_rates[336];

  //rxn 337
  //sp 12
  sp_rates[11] += (fwd_rates[337] - rev_rates[326]) * pres_mod[42];
  //sp 15
  sp_rates[14] += (fwd_rates[337] - rev_rates[326]) * pres_mod[42];
  //sp 56
  sp_rates[55] -= (fwd_rates[337] - rev_rates[326]) * pres_mod[42];

  //rxn 338
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[338] - rev_rates[327]) * pres_mod[43];
  //sp 6
  sp_rates[5] += (fwd_rates[338] - rev_rates[327]) * pres_mod[43];
  //sp 56
  sp_rates[55] -= (fwd_rates[338] - rev_rates[327]) * pres_mod[43];

  //rxn 339
  //sp 57
  sp_rates[56] = (fwd_rates[339] - rev_rates[328]);
  //sp 10
  sp_rates[9] += (fwd_rates[339] - rev_rates[328]);
  //sp 5
  sp_rates[4] -= (fwd_rates[339] - rev_rates[328]);
  //sp 56
  sp_rates[55] -= (fwd_rates[339] - rev_rates[328]);

  //rxn 340
  //sp 57
  sp_rates[56] += (fwd_rates[340] - rev_rates[329]);
  //sp 4
  sp_rates[3] -= (fwd_rates[340] - rev_rates[329]);
  //sp 5
  sp_rates[4] += (fwd_rates[340] - rev_rates[329]);
  //sp 56
  sp_rates[55] -= (fwd_rates[340] - rev_rates[329]);

  //rxn 341
  //sp 3
  sp_rates[2] -= (fwd_rates[341] - rev_rates[330]);
  //sp 14
  sp_rates[13] += (fwd_rates[341] - rev_rates[330]);
  //sp 15
  sp_rates[14] += (fwd_rates[341] - rev_rates[330]);
  //sp 56
  sp_rates[55] -= (fwd_rates[341] - rev_rates[330]);

  //rxn 342
  //sp 57
  sp_rates[56] += (fwd_rates[342] - rev_rates[331]);
  //sp 9
  sp_rates[8] += (fwd_rates[342] - rev_rates[331]);
  //sp 8
  sp_rates[7] -= (fwd_rates[342] - rev_rates[331]);
  //sp 56
  sp_rates[55] -= (fwd_rates[342] - rev_rates[331]);

  //rxn 343
  //sp 57
  sp_rates[56] += (fwd_rates[343] - rev_rates[332]);
  //sp 21
  sp_rates[20] -= (fwd_rates[343] - rev_rates[332]);
  //sp 22
  sp_rates[21] += (fwd_rates[343] - rev_rates[332]);
  //sp 56
  sp_rates[55] -= (fwd_rates[343] - rev_rates[332]);

  //rxn 344
  //sp 7
  sp_rates[6] -= (fwd_rates[344] - rev_rates[333]);
  //sp 8
  sp_rates[7] += (fwd_rates[344] - rev_rates[333]);
  //sp 12
  sp_rates[11] += (fwd_rates[344] - rev_rates[333]);
  //sp 14
  sp_rates[13] += (fwd_rates[344] - rev_rates[333]);
  //sp 56
  sp_rates[55] -= (fwd_rates[344] - rev_rates[333]);

  //rxn 345
  //sp 57
  sp_rates[56] -= (fwd_rates[345] - rev_rates[334]);
  //sp 12
  sp_rates[11] += (fwd_rates[345] - rev_rates[334]);
  //sp 14
  sp_rates[13] += (fwd_rates[345] - rev_rates[334]);

  //rxn 346
  //sp 57
  sp_rates[56] -= (fwd_rates[346] - rev_rates[335]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[346] - rev_rates[335]);
  //sp 7
  sp_rates[6] -= (fwd_rates[346] - rev_rates[335]);
  //sp 8
  sp_rates[7] += (fwd_rates[346] - rev_rates[335]);

  //rxn 347
  //sp 4
  sp_rates[3] -= (fwd_rates[347] - rev_rates[336]);
  //sp 5
  sp_rates[4] += (fwd_rates[347] - rev_rates[336]);
  //sp 12
  sp_rates[11] += (fwd_rates[347] - rev_rates[336]);
  //sp 14
  sp_rates[13] += (fwd_rates[347] - rev_rates[336]);
  //sp 56
  sp_rates[55] -= (fwd_rates[347] - rev_rates[336]);

  //rxn 348
  //sp 8
  sp_rates[7] -= (fwd_rates[348] - rev_rates[337]);
  //sp 9
  sp_rates[8] += (fwd_rates[348] - rev_rates[337]);
  //sp 12
  sp_rates[11] += (fwd_rates[348] - rev_rates[337]);
  //sp 14
  sp_rates[13] += (fwd_rates[348] - rev_rates[337]);
  //sp 56
  sp_rates[55] -= (fwd_rates[348] - rev_rates[337]);

  //rxn 349
  //sp 12
  sp_rates[11] += (fwd_rates[349] - rev_rates[338]);
  //sp 14
  sp_rates[13] += (fwd_rates[349] - rev_rates[338]);
  //sp 21
  sp_rates[20] -= (fwd_rates[349] - rev_rates[338]);
  //sp 22
  sp_rates[21] += (fwd_rates[349] - rev_rates[338]);
  //sp 56
  sp_rates[55] -= (fwd_rates[349] - rev_rates[338]);

  //rxn 350
  //sp 12
  sp_rates[11] += (fwd_rates[350] - rev_rates[339]) * pres_mod[44];
  //sp 21
  sp_rates[20] += (fwd_rates[350] - rev_rates[339]) * pres_mod[44];
  //sp 38
  sp_rates[37] -= (fwd_rates[350] - rev_rates[339]) * pres_mod[44];

  //rxn 351
  //sp 6
  sp_rates[5] += (fwd_rates[351] - rev_rates[340]);
  //sp 3
  sp_rates[2] -= (fwd_rates[351] - rev_rates[340]);
  //sp 28
  sp_rates[27] += (fwd_rates[351] - rev_rates[340]);
  //sp 38
  sp_rates[37] -= (fwd_rates[351] - rev_rates[340]);

  //rxn 352
  //sp 28
  sp_rates[27] += (fwd_rates[352] - rev_rates[341]);
  //sp 4
  sp_rates[3] -= (fwd_rates[352] - rev_rates[341]);
  //sp 5
  sp_rates[4] += (fwd_rates[352] - rev_rates[341]);
  //sp 38
  sp_rates[37] -= (fwd_rates[352] - rev_rates[341]);

  //rxn 353
  //sp 13
  sp_rates[12] += (fwd_rates[353] - rev_rates[342]);
  //sp 4
  sp_rates[3] -= (fwd_rates[353] - rev_rates[342]);
  //sp 21
  sp_rates[20] += (fwd_rates[353] - rev_rates[342]);
  //sp 38
  sp_rates[37] -= (fwd_rates[353] - rev_rates[342]);

  //rxn 354
  //sp 14
  sp_rates[13] += (fwd_rates[354] - rev_rates[343]);
  //sp 4
  sp_rates[3] -= (fwd_rates[354] - rev_rates[343]);
  //sp 38
  sp_rates[37] -= (fwd_rates[354] - rev_rates[343]);
  //sp 15
  sp_rates[14] += (fwd_rates[354] - rev_rates[343]);

  //rxn 355
  //sp 3
  sp_rates[2] -= (fwd_rates[355] - rev_rates[344]) * pres_mod[45];
  //sp 37
  sp_rates[36] += (fwd_rates[355] - rev_rates[344]) * pres_mod[45];
  //sp 38
  sp_rates[37] -= (fwd_rates[355] - rev_rates[344]) * pres_mod[45];

  //rxn 356
  //sp 3
  sp_rates[2] -= (fwd_rates[356] - rev_rates[345]);
  //sp 14
  sp_rates[13] += (fwd_rates[356] - rev_rates[345]);
  //sp 21
  sp_rates[20] += (fwd_rates[356] - rev_rates[345]);
  //sp 38
  sp_rates[37] -= (fwd_rates[356] - rev_rates[345]);

  //rxn 357
  //sp 10
  sp_rates[9] += (fwd_rates[357] - rev_rates[346]);
  //sp 28
  sp_rates[27] += (fwd_rates[357] - rev_rates[346]);
  //sp 5
  sp_rates[4] -= (fwd_rates[357] - rev_rates[346]);
  //sp 38
  sp_rates[37] -= (fwd_rates[357] - rev_rates[346]);

  //rxn 358
  //sp 12
  sp_rates[11] += (fwd_rates[358] - rev_rates[347]) * pres_mod[46];
  //sp 21
  sp_rates[20] += (fwd_rates[358] - rev_rates[347]) * pres_mod[46];
  //sp 38
  sp_rates[37] -= (fwd_rates[358] - rev_rates[347]) * pres_mod[46];

  //rxn 359
  //sp 5
  sp_rates[4] += (fwd_rates[359] - rev_rates[348]);
  //sp 38
  sp_rates[37] -= (fwd_rates[359] - rev_rates[348]);
  //sp 8
  sp_rates[7] -= (fwd_rates[359] - rev_rates[348]);
  //sp 13
  sp_rates[12] += (fwd_rates[359] - rev_rates[348]);
  //sp 21
  sp_rates[20] += (fwd_rates[359] - rev_rates[348]);

  //rxn 360
  //sp 3
  sp_rates[2] += (fwd_rates[360] - rev_rates[349]) * pres_mod[47];
  //sp 28
  sp_rates[27] += (fwd_rates[360] - rev_rates[349]) * pres_mod[47];
  //sp 38
  sp_rates[37] -= (fwd_rates[360] - rev_rates[349]) * pres_mod[47];

  //rxn 361
  //sp 29
  sp_rates[28] += (fwd_rates[361] - rev_rates[350]);
  //sp 12
  sp_rates[11] += (fwd_rates[361] - rev_rates[350]);
  //sp 21
  sp_rates[20] -= (fwd_rates[361] - rev_rates[350]);
  //sp 38
  sp_rates[37] -= (fwd_rates[361] - rev_rates[350]);

  //rxn 362
  //sp 22
  sp_rates[21] += (fwd_rates[362] - rev_rates[351]);
  //sp 28
  sp_rates[27] += (fwd_rates[362] - rev_rates[351]);
  //sp 21
  sp_rates[20] -= (fwd_rates[362] - rev_rates[351]);
  //sp 38
  sp_rates[37] -= (fwd_rates[362] - rev_rates[351]);

  //rxn 363
  //sp 38
  sp_rates[37] -= (fwd_rates[363] - rev_rates[352]);
  //sp 13
  sp_rates[12] += (fwd_rates[363] - rev_rates[352]);
  //sp 21
  sp_rates[20] += (fwd_rates[363] - rev_rates[352]);
  //sp 26
  sp_rates[25] += (fwd_rates[363] - rev_rates[352]);
  //sp 31
  sp_rates[30] -= (fwd_rates[363] - rev_rates[352]);

  //rxn 364
  //sp 3
  sp_rates[2] += (fwd_rates[364] - rev_rates[353]) * pres_mod[48];
  //sp 28
  sp_rates[27] += (fwd_rates[364] - rev_rates[353]) * pres_mod[48];
  //sp 48
  sp_rates[47] -= (fwd_rates[364] - rev_rates[353]) * pres_mod[48];

  //rxn 365
  //sp 12
  sp_rates[11] += (fwd_rates[365] - rev_rates[354]) * pres_mod[49];
  //sp 21
  sp_rates[20] += (fwd_rates[365] - rev_rates[354]) * pres_mod[49];
  //sp 48
  sp_rates[47] -= (fwd_rates[365] - rev_rates[354]) * pres_mod[49];

  //rxn 366
  //sp 3
  sp_rates[2] -= (fwd_rates[366] - rev_rates[355]) * pres_mod[50];
  //sp 37
  sp_rates[36] += (fwd_rates[366] - rev_rates[355]) * pres_mod[50];
  //sp 48
  sp_rates[47] -= (fwd_rates[366] - rev_rates[355]) * pres_mod[50];

  //rxn 367
  //sp 3
  sp_rates[2] -= (fwd_rates[367] - rev_rates[356]);
  //sp 21
  sp_rates[20] += (fwd_rates[367] - rev_rates[356]);
  //sp 14
  sp_rates[13] += (fwd_rates[367] - rev_rates[356]);
  //sp 48
  sp_rates[47] -= (fwd_rates[367] - rev_rates[356]);

  //rxn 368
  //sp 3
  sp_rates[2] -= (fwd_rates[368] - rev_rates[357]);
  //sp 28
  sp_rates[27] += (fwd_rates[368] - rev_rates[357]);
  //sp 6
  sp_rates[5] += (fwd_rates[368] - rev_rates[357]);
  //sp 48
  sp_rates[47] -= (fwd_rates[368] - rev_rates[357]);

  //rxn 369
  //sp 4
  sp_rates[3] -= (fwd_rates[369] - rev_rates[358]);
  //sp 14
  sp_rates[13] += (fwd_rates[369] - rev_rates[358]);
  //sp 15
  sp_rates[14] += (fwd_rates[369] - rev_rates[358]);
  //sp 48
  sp_rates[47] -= (fwd_rates[369] - rev_rates[358]);

  //rxn 370
  //sp 10
  sp_rates[9] += (fwd_rates[370] - rev_rates[359]);
  //sp 28
  sp_rates[27] += (fwd_rates[370] - rev_rates[359]);
  //sp 5
  sp_rates[4] -= (fwd_rates[370] - rev_rates[359]);
  //sp 48
  sp_rates[47] -= (fwd_rates[370] - rev_rates[359]);

  //rxn 371
  //sp 14
  sp_rates[13] += (fwd_rates[371] - rev_rates[360]);
  //sp 5
  sp_rates[4] -= (fwd_rates[371] - rev_rates[360]);
  //sp 24
  sp_rates[23] += (fwd_rates[371] - rev_rates[360]);
  //sp 48
  sp_rates[47] -= (fwd_rates[371] - rev_rates[360]);

  //rxn 372
  //sp 5
  sp_rates[4] += (fwd_rates[372] - rev_rates[361]);
  //sp 8
  sp_rates[7] -= (fwd_rates[372] - rev_rates[361]);
  //sp 14
  sp_rates[13] += (fwd_rates[372] - rev_rates[361]);
  //sp 15
  sp_rates[14] += (fwd_rates[372] - rev_rates[361]);
  //sp 48
  sp_rates[47] -= (fwd_rates[372] - rev_rates[361]);

  //rxn 373
  //sp 8
  sp_rates[7] -= (fwd_rates[373] - rev_rates[362]);
  //sp 7
  sp_rates[6] += (fwd_rates[373] - rev_rates[362]);
  //sp 37
  sp_rates[36] += (fwd_rates[373] - rev_rates[362]);
  //sp 48
  sp_rates[47] -= (fwd_rates[373] - rev_rates[362]);

  //rxn 374
  //sp 28
  sp_rates[27] += (fwd_rates[374] - rev_rates[363]);
  //sp 8
  sp_rates[7] += (fwd_rates[374] - rev_rates[363]);
  //sp 7
  sp_rates[6] -= (fwd_rates[374] - rev_rates[363]);
  //sp 48
  sp_rates[47] -= (fwd_rates[374] - rev_rates[363]);

  //rxn 375
  //sp 5
  sp_rates[4] += fwd_rates[375];
  //sp 7
  sp_rates[6] -= fwd_rates[375];
  //sp 12
  sp_rates[11] += fwd_rates[375];
  //sp 15
  sp_rates[14] += fwd_rates[375];
  //sp 48
  sp_rates[47] -= fwd_rates[375];

  //rxn 376
  //sp 34
  sp_rates[33] += (fwd_rates[376] - rev_rates[364]);
  //sp 3
  sp_rates[2] += (fwd_rates[376] - rev_rates[364]);
  //sp 12
  sp_rates[11] += (fwd_rates[376] - rev_rates[364]);
  //sp 48
  sp_rates[47] -= (fwd_rates[376] - rev_rates[364]);
  //sp 21
  sp_rates[20] -= (fwd_rates[376] - rev_rates[364]);

  //rxn 377
  //sp 18
  sp_rates[17] -= (fwd_rates[377] - rev_rates[365]);
  //sp 27
  sp_rates[26] += (fwd_rates[377] - rev_rates[365]);
  //sp 14
  sp_rates[13] += (fwd_rates[377] - rev_rates[365]);
  //sp 48
  sp_rates[47] -= (fwd_rates[377] - rev_rates[365]);

  //rxn 378
  //sp 33
  sp_rates[32] += (fwd_rates[378] - rev_rates[366]);
  //sp 14
  sp_rates[13] += (fwd_rates[378] - rev_rates[366]);
  //sp 48
  sp_rates[47] -= (fwd_rates[378] - rev_rates[366]);
  //sp 16
  sp_rates[15] -= (fwd_rates[378] - rev_rates[366]);

  //rxn 379
  //sp 58
  sp_rates[57] = -(fwd_rates[379] - rev_rates[367]) * pres_mod[51];
  //sp 21
  sp_rates[20] += (fwd_rates[379] - rev_rates[367]) * pres_mod[51];
  //sp 13
  sp_rates[12] += (fwd_rates[379] - rev_rates[367]) * pres_mod[51];

  //rxn 380
  //sp 59
  sp_rates[58] = -(fwd_rates[380] - rev_rates[368]);
  //sp 38
  sp_rates[37] += (fwd_rates[380] - rev_rates[368]);
  //sp 7
  sp_rates[6] += (fwd_rates[380] - rev_rates[368]);

  //rxn 381
  //sp 59
  sp_rates[58] -= (fwd_rates[381] - rev_rates[369]);
  //sp 60
  sp_rates[59] = (fwd_rates[381] - rev_rates[369]);
  //sp 7
  sp_rates[6] += (fwd_rates[381] - rev_rates[369]);
  //sp 8
  sp_rates[7] -= (fwd_rates[381] - rev_rates[369]);

  //rxn 382
  //sp 9
  sp_rates[8] -= (fwd_rates[382] - rev_rates[370]);
  //sp 59
  sp_rates[58] -= (fwd_rates[382] - rev_rates[370]);
  //sp 60
  sp_rates[59] += (fwd_rates[382] - rev_rates[370]);
  //sp 8
  sp_rates[7] += (fwd_rates[382] - rev_rates[370]);

  //rxn 383
  //sp 59
  sp_rates[58] -= (fwd_rates[383] - rev_rates[371]);
  //sp 60
  sp_rates[59] += (fwd_rates[383] - rev_rates[371]);
  //sp 14
  sp_rates[13] += (fwd_rates[383] - rev_rates[371]);
  //sp 15
  sp_rates[14] -= (fwd_rates[383] - rev_rates[371]);

  //rxn 384
  //sp 59
  sp_rates[58] -= (fwd_rates[384] - rev_rates[372]);
  //sp 60
  sp_rates[59] += (fwd_rates[384] - rev_rates[372]);
  //sp 21
  sp_rates[20] += (fwd_rates[384] - rev_rates[372]);
  //sp 22
  sp_rates[21] -= (fwd_rates[384] - rev_rates[372]);

  //rxn 385
  //sp 33
  sp_rates[32] += (fwd_rates[385] - rev_rates[373]);
  //sp 27
  sp_rates[26] -= (fwd_rates[385] - rev_rates[373]);
  //sp 60
  sp_rates[59] += (fwd_rates[385] - rev_rates[373]);
  //sp 59
  sp_rates[58] -= (fwd_rates[385] - rev_rates[373]);

  //rxn 386
  //sp 34
  sp_rates[33] += (fwd_rates[386] - rev_rates[374]);
  //sp 59
  sp_rates[58] -= (fwd_rates[386] - rev_rates[374]);
  //sp 60
  sp_rates[59] += (fwd_rates[386] - rev_rates[374]);
  //sp 29
  sp_rates[28] -= (fwd_rates[386] - rev_rates[374]);

  //rxn 387
  //sp 59
  sp_rates[58] -= (fwd_rates[387] - rev_rates[375]);
  //sp 60
  sp_rates[59] += (fwd_rates[387] - rev_rates[375]);
  //sp 37
  sp_rates[36] -= (fwd_rates[387] - rev_rates[375]);
  //sp 38
  sp_rates[37] += (fwd_rates[387] - rev_rates[375]);

  //rxn 388
  //sp 58
  sp_rates[57] += (fwd_rates[388] - rev_rates[376]);
  //sp 60
  sp_rates[59] -= (fwd_rates[388] - rev_rates[376]);
  //sp 5
  sp_rates[4] += (fwd_rates[388] - rev_rates[376]);

  //rxn 389
  //sp 27
  sp_rates[26] -= (fwd_rates[389] - rev_rates[377]) * pres_mod[52];
  //sp 44
  sp_rates[43] += (fwd_rates[389] - rev_rates[377]) * pres_mod[52];
  //sp 6
  sp_rates[5] += (fwd_rates[389] - rev_rates[377]) * pres_mod[52];

  //rxn 390
  //sp 33
  sp_rates[32] -= (fwd_rates[390] - rev_rates[378]) * pres_mod[53];
  //sp 3
  sp_rates[2] -= (fwd_rates[390] - rev_rates[378]) * pres_mod[53];
  //sp 27
  sp_rates[26] += (fwd_rates[390] - rev_rates[378]) * pres_mod[53];

  //rxn 391
  //sp 33
  sp_rates[32] += (fwd_rates[391] - rev_rates[379]);
  //sp 27
  sp_rates[26] -= (fwd_rates[391] - rev_rates[379]);
  //sp 3
  sp_rates[2] -= (fwd_rates[391] - rev_rates[379]);
  //sp 6
  sp_rates[5] += (fwd_rates[391] - rev_rates[379]);

  //rxn 392
  //sp 27
  sp_rates[26] -= (fwd_rates[392] - rev_rates[380]);
  //sp 4
  sp_rates[3] -= (fwd_rates[392] - rev_rates[380]);
  //sp 21
  sp_rates[20] += (fwd_rates[392] - rev_rates[380]);
  //sp 14
  sp_rates[13] += (fwd_rates[392] - rev_rates[380]);

  //rxn 393
  //sp 3
  sp_rates[2] += (fwd_rates[393] - rev_rates[381]);
  //sp 27
  sp_rates[26] -= (fwd_rates[393] - rev_rates[381]);
  //sp 4
  sp_rates[3] -= (fwd_rates[393] - rev_rates[381]);
  //sp 48
  sp_rates[47] += (fwd_rates[393] - rev_rates[381]);

  //rxn 394
  //sp 28
  sp_rates[27] += (fwd_rates[394] - rev_rates[382]);
  //sp 27
  sp_rates[26] -= (fwd_rates[394] - rev_rates[382]);
  //sp 4
  sp_rates[3] -= (fwd_rates[394] - rev_rates[382]);
  //sp 6
  sp_rates[5] += (fwd_rates[394] - rev_rates[382]);

  //rxn 395
  //sp 3
  sp_rates[2] += (fwd_rates[395] - rev_rates[383]);
  //sp 19
  sp_rates[18] -= (fwd_rates[395] - rev_rates[383]);
  //sp 27
  sp_rates[26] += (fwd_rates[395] - rev_rates[383]);
  //sp 21
  sp_rates[20] -= (fwd_rates[395] - rev_rates[383]);

  //rxn 396
  //sp 33
  sp_rates[32] += (fwd_rates[396] - rev_rates[384]);
  //sp 34
  sp_rates[33] += (fwd_rates[396] - rev_rates[384]);
  //sp 27
  sp_rates[26] -= 2.0 * (fwd_rates[396] - rev_rates[384]);

  //rxn 397
  //sp 65
  sp_rates[64] = (fwd_rates[397] - rev_rates[385]);
  //sp 27
  sp_rates[26] -= (fwd_rates[397] - rev_rates[385]);
  //sp 3
  sp_rates[2] += (fwd_rates[397] - rev_rates[385]);
  //sp 16
  sp_rates[15] -= (fwd_rates[397] - rev_rates[385]);

  //rxn 398
  //sp 66
  sp_rates[65] = (fwd_rates[398] - rev_rates[386]);
  //sp 27
  sp_rates[26] -= (fwd_rates[398] - rev_rates[386]);
  //sp 3
  sp_rates[2] += (fwd_rates[398] - rev_rates[386]);
  //sp 16
  sp_rates[15] -= (fwd_rates[398] - rev_rates[386]);

  //rxn 399
  //sp 33
  sp_rates[32] += (fwd_rates[399] - rev_rates[387]);
  //sp 27
  sp_rates[26] -= (fwd_rates[399] - rev_rates[387]);
  //sp 21
  sp_rates[20] -= (fwd_rates[399] - rev_rates[387]);
  //sp 22
  sp_rates[21] += (fwd_rates[399] - rev_rates[387]);

  //rxn 400
  //sp 33
  sp_rates[32] += (fwd_rates[400] - rev_rates[388]);
  //sp 9
  sp_rates[8] += (fwd_rates[400] - rev_rates[388]);
  //sp 27
  sp_rates[26] -= (fwd_rates[400] - rev_rates[388]);
  //sp 8
  sp_rates[7] -= (fwd_rates[400] - rev_rates[388]);

  //rxn 401
  //sp 18
  sp_rates[17] += (fwd_rates[401] - rev_rates[389]);
  //sp 19
  sp_rates[18] -= (fwd_rates[401] - rev_rates[389]);

  //rxn 402
  //sp 33
  sp_rates[32] += (fwd_rates[402] - rev_rates[390]);
  //sp 27
  sp_rates[26] -= (fwd_rates[402] - rev_rates[390]);
  //sp 7
  sp_rates[6] -= (fwd_rates[402] - rev_rates[390]);
  //sp 8
  sp_rates[7] += (fwd_rates[402] - rev_rates[390]);

  //rxn 403
  //sp 33
  sp_rates[32] += (fwd_rates[403] - rev_rates[391]);
  //sp 26
  sp_rates[25] -= (fwd_rates[403] - rev_rates[391]);
  //sp 27
  sp_rates[26] -= (fwd_rates[403] - rev_rates[391]);
  //sp 30
  sp_rates[29] += (fwd_rates[403] - rev_rates[391]);

  //rxn 404
  //sp 27
  sp_rates[26] -= (fwd_rates[404] - rev_rates[392]);
  //sp 37
  sp_rates[36] += (fwd_rates[404] - rev_rates[392]);
  //sp 5
  sp_rates[4] += (fwd_rates[404] - rev_rates[392]);
  //sp 8
  sp_rates[7] -= (fwd_rates[404] - rev_rates[392]);

  //rxn 405
  //sp 33
  sp_rates[32] += (fwd_rates[405] - rev_rates[393]);
  //sp 10
  sp_rates[9] += (fwd_rates[405] - rev_rates[393]);
  //sp 27
  sp_rates[26] -= (fwd_rates[405] - rev_rates[393]);
  //sp 5
  sp_rates[4] -= (fwd_rates[405] - rev_rates[393]);

  //rxn 406
  //sp 27
  sp_rates[26] -= (fwd_rates[406] - rev_rates[394]);
  //sp 5
  sp_rates[4] -= (fwd_rates[406] - rev_rates[394]);
  //sp 21
  sp_rates[20] += (fwd_rates[406] - rev_rates[394]);
  //sp 15
  sp_rates[14] += (fwd_rates[406] - rev_rates[394]);

  //rxn 407
  //sp 27
  sp_rates[26] -= (fwd_rates[407] - rev_rates[395]);
  //sp 3
  sp_rates[2] += (fwd_rates[407] - rev_rates[395]);
  //sp 5
  sp_rates[4] -= (fwd_rates[407] - rev_rates[395]);
  //sp 37
  sp_rates[36] += (fwd_rates[407] - rev_rates[395]);

  //rxn 408
  //sp 27
  sp_rates[26] -= (fwd_rates[408] - rev_rates[396]);
  //sp 3
  sp_rates[2] += (fwd_rates[408] - rev_rates[396]);
  //sp 5
  sp_rates[4] -= (fwd_rates[408] - rev_rates[396]);
  //sp 67
  sp_rates[66] = (fwd_rates[408] - rev_rates[396]);

  //rxn 409
  //sp 14
  sp_rates[13] += (fwd_rates[409] - rev_rates[397]) * pres_mod[54];
  //sp 37
  sp_rates[36] -= (fwd_rates[409] - rev_rates[397]) * pres_mod[54];
  //sp 21
  sp_rates[20] += (fwd_rates[409] - rev_rates[397]) * pres_mod[54];

  //rxn 410
  //sp 3
  sp_rates[2] -= (fwd_rates[410] - rev_rates[398]);
  //sp 6
  sp_rates[5] += (fwd_rates[410] - rev_rates[398]);
  //sp 37
  sp_rates[36] -= (fwd_rates[410] - rev_rates[398]);
  //sp 38
  sp_rates[37] += (fwd_rates[410] - rev_rates[398]);

  //rxn 411
  //sp 3
  sp_rates[2] -= (fwd_rates[411] - rev_rates[399]);
  //sp 37
  sp_rates[36] -= (fwd_rates[411] - rev_rates[399]);
  //sp 6
  sp_rates[5] += (fwd_rates[411] - rev_rates[399]);
  //sp 48
  sp_rates[47] += (fwd_rates[411] - rev_rates[399]);

  //rxn 412
  //sp 5
  sp_rates[4] += (fwd_rates[412] - rev_rates[400]);
  //sp 4
  sp_rates[3] -= (fwd_rates[412] - rev_rates[400]);
  //sp 37
  sp_rates[36] -= (fwd_rates[412] - rev_rates[400]);
  //sp 38
  sp_rates[37] += (fwd_rates[412] - rev_rates[400]);

  //rxn 413
  //sp 4
  sp_rates[3] -= (fwd_rates[413] - rev_rates[401]);
  //sp 37
  sp_rates[36] -= (fwd_rates[413] - rev_rates[401]);
  //sp 5
  sp_rates[4] += (fwd_rates[413] - rev_rates[401]);
  //sp 48
  sp_rates[47] += (fwd_rates[413] - rev_rates[401]);

  //rxn 414
  //sp 10
  sp_rates[9] += (fwd_rates[414] - rev_rates[402]);
  //sp 38
  sp_rates[37] += (fwd_rates[414] - rev_rates[402]);
  //sp 37
  sp_rates[36] -= (fwd_rates[414] - rev_rates[402]);
  //sp 5
  sp_rates[4] -= (fwd_rates[414] - rev_rates[402]);

  //rxn 415
  //sp 10
  sp_rates[9] += (fwd_rates[415] - rev_rates[403]);
  //sp 37
  sp_rates[36] -= (fwd_rates[415] - rev_rates[403]);
  //sp 5
  sp_rates[4] -= (fwd_rates[415] - rev_rates[403]);
  //sp 48
  sp_rates[47] += (fwd_rates[415] - rev_rates[403]);

  //rxn 416
  //sp 9
  sp_rates[8] += (fwd_rates[416] - rev_rates[404]);
  //sp 37
  sp_rates[36] -= (fwd_rates[416] - rev_rates[404]);
  //sp 38
  sp_rates[37] += (fwd_rates[416] - rev_rates[404]);
  //sp 8
  sp_rates[7] -= (fwd_rates[416] - rev_rates[404]);

  //rxn 417
  //sp 9
  sp_rates[8] += (fwd_rates[417] - rev_rates[405]);
  //sp 37
  sp_rates[36] -= (fwd_rates[417] - rev_rates[405]);
  //sp 48
  sp_rates[47] += (fwd_rates[417] - rev_rates[405]);
  //sp 8
  sp_rates[7] -= (fwd_rates[417] - rev_rates[405]);

  //rxn 418
  //sp 22
  sp_rates[21] += (fwd_rates[418] - rev_rates[406]);
  //sp 21
  sp_rates[20] -= (fwd_rates[418] - rev_rates[406]);
  //sp 38
  sp_rates[37] += (fwd_rates[418] - rev_rates[406]);
  //sp 37
  sp_rates[36] -= (fwd_rates[418] - rev_rates[406]);

  //rxn 419
  //sp 21
  sp_rates[20] -= (fwd_rates[419] - rev_rates[407]);
  //sp 22
  sp_rates[21] += (fwd_rates[419] - rev_rates[407]);
  //sp 37
  sp_rates[36] -= (fwd_rates[419] - rev_rates[407]);
  //sp 48
  sp_rates[47] += (fwd_rates[419] - rev_rates[407]);

  //rxn 420
  //sp 37
  sp_rates[36] -= (fwd_rates[420] - rev_rates[408]);
  //sp 38
  sp_rates[37] += (fwd_rates[420] - rev_rates[408]);
  //sp 7
  sp_rates[6] -= (fwd_rates[420] - rev_rates[408]);
  //sp 8
  sp_rates[7] += (fwd_rates[420] - rev_rates[408]);

  //rxn 421
  //sp 12
  sp_rates[11] += (fwd_rates[421] - rev_rates[409]) * pres_mod[55];
  //sp 37
  sp_rates[36] -= (fwd_rates[421] - rev_rates[409]) * pres_mod[55];
  //sp 22
  sp_rates[21] += (fwd_rates[421] - rev_rates[409]) * pres_mod[55];

  //rxn 422
  //sp 21
  sp_rates[20] += (fwd_rates[422] - rev_rates[410]);
  //sp 37
  sp_rates[36] -= (fwd_rates[422] - rev_rates[410]);
  //sp 5
  sp_rates[4] -= (fwd_rates[422] - rev_rates[410]);
  //sp 23
  sp_rates[22] += (fwd_rates[422] - rev_rates[410]);

  //rxn 423
  //sp 18
  sp_rates[17] -= (fwd_rates[423] - rev_rates[411]);
  //sp 38
  sp_rates[37] += (fwd_rates[423] - rev_rates[411]);
  //sp 37
  sp_rates[36] -= (fwd_rates[423] - rev_rates[411]);
  //sp 21
  sp_rates[20] += (fwd_rates[423] - rev_rates[411]);

  //rxn 424
  //sp 38
  sp_rates[37] += (fwd_rates[424] - rev_rates[412]);
  //sp 37
  sp_rates[36] -= (fwd_rates[424] - rev_rates[412]);
  //sp 14
  sp_rates[13] -= (fwd_rates[424] - rev_rates[412]);
  //sp 15
  sp_rates[14] += (fwd_rates[424] - rev_rates[412]);

  //rxn 425
  //sp 26
  sp_rates[25] -= (fwd_rates[425] - rev_rates[413]);
  //sp 30
  sp_rates[29] += (fwd_rates[425] - rev_rates[413]);
  //sp 37
  sp_rates[36] -= (fwd_rates[425] - rev_rates[413]);
  //sp 38
  sp_rates[37] += (fwd_rates[425] - rev_rates[413]);

  //rxn 426
  //sp 67
  sp_rates[66] -= (fwd_rates[426] - rev_rates[414]);
  //sp 3
  sp_rates[2] -= (fwd_rates[426] - rev_rates[414]);
  //sp 6
  sp_rates[5] += (fwd_rates[426] - rev_rates[414]);
  //sp 48
  sp_rates[47] += (fwd_rates[426] - rev_rates[414]);

  //rxn 427
  //sp 67
  sp_rates[66] -= (fwd_rates[427] - rev_rates[415]);
  //sp 4
  sp_rates[3] -= (fwd_rates[427] - rev_rates[415]);
  //sp 5
  sp_rates[4] += (fwd_rates[427] - rev_rates[415]);
  //sp 48
  sp_rates[47] += (fwd_rates[427] - rev_rates[415]);

  //rxn 428
  //sp 10
  sp_rates[9] += (fwd_rates[428] - rev_rates[416]);
  //sp 67
  sp_rates[66] -= (fwd_rates[428] - rev_rates[416]);
  //sp 5
  sp_rates[4] -= (fwd_rates[428] - rev_rates[416]);
  //sp 48
  sp_rates[47] += (fwd_rates[428] - rev_rates[416]);

  //rxn 429
  //sp 9
  sp_rates[8] += (fwd_rates[429] - rev_rates[417]);
  //sp 67
  sp_rates[66] -= (fwd_rates[429] - rev_rates[417]);
  //sp 48
  sp_rates[47] += (fwd_rates[429] - rev_rates[417]);
  //sp 8
  sp_rates[7] -= (fwd_rates[429] - rev_rates[417]);

  //rxn 430
  //sp 67
  sp_rates[66] -= (fwd_rates[430] - rev_rates[418]);
  //sp 21
  sp_rates[20] -= (fwd_rates[430] - rev_rates[418]);
  //sp 22
  sp_rates[21] += (fwd_rates[430] - rev_rates[418]);
  //sp 48
  sp_rates[47] += (fwd_rates[430] - rev_rates[418]);

  //rxn 431
  //sp 67
  sp_rates[66] -= (fwd_rates[431] - rev_rates[419]);
  //sp 37
  sp_rates[36] += (fwd_rates[431] - rev_rates[419]);

  //rxn 432
  //sp 67
  sp_rates[66] -= (fwd_rates[432] - rev_rates[420]);
  //sp 5
  sp_rates[4] -= (fwd_rates[432] - rev_rates[420]);
  //sp 6
  sp_rates[5] += (fwd_rates[432] - rev_rates[420]);
  //sp 13
  sp_rates[12] += (fwd_rates[432] - rev_rates[420]);
  //sp 21
  sp_rates[20] += (fwd_rates[432] - rev_rates[420]);

  //rxn 433
  //sp 67
  sp_rates[66] -= (fwd_rates[433] - rev_rates[421]);
  //sp 8
  sp_rates[7] += (fwd_rates[433] - rev_rates[421]);
  //sp 7
  sp_rates[6] -= (fwd_rates[433] - rev_rates[421]);
  //sp 48
  sp_rates[47] += (fwd_rates[433] - rev_rates[421]);

  //rxn 434
  //sp 67
  sp_rates[66] -= (fwd_rates[434] - rev_rates[422]);
  //sp 32
  sp_rates[31] += (fwd_rates[434] - rev_rates[422]);
  //sp 31
  sp_rates[30] -= (fwd_rates[434] - rev_rates[422]);
  //sp 48
  sp_rates[47] += (fwd_rates[434] - rev_rates[422]);

  //rxn 435
  //sp 67
  sp_rates[66] -= (fwd_rates[435] - rev_rates[423]);
  //sp 37
  sp_rates[36] += (fwd_rates[435] - rev_rates[423]);

  //rxn 436
  //sp 67
  sp_rates[66] -= (fwd_rates[436] - rev_rates[424]);
  //sp 37
  sp_rates[36] += (fwd_rates[436] - rev_rates[424]);

  //rxn 437
  //sp 52
  sp_rates[51] += (fwd_rates[437] - rev_rates[425]);
  //sp 67
  sp_rates[66] -= (fwd_rates[437] - rev_rates[425]);
  //sp 3
  sp_rates[2] -= (fwd_rates[437] - rev_rates[425]);
  //sp 6
  sp_rates[5] += (fwd_rates[437] - rev_rates[425]);

  //rxn 438
  //sp 34
  sp_rates[33] += (fwd_rates[438] - rev_rates[426]) * pres_mod[56];
  //sp 27
  sp_rates[26] -= (fwd_rates[438] - rev_rates[426]) * pres_mod[56];
  //sp 3
  sp_rates[2] -= (fwd_rates[438] - rev_rates[426]) * pres_mod[56];

  //rxn 439
  //sp 34
  sp_rates[33] -= (fwd_rates[439] - rev_rates[427]);
  //sp 3
  sp_rates[2] -= (fwd_rates[439] - rev_rates[427]);
  //sp 27
  sp_rates[26] += (fwd_rates[439] - rev_rates[427]);
  //sp 6
  sp_rates[5] += (fwd_rates[439] - rev_rates[427]);

  //rxn 440
  //sp 34
  sp_rates[33] -= (fwd_rates[440] - rev_rates[428]);
  //sp 27
  sp_rates[26] += (fwd_rates[440] - rev_rates[428]);
  //sp 5
  sp_rates[4] -= (fwd_rates[440] - rev_rates[428]);
  //sp 10
  sp_rates[9] += (fwd_rates[440] - rev_rates[428]);

  //rxn 441
  //sp 34
  sp_rates[33] -= (fwd_rates[441] - rev_rates[429]);
  //sp 5
  sp_rates[4] -= (fwd_rates[441] - rev_rates[429]);
  //sp 21
  sp_rates[20] += (fwd_rates[441] - rev_rates[429]);
  //sp 24
  sp_rates[23] += (fwd_rates[441] - rev_rates[429]);

  //rxn 442
  //sp 34
  sp_rates[33] -= (fwd_rates[442] - rev_rates[430]);
  //sp 3
  sp_rates[2] += (fwd_rates[442] - rev_rates[430]);
  //sp 4
  sp_rates[3] -= (fwd_rates[442] - rev_rates[430]);
  //sp 37
  sp_rates[36] += (fwd_rates[442] - rev_rates[430]);

  //rxn 443
  //sp 34
  sp_rates[33] -= (fwd_rates[443] - rev_rates[431]);
  //sp 4
  sp_rates[3] -= (fwd_rates[443] - rev_rates[431]);
  //sp 21
  sp_rates[20] += (fwd_rates[443] - rev_rates[431]);
  //sp 15
  sp_rates[14] += (fwd_rates[443] - rev_rates[431]);

  //rxn 444
  //sp 34
  sp_rates[33] -= (fwd_rates[444] - rev_rates[432]);
  //sp 27
  sp_rates[26] += (fwd_rates[444] - rev_rates[432]);
  //sp 4
  sp_rates[3] -= (fwd_rates[444] - rev_rates[432]);
  //sp 5
  sp_rates[4] += (fwd_rates[444] - rev_rates[432]);

  //rxn 445
  //sp 34
  sp_rates[33] -= (fwd_rates[445] - rev_rates[433]);
  //sp 36
  sp_rates[35] += (fwd_rates[445] - rev_rates[433]);
  //sp 5
  sp_rates[4] += (fwd_rates[445] - rev_rates[433]);
  //sp 8
  sp_rates[7] -= (fwd_rates[445] - rev_rates[433]);

  //rxn 446
  //sp 34
  sp_rates[33] += (fwd_rates[446] - rev_rates[434]);
  //sp 3
  sp_rates[2] += (fwd_rates[446] - rev_rates[434]);
  //sp 21
  sp_rates[20] -= 2.0 * (fwd_rates[446] - rev_rates[434]);

  //rxn 447
  //sp 34
  sp_rates[33] -= (fwd_rates[447] - rev_rates[435]);
  //sp 27
  sp_rates[26] += (fwd_rates[447] - rev_rates[435]);
  //sp 21
  sp_rates[20] -= (fwd_rates[447] - rev_rates[435]);
  //sp 22
  sp_rates[21] += (fwd_rates[447] - rev_rates[435]);

  //rxn 448
  //sp 34
  sp_rates[33] -= 2.0 * (fwd_rates[448] - rev_rates[436]);
  //sp 27
  sp_rates[26] += (fwd_rates[448] - rev_rates[436]);
  //sp 29
  sp_rates[28] += (fwd_rates[448] - rev_rates[436]);

  //rxn 449
  //sp 34
  sp_rates[33] -= (fwd_rates[449] - rev_rates[437]);
  //sp 37
  sp_rates[36] += (fwd_rates[449] - rev_rates[437]);
  //sp 5
  sp_rates[4] += (fwd_rates[449] - rev_rates[437]);
  //sp 7
  sp_rates[6] -= (fwd_rates[449] - rev_rates[437]);

  //rxn 450
  //sp 34
  sp_rates[33] -= (fwd_rates[450] - rev_rates[438]) * pres_mod[57];
  //sp 27
  sp_rates[26] += (fwd_rates[450] - rev_rates[438]) * pres_mod[57];
  //sp 7
  sp_rates[6] -= (fwd_rates[450] - rev_rates[438]) * pres_mod[57];
  //sp 8
  sp_rates[7] += (fwd_rates[450] - rev_rates[438]) * pres_mod[57];

  //rxn 451
  //sp 34
  sp_rates[33] -= (fwd_rates[451] - rev_rates[439]);
  //sp 3
  sp_rates[2] += (fwd_rates[451] - rev_rates[439]);
  //sp 4
  sp_rates[3] -= (fwd_rates[451] - rev_rates[439]);
  //sp 12
  sp_rates[11] += (fwd_rates[451] - rev_rates[439]);
  //sp 22
  sp_rates[21] += (fwd_rates[451] - rev_rates[439]);

  //rxn 452
  //sp 34
  sp_rates[33] -= (fwd_rates[452] - rev_rates[440]);
  //sp 4
  sp_rates[3] -= (fwd_rates[452] - rev_rates[440]);
  //sp 6
  sp_rates[5] += (fwd_rates[452] - rev_rates[440]);
  //sp 12
  sp_rates[11] += (fwd_rates[452] - rev_rates[440]);
  //sp 21
  sp_rates[20] += (fwd_rates[452] - rev_rates[440]);

  //rxn 453
  //sp 3
  sp_rates[2] += (fwd_rates[453] - rev_rates[441]) * pres_mod[58];
  //sp 36
  sp_rates[35] -= (fwd_rates[453] - rev_rates[441]) * pres_mod[58];
  //sp 37
  sp_rates[36] += (fwd_rates[453] - rev_rates[441]) * pres_mod[58];

  //rxn 454
  //sp 36
  sp_rates[35] -= (fwd_rates[454] - rev_rates[442]) * pres_mod[59];
  //sp 21
  sp_rates[20] += (fwd_rates[454] - rev_rates[442]) * pres_mod[59];
  //sp 15
  sp_rates[14] += (fwd_rates[454] - rev_rates[442]) * pres_mod[59];

  //rxn 455
  //sp 36
  sp_rates[35] -= (fwd_rates[455] - rev_rates[443]);
  //sp 37
  sp_rates[36] += (fwd_rates[455] - rev_rates[443]);
  //sp 7
  sp_rates[6] -= (fwd_rates[455] - rev_rates[443]);
  //sp 8
  sp_rates[7] += (fwd_rates[455] - rev_rates[443]);

  //rxn 456
  //sp 10
  sp_rates[9] += (fwd_rates[456] - rev_rates[444]);
  //sp 36
  sp_rates[35] -= (fwd_rates[456] - rev_rates[444]);
  //sp 5
  sp_rates[4] -= (fwd_rates[456] - rev_rates[444]);
  //sp 37
  sp_rates[36] += (fwd_rates[456] - rev_rates[444]);

  //rxn 457
  //sp 3
  sp_rates[2] -= (fwd_rates[457] - rev_rates[445]);
  //sp 36
  sp_rates[35] -= (fwd_rates[457] - rev_rates[445]);
  //sp 37
  sp_rates[36] += (fwd_rates[457] - rev_rates[445]);
  //sp 6
  sp_rates[5] += (fwd_rates[457] - rev_rates[445]);

  //rxn 458
  //sp 3
  sp_rates[2] -= (fwd_rates[458] - rev_rates[446]);
  //sp 36
  sp_rates[35] -= (fwd_rates[458] - rev_rates[446]);
  //sp 21
  sp_rates[20] += (fwd_rates[458] - rev_rates[446]);
  //sp 24
  sp_rates[23] += (fwd_rates[458] - rev_rates[446]);

  //rxn 459
  //sp 27
  sp_rates[26] += (fwd_rates[459] - rev_rates[447]);
  //sp 10
  sp_rates[9] += (fwd_rates[459] - rev_rates[447]);
  //sp 3
  sp_rates[2] -= (fwd_rates[459] - rev_rates[447]);
  //sp 36
  sp_rates[35] -= (fwd_rates[459] - rev_rates[447]);

  //rxn 460
  //sp 4
  sp_rates[3] -= (fwd_rates[460] - rev_rates[448]);
  //sp 36
  sp_rates[35] -= (fwd_rates[460] - rev_rates[448]);
  //sp 37
  sp_rates[36] += (fwd_rates[460] - rev_rates[448]);
  //sp 5
  sp_rates[4] += (fwd_rates[460] - rev_rates[448]);

  //rxn 461
  //sp 12
  sp_rates[11] -= (fwd_rates[461] - rev_rates[449]);
  //sp 34
  sp_rates[33] += (fwd_rates[461] - rev_rates[449]);
  //sp 36
  sp_rates[35] -= (fwd_rates[461] - rev_rates[449]);
  //sp 13
  sp_rates[12] += (fwd_rates[461] - rev_rates[449]);

  //rxn 462
  //sp 3
  sp_rates[2] += (fwd_rates[462] - rev_rates[450]);
  //sp 68
  sp_rates[67] = -(fwd_rates[462] - rev_rates[450]);
  //sp 37
  sp_rates[36] += (fwd_rates[462] - rev_rates[450]);

  //rxn 463
  //sp 68
  sp_rates[67] -= (fwd_rates[463] - rev_rates[451]);
  //sp 21
  sp_rates[20] += (fwd_rates[463] - rev_rates[451]);
  //sp 15
  sp_rates[14] += (fwd_rates[463] - rev_rates[451]);

  //rxn 464
  //sp 3
  sp_rates[2] -= (fwd_rates[464] - rev_rates[452]);
  //sp 68
  sp_rates[67] -= (fwd_rates[464] - rev_rates[452]);
  //sp 37
  sp_rates[36] += (fwd_rates[464] - rev_rates[452]);
  //sp 6
  sp_rates[5] += (fwd_rates[464] - rev_rates[452]);

  //rxn 465
  //sp 3
  sp_rates[2] -= (fwd_rates[465] - rev_rates[453]);
  //sp 68
  sp_rates[67] -= (fwd_rates[465] - rev_rates[453]);
  //sp 21
  sp_rates[20] += (fwd_rates[465] - rev_rates[453]);
  //sp 24
  sp_rates[23] += (fwd_rates[465] - rev_rates[453]);

  //rxn 466
  //sp 27
  sp_rates[26] += (fwd_rates[466] - rev_rates[454]);
  //sp 10
  sp_rates[9] += (fwd_rates[466] - rev_rates[454]);
  //sp 3
  sp_rates[2] -= (fwd_rates[466] - rev_rates[454]);
  //sp 68
  sp_rates[67] -= (fwd_rates[466] - rev_rates[454]);

  //rxn 467
  //sp 10
  sp_rates[9] += (fwd_rates[467] - rev_rates[455]);
  //sp 68
  sp_rates[67] -= (fwd_rates[467] - rev_rates[455]);
  //sp 5
  sp_rates[4] -= (fwd_rates[467] - rev_rates[455]);
  //sp 37
  sp_rates[36] += (fwd_rates[467] - rev_rates[455]);

  //rxn 468
  //sp 4
  sp_rates[3] -= (fwd_rates[468] - rev_rates[456]);
  //sp 68
  sp_rates[67] -= (fwd_rates[468] - rev_rates[456]);
  //sp 37
  sp_rates[36] += (fwd_rates[468] - rev_rates[456]);
  //sp 5
  sp_rates[4] += (fwd_rates[468] - rev_rates[456]);

  //rxn 469
  //sp 68
  sp_rates[67] -= (fwd_rates[469] - rev_rates[457]);
  //sp 37
  sp_rates[36] += (fwd_rates[469] - rev_rates[457]);
  //sp 7
  sp_rates[6] -= (fwd_rates[469] - rev_rates[457]);
  //sp 8
  sp_rates[7] += (fwd_rates[469] - rev_rates[457]);

  //rxn 470
  //sp 68
  sp_rates[67] -= (fwd_rates[470] - rev_rates[458]);
  //sp 37
  sp_rates[36] += (fwd_rates[470] - rev_rates[458]);
  //sp 5
  sp_rates[4] += 2.0 * (fwd_rates[470] - rev_rates[458]);
  //sp 8
  sp_rates[7] -= (fwd_rates[470] - rev_rates[458]);

  //rxn 471
  //sp 3
  sp_rates[2] += (fwd_rates[471] - rev_rates[459]);
  //sp 67
  sp_rates[66] += (fwd_rates[471] - rev_rates[459]);
  //sp 68
  sp_rates[67] -= (fwd_rates[471] - rev_rates[459]);

  //rxn 472
  //sp 6
  sp_rates[5] += (fwd_rates[472] - rev_rates[460]);
  //sp 3
  sp_rates[2] -= (fwd_rates[472] - rev_rates[460]);
  //sp 68
  sp_rates[67] -= (fwd_rates[472] - rev_rates[460]);
  //sp 67
  sp_rates[66] += (fwd_rates[472] - rev_rates[460]);

  //rxn 473
  //sp 67
  sp_rates[66] += (fwd_rates[473] - rev_rates[461]);
  //sp 68
  sp_rates[67] -= (fwd_rates[473] - rev_rates[461]);
  //sp 7
  sp_rates[6] -= (fwd_rates[473] - rev_rates[461]);
  //sp 8
  sp_rates[7] += (fwd_rates[473] - rev_rates[461]);

  //rxn 474
  //sp 27
  sp_rates[26] -= (fwd_rates[474] - rev_rates[462]);
  //sp 5
  sp_rates[4] -= (fwd_rates[474] - rev_rates[462]);
  //sp 69
  sp_rates[68] = (fwd_rates[474] - rev_rates[462]);

  //rxn 475
  //sp 68
  sp_rates[67] -= (fwd_rates[475] - rev_rates[463]);
  //sp 69
  sp_rates[68] += (fwd_rates[475] - rev_rates[463]);

  //rxn 476
  //sp 3
  sp_rates[2] -= (fwd_rates[476] - rev_rates[464]);
  //sp 69
  sp_rates[68] -= (fwd_rates[476] - rev_rates[464]);
  //sp 6
  sp_rates[5] += (fwd_rates[476] - rev_rates[464]);
  //sp 37
  sp_rates[36] += (fwd_rates[476] - rev_rates[464]);

  //rxn 477
  //sp 67
  sp_rates[66] += (fwd_rates[477] - rev_rates[465]);
  //sp 69
  sp_rates[68] -= (fwd_rates[477] - rev_rates[465]);
  //sp 7
  sp_rates[6] -= (fwd_rates[477] - rev_rates[465]);
  //sp 8
  sp_rates[7] += (fwd_rates[477] - rev_rates[465]);

  //rxn 478
  //sp 67
  sp_rates[66] += (fwd_rates[478] - rev_rates[466]);
  //sp 3
  sp_rates[2] += (fwd_rates[478] - rev_rates[466]);
  //sp 69
  sp_rates[68] -= (fwd_rates[478] - rev_rates[466]);

  //rxn 479
  //sp 3
  sp_rates[2] -= (fwd_rates[479] - rev_rates[467]);
  //sp 6
  sp_rates[5] += (fwd_rates[479] - rev_rates[467]);
  //sp 69
  sp_rates[68] -= (fwd_rates[479] - rev_rates[467]);
  //sp 67
  sp_rates[66] += (fwd_rates[479] - rev_rates[467]);

  //rxn 480
  //sp 67
  sp_rates[66] += (fwd_rates[480] - rev_rates[468]);
  //sp 4
  sp_rates[3] -= (fwd_rates[480] - rev_rates[468]);
  //sp 69
  sp_rates[68] -= (fwd_rates[480] - rev_rates[468]);
  //sp 5
  sp_rates[4] += (fwd_rates[480] - rev_rates[468]);

  //rxn 481
  //sp 10
  sp_rates[9] += (fwd_rates[481] - rev_rates[469]);
  //sp 67
  sp_rates[66] += (fwd_rates[481] - rev_rates[469]);
  //sp 69
  sp_rates[68] -= (fwd_rates[481] - rev_rates[469]);
  //sp 5
  sp_rates[4] -= (fwd_rates[481] - rev_rates[469]);

  //rxn 482
  //sp 69
  sp_rates[68] -= (fwd_rates[482] - rev_rates[470]);
  //sp 8
  sp_rates[7] += (fwd_rates[482] - rev_rates[470]);
  //sp 7
  sp_rates[6] -= (fwd_rates[482] - rev_rates[470]);
  //sp 37
  sp_rates[36] += (fwd_rates[482] - rev_rates[470]);

  //rxn 483
  //sp 69
  sp_rates[68] -= (fwd_rates[483] - rev_rates[471]);
  //sp 8
  sp_rates[7] += (fwd_rates[483] - rev_rates[471]);
  //sp 7
  sp_rates[6] -= (fwd_rates[483] - rev_rates[471]);
  //sp 37
  sp_rates[36] += (fwd_rates[483] - rev_rates[471]);

  //rxn 484
  //sp 69
  sp_rates[68] -= fwd_rates[484];
  //sp 5
  sp_rates[4] += fwd_rates[484];
  //sp 7
  sp_rates[6] -= fwd_rates[484];
  //sp 15
  sp_rates[14] += 2.0 * fwd_rates[484];

  //rxn 485
  //sp 10
  sp_rates[9] += fwd_rates[485];
  //sp 69
  sp_rates[68] -= fwd_rates[485];
  //sp 15
  sp_rates[14] += 2.0 * fwd_rates[485];
  //sp 8
  sp_rates[7] -= fwd_rates[485];

  //rxn 486
  //sp 34
  sp_rates[33] -= (fwd_rates[486] - rev_rates[472]) * pres_mod[60];
  //sp 70
  sp_rates[69] = (fwd_rates[486] - rev_rates[472]) * pres_mod[60];
  //sp 7
  sp_rates[6] -= (fwd_rates[486] - rev_rates[472]) * pres_mod[60];

  //rxn 487
  //sp 5
  sp_rates[4] += (fwd_rates[487] - rev_rates[473]);
  //sp 37
  sp_rates[36] += (fwd_rates[487] - rev_rates[473]);
  //sp 70
  sp_rates[69] -= (fwd_rates[487] - rev_rates[473]);

  //rxn 488
  //sp 27
  sp_rates[26] += (fwd_rates[488] - rev_rates[474]) * pres_mod[61];
  //sp 70
  sp_rates[69] -= (fwd_rates[488] - rev_rates[474]) * pres_mod[61];
  //sp 8
  sp_rates[7] += (fwd_rates[488] - rev_rates[474]) * pres_mod[61];

  //rxn 489
  //sp 5
  sp_rates[4] += (fwd_rates[489] - rev_rates[475]);
  //sp 70
  sp_rates[69] -= (fwd_rates[489] - rev_rates[475]);
  //sp 71
  sp_rates[70] = (fwd_rates[489] - rev_rates[475]);

  //rxn 490
  //sp 7
  sp_rates[6] += (fwd_rates[490] - rev_rates[476]);
  //sp 70
  sp_rates[69] -= (fwd_rates[490] - rev_rates[476]);
  //sp 72
  sp_rates[71] = (fwd_rates[490] - rev_rates[476]);
  //sp 8
  sp_rates[7] -= (fwd_rates[490] - rev_rates[476]);

  //rxn 491
  //sp 22
  sp_rates[21] -= (fwd_rates[491] - rev_rates[477]);
  //sp 21
  sp_rates[20] += (fwd_rates[491] - rev_rates[477]);
  //sp 70
  sp_rates[69] -= (fwd_rates[491] - rev_rates[477]);
  //sp 72
  sp_rates[71] += (fwd_rates[491] - rev_rates[477]);

  //rxn 492
  //sp 14
  sp_rates[13] += (fwd_rates[492] - rev_rates[478]);
  //sp 70
  sp_rates[69] -= (fwd_rates[492] - rev_rates[478]);
  //sp 15
  sp_rates[14] -= (fwd_rates[492] - rev_rates[478]);
  //sp 72
  sp_rates[71] += (fwd_rates[492] - rev_rates[478]);

  //rxn 493
  //sp 30
  sp_rates[29] -= (fwd_rates[493] - rev_rates[479]);
  //sp 24
  sp_rates[23] += (fwd_rates[493] - rev_rates[479]);
  //sp 70
  sp_rates[69] -= (fwd_rates[493] - rev_rates[479]);
  //sp 72
  sp_rates[71] += (fwd_rates[493] - rev_rates[479]);

  //rxn 494
  //sp 33
  sp_rates[32] += (fwd_rates[494] - rev_rates[480]);
  //sp 27
  sp_rates[26] -= (fwd_rates[494] - rev_rates[480]);
  //sp 70
  sp_rates[69] -= (fwd_rates[494] - rev_rates[480]);
  //sp 72
  sp_rates[71] += (fwd_rates[494] - rev_rates[480]);

  //rxn 495
  //sp 34
  sp_rates[33] += (fwd_rates[495] - rev_rates[481]);
  //sp 29
  sp_rates[28] -= (fwd_rates[495] - rev_rates[481]);
  //sp 70
  sp_rates[69] -= (fwd_rates[495] - rev_rates[481]);
  //sp 72
  sp_rates[71] += (fwd_rates[495] - rev_rates[481]);

  //rxn 496
  //sp 61
  sp_rates[60] = -(fwd_rates[496] - rev_rates[482]);
  //sp 70
  sp_rates[69] -= (fwd_rates[496] - rev_rates[482]);
  //sp 47
  sp_rates[46] += (fwd_rates[496] - rev_rates[482]);
  //sp 72
  sp_rates[71] += (fwd_rates[496] - rev_rates[482]);

  //rxn 497
  //sp 73
  sp_rates[72] = -(fwd_rates[497] - rev_rates[483]);
  //sp 74
  sp_rates[73] = (fwd_rates[497] - rev_rates[483]);
  //sp 70
  sp_rates[69] -= (fwd_rates[497] - rev_rates[483]);
  //sp 72
  sp_rates[71] += (fwd_rates[497] - rev_rates[483]);

  //rxn 498
  //sp 75
  sp_rates[74] = (fwd_rates[498] - rev_rates[484]);
  //sp 70
  sp_rates[69] -= (fwd_rates[498] - rev_rates[484]);

  //rxn 499
  //sp 34
  sp_rates[33] -= (fwd_rates[499] - rev_rates[485]);
  //sp 75
  sp_rates[74] += (fwd_rates[499] - rev_rates[485]);
  //sp 7
  sp_rates[6] -= (fwd_rates[499] - rev_rates[485]);

  //rxn 500
  //sp 27
  sp_rates[26] -= (fwd_rates[500] - rev_rates[486]);
  //sp 75
  sp_rates[74] += (fwd_rates[500] - rev_rates[486]);
  //sp 8
  sp_rates[7] -= (fwd_rates[500] - rev_rates[486]);

  //rxn 501
  //sp 75
  sp_rates[74] -= (fwd_rates[501] - rev_rates[487]);
  //sp 37
  sp_rates[36] += (fwd_rates[501] - rev_rates[487]);
  //sp 5
  sp_rates[4] += (fwd_rates[501] - rev_rates[487]);

  //rxn 502
  //sp 75
  sp_rates[74] -= (fwd_rates[502] - rev_rates[488]);
  //sp 5
  sp_rates[4] += (fwd_rates[502] - rev_rates[488]);
  //sp 71
  sp_rates[70] += (fwd_rates[502] - rev_rates[488]);

  //rxn 503
  //sp 34
  sp_rates[33] -= (fwd_rates[503] - rev_rates[489]);
  //sp 5
  sp_rates[4] += (fwd_rates[503] - rev_rates[489]);
  //sp 7
  sp_rates[6] -= (fwd_rates[503] - rev_rates[489]);
  //sp 71
  sp_rates[70] += (fwd_rates[503] - rev_rates[489]);

  //rxn 504
  //sp 21
  sp_rates[20] += (fwd_rates[504] - rev_rates[490]);
  //sp 14
  sp_rates[13] += (fwd_rates[504] - rev_rates[490]);
  //sp 71
  sp_rates[70] -= (fwd_rates[504] - rev_rates[490]);

  //rxn 505
  //sp 37
  sp_rates[36] += (fwd_rates[505] - rev_rates[491]);
  //sp 71
  sp_rates[70] -= (fwd_rates[505] - rev_rates[491]);

  //rxn 506
  //sp 10
  sp_rates[9] += (fwd_rates[506] - rev_rates[492]);
  //sp 76
  sp_rates[75] = (fwd_rates[506] - rev_rates[492]);
  //sp 5
  sp_rates[4] -= (fwd_rates[506] - rev_rates[492]);
  //sp 71
  sp_rates[70] -= (fwd_rates[506] - rev_rates[492]);

  //rxn 507
  //sp 3
  sp_rates[2] -= (fwd_rates[507] - rev_rates[493]);
  //sp 76
  sp_rates[75] += (fwd_rates[507] - rev_rates[493]);
  //sp 6
  sp_rates[5] += (fwd_rates[507] - rev_rates[493]);
  //sp 71
  sp_rates[70] -= (fwd_rates[507] - rev_rates[493]);

  //rxn 508
  //sp 76
  sp_rates[75] += (fwd_rates[508] - rev_rates[494]);
  //sp 32
  sp_rates[31] += (fwd_rates[508] - rev_rates[494]);
  //sp 71
  sp_rates[70] -= (fwd_rates[508] - rev_rates[494]);
  //sp 31
  sp_rates[30] -= (fwd_rates[508] - rev_rates[494]);

  //rxn 509
  //sp 76
  sp_rates[75] += (fwd_rates[509] - rev_rates[495]);
  //sp 70
  sp_rates[69] -= (fwd_rates[509] - rev_rates[495]);
  //sp 71
  sp_rates[70] -= (fwd_rates[509] - rev_rates[495]);
  //sp 72
  sp_rates[71] += (fwd_rates[509] - rev_rates[495]);

  //rxn 510
  //sp 76
  sp_rates[75] += (fwd_rates[510] - rev_rates[496]);
  //sp 21
  sp_rates[20] -= (fwd_rates[510] - rev_rates[496]);
  //sp 22
  sp_rates[21] += (fwd_rates[510] - rev_rates[496]);
  //sp 71
  sp_rates[70] -= (fwd_rates[510] - rev_rates[496]);

  //rxn 511
  //sp 26
  sp_rates[25] -= (fwd_rates[511] - rev_rates[497]);
  //sp 76
  sp_rates[75] += (fwd_rates[511] - rev_rates[497]);
  //sp 30
  sp_rates[29] += (fwd_rates[511] - rev_rates[497]);
  //sp 71
  sp_rates[70] -= (fwd_rates[511] - rev_rates[497]);

  //rxn 512
  //sp 27
  sp_rates[26] -= (fwd_rates[512] - rev_rates[498]);
  //sp 5
  sp_rates[4] += (fwd_rates[512] - rev_rates[498]);
  //sp 71
  sp_rates[70] += (fwd_rates[512] - rev_rates[498]);
  //sp 8
  sp_rates[7] -= (fwd_rates[512] - rev_rates[498]);

  //rxn 513
  //sp 9
  sp_rates[8] += (fwd_rates[513] - rev_rates[499]);
  //sp 76
  sp_rates[75] += (fwd_rates[513] - rev_rates[499]);
  //sp 71
  sp_rates[70] -= (fwd_rates[513] - rev_rates[499]);
  //sp 8
  sp_rates[7] -= (fwd_rates[513] - rev_rates[499]);

  //rxn 514
  //sp 71
  sp_rates[70] += (fwd_rates[514] - rev_rates[500]);
  //sp 27
  sp_rates[26] -= (fwd_rates[514] - rev_rates[500]);
  //sp 31
  sp_rates[30] -= (fwd_rates[514] - rev_rates[500]);
  //sp 26
  sp_rates[25] += (fwd_rates[514] - rev_rates[500]);

  //rxn 515
  //sp 27
  sp_rates[26] -= (fwd_rates[515] - rev_rates[501]);
  //sp 36
  sp_rates[35] += (fwd_rates[515] - rev_rates[501]);
  //sp 70
  sp_rates[69] -= (fwd_rates[515] - rev_rates[501]);
  //sp 71
  sp_rates[70] += (fwd_rates[515] - rev_rates[501]);

  //rxn 516
  //sp 76
  sp_rates[75] -= (fwd_rates[516] - rev_rates[502]);
  //sp 38
  sp_rates[37] += (fwd_rates[516] - rev_rates[502]);

  //rxn 517
  //sp 76
  sp_rates[75] -= (fwd_rates[517] - rev_rates[503]);
  //sp 48
  sp_rates[47] += (fwd_rates[517] - rev_rates[503]);

  //rxn 518
  //sp 36
  sp_rates[35] += (fwd_rates[518] - rev_rates[504]);
  //sp 5
  sp_rates[4] += (fwd_rates[518] - rev_rates[504]);
  //sp 72
  sp_rates[71] -= (fwd_rates[518] - rev_rates[504]);

  //rxn 519
  //sp 34
  sp_rates[33] -= (fwd_rates[519] - rev_rates[505]) * pres_mod[62];
  //sp 3
  sp_rates[2] -= (fwd_rates[519] - rev_rates[505]) * pres_mod[62];
  //sp 29
  sp_rates[28] += (fwd_rates[519] - rev_rates[505]) * pres_mod[62];

  //rxn 520
  //sp 34
  sp_rates[33] += (fwd_rates[520] - rev_rates[506]);
  //sp 3
  sp_rates[2] -= (fwd_rates[520] - rev_rates[506]);
  //sp 29
  sp_rates[28] -= (fwd_rates[520] - rev_rates[506]);
  //sp 6
  sp_rates[5] += (fwd_rates[520] - rev_rates[506]);

  //rxn 521
  //sp 34
  sp_rates[33] += (fwd_rates[521] - rev_rates[507]);
  //sp 4
  sp_rates[3] -= (fwd_rates[521] - rev_rates[507]);
  //sp 29
  sp_rates[28] -= (fwd_rates[521] - rev_rates[507]);
  //sp 5
  sp_rates[4] += (fwd_rates[521] - rev_rates[507]);

  //rxn 522
  //sp 34
  sp_rates[33] += (fwd_rates[522] - rev_rates[508]);
  //sp 29
  sp_rates[28] -= (fwd_rates[522] - rev_rates[508]);
  //sp 5
  sp_rates[4] -= (fwd_rates[522] - rev_rates[508]);
  //sp 10
  sp_rates[9] += (fwd_rates[522] - rev_rates[508]);

  //rxn 523
  //sp 34
  sp_rates[33] += (fwd_rates[523] - rev_rates[509]);
  //sp 29
  sp_rates[28] -= (fwd_rates[523] - rev_rates[509]);
  //sp 7
  sp_rates[6] -= (fwd_rates[523] - rev_rates[509]);
  //sp 8
  sp_rates[7] += (fwd_rates[523] - rev_rates[509]);

  //rxn 524
  //sp 9
  sp_rates[8] += (fwd_rates[524] - rev_rates[510]);
  //sp 34
  sp_rates[33] += (fwd_rates[524] - rev_rates[510]);
  //sp 29
  sp_rates[28] -= (fwd_rates[524] - rev_rates[510]);
  //sp 8
  sp_rates[7] -= (fwd_rates[524] - rev_rates[510]);

  //rxn 525
  //sp 34
  sp_rates[33] += (fwd_rates[525] - rev_rates[511]);
  //sp 22
  sp_rates[21] += (fwd_rates[525] - rev_rates[511]);
  //sp 29
  sp_rates[28] -= (fwd_rates[525] - rev_rates[511]);
  //sp 21
  sp_rates[20] -= (fwd_rates[525] - rev_rates[511]);

  //rxn 526
  //sp 34
  sp_rates[33] += (fwd_rates[526] - rev_rates[512]);
  //sp 22
  sp_rates[21] += (fwd_rates[526] - rev_rates[512]);
  //sp 29
  sp_rates[28] -= (fwd_rates[526] - rev_rates[512]);
  //sp 21
  sp_rates[20] -= (fwd_rates[526] - rev_rates[512]);

  //rxn 527
  //sp 26
  sp_rates[25] -= (fwd_rates[527] - rev_rates[513]);
  //sp 29
  sp_rates[28] -= (fwd_rates[527] - rev_rates[513]);
  //sp 30
  sp_rates[29] += (fwd_rates[527] - rev_rates[513]);
  //sp 34
  sp_rates[33] += (fwd_rates[527] - rev_rates[513]);

  //rxn 528
  //sp 34
  sp_rates[33] += (fwd_rates[528] - rev_rates[514]);
  //sp 29
  sp_rates[28] -= (fwd_rates[528] - rev_rates[514]);
  //sp 18
  sp_rates[17] += (fwd_rates[528] - rev_rates[514]);
  //sp 16
  sp_rates[15] -= (fwd_rates[528] - rev_rates[514]);

  //rxn 529
  //sp 34
  sp_rates[33] += (fwd_rates[529] - rev_rates[515]);
  //sp 19
  sp_rates[18] -= (fwd_rates[529] - rev_rates[515]);
  //sp 29
  sp_rates[28] -= (fwd_rates[529] - rev_rates[515]);
  //sp 21
  sp_rates[20] += (fwd_rates[529] - rev_rates[515]);

  //rxn 530
  //sp 34
  sp_rates[33] += (fwd_rates[530] - rev_rates[516]);
  //sp 29
  sp_rates[28] -= (fwd_rates[530] - rev_rates[516]);
  //sp 14
  sp_rates[13] -= (fwd_rates[530] - rev_rates[516]);
  //sp 15
  sp_rates[14] += (fwd_rates[530] - rev_rates[516]);

  //rxn 531
  //sp 35
  sp_rates[34] = -(fwd_rates[531] - rev_rates[517]) * pres_mod[63];
  //sp 21
  sp_rates[20] += (fwd_rates[531] - rev_rates[517]) * pres_mod[63];
  //sp 24
  sp_rates[23] += (fwd_rates[531] - rev_rates[517]) * pres_mod[63];

  //rxn 532
  //sp 10
  sp_rates[9] += (fwd_rates[532] - rev_rates[518]) * pres_mod[64];
  //sp 35
  sp_rates[34] -= (fwd_rates[532] - rev_rates[518]) * pres_mod[64];
  //sp 27
  sp_rates[26] += (fwd_rates[532] - rev_rates[518]) * pres_mod[64];

  //rxn 533
  //sp 35
  sp_rates[34] -= (fwd_rates[533] - rev_rates[519]) * pres_mod[65];
  //sp 37
  sp_rates[36] += (fwd_rates[533] - rev_rates[519]) * pres_mod[65];
  //sp 6
  sp_rates[5] += (fwd_rates[533] - rev_rates[519]) * pres_mod[65];

  //rxn 534
  //sp 34
  sp_rates[33] += (fwd_rates[534] - rev_rates[520]) * pres_mod[66];
  //sp 35
  sp_rates[34] -= (fwd_rates[534] - rev_rates[520]) * pres_mod[66];
  //sp 5
  sp_rates[4] += (fwd_rates[534] - rev_rates[520]) * pres_mod[66];

  //rxn 535
  //sp 10
  sp_rates[9] += (fwd_rates[535] - rev_rates[521]);
  //sp 35
  sp_rates[34] -= (fwd_rates[535] - rev_rates[521]);
  //sp 5
  sp_rates[4] -= (fwd_rates[535] - rev_rates[521]);
  //sp 69
  sp_rates[68] += (fwd_rates[535] - rev_rates[521]);

  //rxn 536
  //sp 10
  sp_rates[9] += (fwd_rates[536] - rev_rates[522]);
  //sp 35
  sp_rates[34] -= (fwd_rates[536] - rev_rates[522]);
  //sp 68
  sp_rates[67] += (fwd_rates[536] - rev_rates[522]);
  //sp 5
  sp_rates[4] -= (fwd_rates[536] - rev_rates[522]);

  //rxn 537
  //sp 10
  sp_rates[9] += (fwd_rates[537] - rev_rates[523]);
  //sp 35
  sp_rates[34] -= (fwd_rates[537] - rev_rates[523]);
  //sp 36
  sp_rates[35] += (fwd_rates[537] - rev_rates[523]);
  //sp 5
  sp_rates[4] -= (fwd_rates[537] - rev_rates[523]);

  //rxn 538
  //sp 35
  sp_rates[34] -= (fwd_rates[538] - rev_rates[524]);
  //sp 3
  sp_rates[2] -= (fwd_rates[538] - rev_rates[524]);
  //sp 69
  sp_rates[68] += (fwd_rates[538] - rev_rates[524]);
  //sp 6
  sp_rates[5] += (fwd_rates[538] - rev_rates[524]);

  //rxn 539
  //sp 68
  sp_rates[67] += (fwd_rates[539] - rev_rates[525]);
  //sp 35
  sp_rates[34] -= (fwd_rates[539] - rev_rates[525]);
  //sp 3
  sp_rates[2] -= (fwd_rates[539] - rev_rates[525]);
  //sp 6
  sp_rates[5] += (fwd_rates[539] - rev_rates[525]);

  //rxn 540
  //sp 36
  sp_rates[35] += (fwd_rates[540] - rev_rates[526]);
  //sp 35
  sp_rates[34] -= (fwd_rates[540] - rev_rates[526]);
  //sp 3
  sp_rates[2] -= (fwd_rates[540] - rev_rates[526]);
  //sp 6
  sp_rates[5] += (fwd_rates[540] - rev_rates[526]);

  //rxn 541
  //sp 35
  sp_rates[34] -= (fwd_rates[541] - rev_rates[527]);
  //sp 4
  sp_rates[3] -= (fwd_rates[541] - rev_rates[527]);
  //sp 69
  sp_rates[68] += (fwd_rates[541] - rev_rates[527]);
  //sp 5
  sp_rates[4] += (fwd_rates[541] - rev_rates[527]);

  //rxn 542
  //sp 68
  sp_rates[67] += (fwd_rates[542] - rev_rates[528]);
  //sp 35
  sp_rates[34] -= (fwd_rates[542] - rev_rates[528]);
  //sp 4
  sp_rates[3] -= (fwd_rates[542] - rev_rates[528]);
  //sp 5
  sp_rates[4] += (fwd_rates[542] - rev_rates[528]);

  //rxn 543
  //sp 36
  sp_rates[35] += (fwd_rates[543] - rev_rates[529]);
  //sp 35
  sp_rates[34] -= (fwd_rates[543] - rev_rates[529]);
  //sp 4
  sp_rates[3] -= (fwd_rates[543] - rev_rates[529]);
  //sp 5
  sp_rates[4] += (fwd_rates[543] - rev_rates[529]);

  //rxn 544
  //sp 35
  sp_rates[34] -= (fwd_rates[544] - rev_rates[530]);
  //sp 21
  sp_rates[20] -= (fwd_rates[544] - rev_rates[530]);
  //sp 22
  sp_rates[21] += (fwd_rates[544] - rev_rates[530]);
  //sp 69
  sp_rates[68] += (fwd_rates[544] - rev_rates[530]);

  //rxn 545
  //sp 35
  sp_rates[34] -= (fwd_rates[545] - rev_rates[531]);
  //sp 68
  sp_rates[67] += (fwd_rates[545] - rev_rates[531]);
  //sp 21
  sp_rates[20] -= (fwd_rates[545] - rev_rates[531]);
  //sp 22
  sp_rates[21] += (fwd_rates[545] - rev_rates[531]);

  //rxn 546
  //sp 35
  sp_rates[34] -= (fwd_rates[546] - rev_rates[532]);
  //sp 36
  sp_rates[35] += (fwd_rates[546] - rev_rates[532]);
  //sp 21
  sp_rates[20] -= (fwd_rates[546] - rev_rates[532]);
  //sp 22
  sp_rates[21] += (fwd_rates[546] - rev_rates[532]);

  //rxn 547
  //sp 9
  sp_rates[8] += (fwd_rates[547] - rev_rates[533]);
  //sp 35
  sp_rates[34] -= (fwd_rates[547] - rev_rates[533]);
  //sp 69
  sp_rates[68] += (fwd_rates[547] - rev_rates[533]);
  //sp 8
  sp_rates[7] -= (fwd_rates[547] - rev_rates[533]);

  //rxn 548
  //sp 9
  sp_rates[8] += (fwd_rates[548] - rev_rates[534]);
  //sp 35
  sp_rates[34] -= (fwd_rates[548] - rev_rates[534]);
  //sp 68
  sp_rates[67] += (fwd_rates[548] - rev_rates[534]);
  //sp 8
  sp_rates[7] -= (fwd_rates[548] - rev_rates[534]);

  //rxn 549
  //sp 9
  sp_rates[8] += (fwd_rates[549] - rev_rates[535]);
  //sp 35
  sp_rates[34] -= (fwd_rates[549] - rev_rates[535]);
  //sp 36
  sp_rates[35] += (fwd_rates[549] - rev_rates[535]);
  //sp 8
  sp_rates[7] -= (fwd_rates[549] - rev_rates[535]);

  //rxn 550
  //sp 35
  sp_rates[34] -= (fwd_rates[550] - rev_rates[536]);
  //sp 69
  sp_rates[68] += (fwd_rates[550] - rev_rates[536]);
  //sp 7
  sp_rates[6] -= (fwd_rates[550] - rev_rates[536]);
  //sp 8
  sp_rates[7] += (fwd_rates[550] - rev_rates[536]);

  //rxn 551
  //sp 35
  sp_rates[34] -= (fwd_rates[551] - rev_rates[537]);
  //sp 68
  sp_rates[67] += (fwd_rates[551] - rev_rates[537]);
  //sp 7
  sp_rates[6] -= (fwd_rates[551] - rev_rates[537]);
  //sp 8
  sp_rates[7] += (fwd_rates[551] - rev_rates[537]);

  //rxn 552
  //sp 35
  sp_rates[34] -= (fwd_rates[552] - rev_rates[538]);
  //sp 36
  sp_rates[35] += (fwd_rates[552] - rev_rates[538]);
  //sp 7
  sp_rates[6] -= (fwd_rates[552] - rev_rates[538]);
  //sp 8
  sp_rates[7] += (fwd_rates[552] - rev_rates[538]);

  //rxn 553
  //sp 34
  sp_rates[33] -= (fwd_rates[553] - rev_rates[539]);
  //sp 35
  sp_rates[34] -= (fwd_rates[553] - rev_rates[539]);
  //sp 29
  sp_rates[28] += (fwd_rates[553] - rev_rates[539]);
  //sp 69
  sp_rates[68] += (fwd_rates[553] - rev_rates[539]);

  //rxn 554
  //sp 34
  sp_rates[33] -= (fwd_rates[554] - rev_rates[540]);
  //sp 35
  sp_rates[34] -= (fwd_rates[554] - rev_rates[540]);
  //sp 68
  sp_rates[67] += (fwd_rates[554] - rev_rates[540]);
  //sp 29
  sp_rates[28] += (fwd_rates[554] - rev_rates[540]);

  //rxn 555
  //sp 34
  sp_rates[33] -= (fwd_rates[555] - rev_rates[541]);
  //sp 35
  sp_rates[34] -= (fwd_rates[555] - rev_rates[541]);
  //sp 36
  sp_rates[35] += (fwd_rates[555] - rev_rates[541]);
  //sp 29
  sp_rates[28] += (fwd_rates[555] - rev_rates[541]);

  //rxn 556
  //sp 26
  sp_rates[25] += (fwd_rates[556] - rev_rates[542]) * pres_mod[67];
  //sp 77
  sp_rates[76] = -(fwd_rates[556] - rev_rates[542]) * pres_mod[67];
  //sp 21
  sp_rates[20] += (fwd_rates[556] - rev_rates[542]) * pres_mod[67];

  //rxn 557
  //sp 10
  sp_rates[9] += (fwd_rates[557] - rev_rates[543]);
  //sp 78
  sp_rates[77] = (fwd_rates[557] - rev_rates[543]);
  //sp 77
  sp_rates[76] -= (fwd_rates[557] - rev_rates[543]);
  //sp 5
  sp_rates[4] -= (fwd_rates[557] - rev_rates[543]);

  //rxn 558
  //sp 3
  sp_rates[2] -= (fwd_rates[558] - rev_rates[544]);
  //sp 6
  sp_rates[5] += (fwd_rates[558] - rev_rates[544]);
  //sp 77
  sp_rates[76] -= (fwd_rates[558] - rev_rates[544]);
  //sp 78
  sp_rates[77] += (fwd_rates[558] - rev_rates[544]);

  //rxn 559
  //sp 22
  sp_rates[21] += (fwd_rates[559] - rev_rates[545]);
  //sp 21
  sp_rates[20] -= (fwd_rates[559] - rev_rates[545]);
  //sp 78
  sp_rates[77] += (fwd_rates[559] - rev_rates[545]);
  //sp 77
  sp_rates[76] -= (fwd_rates[559] - rev_rates[545]);

  //rxn 560
  //sp 5
  sp_rates[4] += (fwd_rates[560] - rev_rates[546]);
  //sp 4
  sp_rates[3] -= (fwd_rates[560] - rev_rates[546]);
  //sp 77
  sp_rates[76] -= (fwd_rates[560] - rev_rates[546]);
  //sp 78
  sp_rates[77] += (fwd_rates[560] - rev_rates[546]);

  //rxn 561
  //sp 9
  sp_rates[8] += (fwd_rates[561] - rev_rates[547]);
  //sp 77
  sp_rates[76] -= (fwd_rates[561] - rev_rates[547]);
  //sp 78
  sp_rates[77] += (fwd_rates[561] - rev_rates[547]);
  //sp 8
  sp_rates[7] -= (fwd_rates[561] - rev_rates[547]);

  //rxn 562
  //sp 77
  sp_rates[76] -= (fwd_rates[562] - rev_rates[548]);
  //sp 78
  sp_rates[77] += (fwd_rates[562] - rev_rates[548]);
  //sp 7
  sp_rates[6] -= (fwd_rates[562] - rev_rates[548]);
  //sp 8
  sp_rates[7] += (fwd_rates[562] - rev_rates[548]);

  //rxn 563
  //sp 26
  sp_rates[25] -= (fwd_rates[563] - rev_rates[549]);
  //sp 30
  sp_rates[29] += (fwd_rates[563] - rev_rates[549]);
  //sp 77
  sp_rates[76] -= (fwd_rates[563] - rev_rates[549]);
  //sp 78
  sp_rates[77] += (fwd_rates[563] - rev_rates[549]);

  //rxn 564
  //sp 77
  sp_rates[76] -= (fwd_rates[564] - rev_rates[550]);
  //sp 78
  sp_rates[77] += (fwd_rates[564] - rev_rates[550]);
  //sp 31
  sp_rates[30] -= (fwd_rates[564] - rev_rates[550]);
  //sp 32
  sp_rates[31] += (fwd_rates[564] - rev_rates[550]);

  //rxn 565
  //sp 77
  sp_rates[76] -= (fwd_rates[565] - rev_rates[551]);
  //sp 78
  sp_rates[77] += (fwd_rates[565] - rev_rates[551]);
  //sp 79
  sp_rates[78] = -(fwd_rates[565] - rev_rates[551]);
  //sp 80
  sp_rates[79] = (fwd_rates[565] - rev_rates[551]);

  //rxn 566
  //sp 25
  sp_rates[24] -= (fwd_rates[566] - rev_rates[552]);
  //sp 77
  sp_rates[76] -= (fwd_rates[566] - rev_rates[552]);
  //sp 78
  sp_rates[77] += (fwd_rates[566] - rev_rates[552]);
  //sp 23
  sp_rates[22] += (fwd_rates[566] - rev_rates[552]);

  //rxn 567
  //sp 21
  sp_rates[20] += (fwd_rates[567] - rev_rates[553]);
  //sp 78
  sp_rates[77] -= (fwd_rates[567] - rev_rates[553]);
  //sp 15
  sp_rates[14] += (fwd_rates[567] - rev_rates[553]);

  //rxn 568
  //sp 26
  sp_rates[25] -= (fwd_rates[568] - rev_rates[554]);
  //sp 77
  sp_rates[76] += (fwd_rates[568] - rev_rates[554]);
  //sp 78
  sp_rates[77] -= (fwd_rates[568] - rev_rates[554]);
  //sp 15
  sp_rates[14] += (fwd_rates[568] - rev_rates[554]);

  //rxn 569
  //sp 14
  sp_rates[13] += (fwd_rates[569] - rev_rates[555]);
  //sp 77
  sp_rates[76] += (fwd_rates[569] - rev_rates[555]);
  //sp 78
  sp_rates[77] -= (fwd_rates[569] - rev_rates[555]);
  //sp 15
  sp_rates[14] -= (fwd_rates[569] - rev_rates[555]);

  //rxn 570
  //sp 81
  sp_rates[80] = (fwd_rates[570] - rev_rates[556]);
  //sp 5
  sp_rates[4] += (fwd_rates[570] - rev_rates[556]);
  //sp 78
  sp_rates[77] -= (fwd_rates[570] - rev_rates[556]);
  //sp 8
  sp_rates[7] -= (fwd_rates[570] - rev_rates[556]);

  //rxn 571
  //sp 81
  sp_rates[80] -= (fwd_rates[571] - rev_rates[557]);
  //sp 82
  sp_rates[81] = (fwd_rates[571] - rev_rates[557]);
  //sp 3
  sp_rates[2] += (fwd_rates[571] - rev_rates[557]);

  //rxn 572
  //sp 25
  sp_rates[24] -= (fwd_rates[572] - rev_rates[558]) * pres_mod[68];
  //sp 3
  sp_rates[2] += (fwd_rates[572] - rev_rates[558]) * pres_mod[68];
  //sp 13
  sp_rates[12] += (fwd_rates[572] - rev_rates[558]) * pres_mod[68];

  //rxn 573
  //sp 25
  sp_rates[24] -= (fwd_rates[573] - rev_rates[559]);
  //sp 14
  sp_rates[13] += (fwd_rates[573] - rev_rates[559]);
  //sp 15
  sp_rates[14] -= (fwd_rates[573] - rev_rates[559]);
  //sp 23
  sp_rates[22] += (fwd_rates[573] - rev_rates[559]);

  //rxn 574
  //sp 82
  sp_rates[81] -= (fwd_rates[574] - rev_rates[560]) * pres_mod[69];
  //sp 12
  sp_rates[11] += (fwd_rates[574] - rev_rates[560]) * pres_mod[69];
  //sp 30
  sp_rates[29] += (fwd_rates[574] - rev_rates[560]) * pres_mod[69];

  //rxn 575
  //sp 82
  sp_rates[81] -= (fwd_rates[575] - rev_rates[561]) * pres_mod[70];
  //sp 13
  sp_rates[12] += (fwd_rates[575] - rev_rates[561]) * pres_mod[70];
  //sp 22
  sp_rates[21] += (fwd_rates[575] - rev_rates[561]) * pres_mod[70];

  //rxn 576
  //sp 82
  sp_rates[81] -= (fwd_rates[576] - rev_rates[562]) * pres_mod[71];
  //sp 15
  sp_rates[14] += 2.0 * (fwd_rates[576] - rev_rates[562]) * pres_mod[71];

  //rxn 577
  //sp 25
  sp_rates[24] += (fwd_rates[577] - rev_rates[563]) * pres_mod[72];
  //sp 82
  sp_rates[81] -= (fwd_rates[577] - rev_rates[563]) * pres_mod[72];
  //sp 21
  sp_rates[20] += (fwd_rates[577] - rev_rates[563]) * pres_mod[72];

  //rxn 578
  //sp 82
  sp_rates[81] -= (fwd_rates[578] - rev_rates[564]) * pres_mod[73];
  //sp 14
  sp_rates[13] += (fwd_rates[578] - rev_rates[564]) * pres_mod[73];
  //sp 26
  sp_rates[25] += (fwd_rates[578] - rev_rates[564]) * pres_mod[73];

  //rxn 579
  //sp 82
  sp_rates[81] += (fwd_rates[579] - rev_rates[565]);
  //sp 83
  sp_rates[82] = -(fwd_rates[579] - rev_rates[565]);
  //sp 3
  sp_rates[2] -= (fwd_rates[579] - rev_rates[565]);

  //rxn 580
  //sp 82
  sp_rates[81] += (fwd_rates[580] - rev_rates[566]);
  //sp 3
  sp_rates[2] -= (fwd_rates[580] - rev_rates[566]);
  //sp 84
  sp_rates[83] = -(fwd_rates[580] - rev_rates[566]);

  //rxn 581
  //sp 82
  sp_rates[81] -= (fwd_rates[581] - rev_rates[567]);
  //sp 84
  sp_rates[83] += (fwd_rates[581] - rev_rates[567]);
  //sp 7
  sp_rates[6] -= (fwd_rates[581] - rev_rates[567]);
  //sp 8
  sp_rates[7] += (fwd_rates[581] - rev_rates[567]);

  //rxn 582
  //sp 82
  sp_rates[81] -= (fwd_rates[582] - rev_rates[568]);
  //sp 83
  sp_rates[82] += (fwd_rates[582] - rev_rates[568]);
  //sp 7
  sp_rates[6] -= (fwd_rates[582] - rev_rates[568]);
  //sp 8
  sp_rates[7] += (fwd_rates[582] - rev_rates[568]);

  //rxn 583
  //sp 82
  sp_rates[81] -= (fwd_rates[583] - rev_rates[569]);
  //sp 84
  sp_rates[83] += (fwd_rates[583] - rev_rates[569]);
  //sp 5
  sp_rates[4] -= (fwd_rates[583] - rev_rates[569]);
  //sp 10
  sp_rates[9] += (fwd_rates[583] - rev_rates[569]);

  //rxn 584
  //sp 82
  sp_rates[81] -= (fwd_rates[584] - rev_rates[570]);
  //sp 83
  sp_rates[82] += (fwd_rates[584] - rev_rates[570]);
  //sp 5
  sp_rates[4] -= (fwd_rates[584] - rev_rates[570]);
  //sp 10
  sp_rates[9] += (fwd_rates[584] - rev_rates[570]);

  //rxn 585
  //sp 9
  sp_rates[8] += (fwd_rates[585] - rev_rates[571]);
  //sp 82
  sp_rates[81] -= (fwd_rates[585] - rev_rates[571]);
  //sp 84
  sp_rates[83] += (fwd_rates[585] - rev_rates[571]);
  //sp 8
  sp_rates[7] -= (fwd_rates[585] - rev_rates[571]);

  //rxn 586
  //sp 9
  sp_rates[8] += (fwd_rates[586] - rev_rates[572]);
  //sp 82
  sp_rates[81] -= (fwd_rates[586] - rev_rates[572]);
  //sp 83
  sp_rates[82] += (fwd_rates[586] - rev_rates[572]);
  //sp 8
  sp_rates[7] -= (fwd_rates[586] - rev_rates[572]);

  //rxn 587
  //sp 82
  sp_rates[81] -= (fwd_rates[587] - rev_rates[573]);
  //sp 84
  sp_rates[83] += (fwd_rates[587] - rev_rates[573]);
  //sp 4
  sp_rates[3] -= (fwd_rates[587] - rev_rates[573]);
  //sp 5
  sp_rates[4] += (fwd_rates[587] - rev_rates[573]);

  //rxn 588
  //sp 82
  sp_rates[81] -= (fwd_rates[588] - rev_rates[574]);
  //sp 83
  sp_rates[82] += (fwd_rates[588] - rev_rates[574]);
  //sp 4
  sp_rates[3] -= (fwd_rates[588] - rev_rates[574]);
  //sp 5
  sp_rates[4] += (fwd_rates[588] - rev_rates[574]);

  //rxn 589
  //sp 82
  sp_rates[81] -= (fwd_rates[589] - rev_rates[575]);
  //sp 3
  sp_rates[2] -= (fwd_rates[589] - rev_rates[575]);
  //sp 84
  sp_rates[83] += (fwd_rates[589] - rev_rates[575]);
  //sp 6
  sp_rates[5] += (fwd_rates[589] - rev_rates[575]);

  //rxn 590
  //sp 82
  sp_rates[81] -= (fwd_rates[590] - rev_rates[576]);
  //sp 3
  sp_rates[2] -= (fwd_rates[590] - rev_rates[576]);
  //sp 6
  sp_rates[5] += (fwd_rates[590] - rev_rates[576]);
  //sp 83
  sp_rates[82] += (fwd_rates[590] - rev_rates[576]);

  //rxn 591
  //sp 82
  sp_rates[81] -= (fwd_rates[591] - rev_rates[577]);
  //sp 84
  sp_rates[83] += (fwd_rates[591] - rev_rates[577]);
  //sp 21
  sp_rates[20] -= (fwd_rates[591] - rev_rates[577]);
  //sp 22
  sp_rates[21] += (fwd_rates[591] - rev_rates[577]);

  //rxn 592
  //sp 82
  sp_rates[81] -= (fwd_rates[592] - rev_rates[578]);
  //sp 83
  sp_rates[82] += (fwd_rates[592] - rev_rates[578]);
  //sp 21
  sp_rates[20] -= (fwd_rates[592] - rev_rates[578]);
  //sp 22
  sp_rates[21] += (fwd_rates[592] - rev_rates[578]);

  //rxn 593
  //sp 82
  sp_rates[81] -= (fwd_rates[593] - rev_rates[579]);
  //sp 26
  sp_rates[25] -= (fwd_rates[593] - rev_rates[579]);
  //sp 84
  sp_rates[83] += (fwd_rates[593] - rev_rates[579]);
  //sp 30
  sp_rates[29] += (fwd_rates[593] - rev_rates[579]);

  //rxn 594
  //sp 82
  sp_rates[81] -= (fwd_rates[594] - rev_rates[580]);
  //sp 26
  sp_rates[25] -= (fwd_rates[594] - rev_rates[580]);
  //sp 83
  sp_rates[82] += (fwd_rates[594] - rev_rates[580]);
  //sp 30
  sp_rates[29] += (fwd_rates[594] - rev_rates[580]);

  //rxn 595
  //sp 82
  sp_rates[81] -= (fwd_rates[595] - rev_rates[581]);
  //sp 84
  sp_rates[83] += (fwd_rates[595] - rev_rates[581]);
  //sp 31
  sp_rates[30] -= (fwd_rates[595] - rev_rates[581]);
  //sp 32
  sp_rates[31] += (fwd_rates[595] - rev_rates[581]);

  //rxn 596
  //sp 82
  sp_rates[81] -= (fwd_rates[596] - rev_rates[582]);
  //sp 83
  sp_rates[82] += (fwd_rates[596] - rev_rates[582]);
  //sp 31
  sp_rates[30] -= (fwd_rates[596] - rev_rates[582]);
  //sp 32
  sp_rates[31] += (fwd_rates[596] - rev_rates[582]);

  //rxn 597
  //sp 82
  sp_rates[81] -= (fwd_rates[597] - rev_rates[583]);
  //sp 84
  sp_rates[83] += (fwd_rates[597] - rev_rates[583]);
  //sp 14
  sp_rates[13] -= (fwd_rates[597] - rev_rates[583]);
  //sp 15
  sp_rates[14] += (fwd_rates[597] - rev_rates[583]);

  //rxn 598
  //sp 82
  sp_rates[81] -= (fwd_rates[598] - rev_rates[584]);
  //sp 83
  sp_rates[82] += (fwd_rates[598] - rev_rates[584]);
  //sp 14
  sp_rates[13] -= (fwd_rates[598] - rev_rates[584]);
  //sp 15
  sp_rates[14] += (fwd_rates[598] - rev_rates[584]);

  //rxn 599
  //sp 83
  sp_rates[82] += (fwd_rates[599] - rev_rates[585]);
  //sp 84
  sp_rates[83] -= (fwd_rates[599] - rev_rates[585]);

  //rxn 600
  //sp 84
  sp_rates[83] -= (fwd_rates[600] - rev_rates[586]);
  //sp 21
  sp_rates[20] += (fwd_rates[600] - rev_rates[586]);
  //sp 13
  sp_rates[12] += (fwd_rates[600] - rev_rates[586]);

  //rxn 601
  //sp 12
  sp_rates[11] += (fwd_rates[601] - rev_rates[587]);
  //sp 26
  sp_rates[25] += (fwd_rates[601] - rev_rates[587]);
  //sp 84
  sp_rates[83] -= (fwd_rates[601] - rev_rates[587]);

  //rxn 602
  //sp 83
  sp_rates[82] += (fwd_rates[602] - rev_rates[588]);
  //sp 14
  sp_rates[13] -= (fwd_rates[602] - rev_rates[588]);
  //sp 15
  sp_rates[14] -= (fwd_rates[602] - rev_rates[588]);

  //rxn 603
  //sp 79
  sp_rates[78] += (fwd_rates[603] - rev_rates[589]);
  //sp 78
  sp_rates[77] -= (fwd_rates[603] - rev_rates[589]);
  //sp 7
  sp_rates[6] -= (fwd_rates[603] - rev_rates[589]);

  //rxn 604
  //sp 5
  sp_rates[4] += fwd_rates[604];
  //sp 78
  sp_rates[77] -= fwd_rates[604];
  //sp 7
  sp_rates[6] -= fwd_rates[604];
  //sp 15
  sp_rates[14] += 2.0 * fwd_rates[604];

  //rxn 605
  //sp 85
  sp_rates[84] = (fwd_rates[605] - rev_rates[590]);
  //sp 78
  sp_rates[77] -= (fwd_rates[605] - rev_rates[590]);
  //sp 7
  sp_rates[6] -= (fwd_rates[605] - rev_rates[590]);

  //rxn 606
  //sp 7
  sp_rates[6] -= fwd_rates[606];
  //sp 10
  sp_rates[9] += fwd_rates[606];
  //sp 14
  sp_rates[13] += fwd_rates[606];
  //sp 78
  sp_rates[77] -= fwd_rates[606];
  //sp 15
  sp_rates[14] += fwd_rates[606];

  //rxn 607
  //sp 38
  sp_rates[37] += (fwd_rates[607] - rev_rates[591]);
  //sp 37
  sp_rates[36] -= (fwd_rates[607] - rev_rates[591]);
  //sp 78
  sp_rates[77] -= (fwd_rates[607] - rev_rates[591]);
  //sp 77
  sp_rates[76] += (fwd_rates[607] - rev_rates[591]);

  //rxn 608
  //sp 79
  sp_rates[78] -= (fwd_rates[608] - rev_rates[592]);
  //sp 14
  sp_rates[13] += (fwd_rates[608] - rev_rates[592]);
  //sp 15
  sp_rates[14] -= (fwd_rates[608] - rev_rates[592]);
  //sp 80
  sp_rates[79] += (fwd_rates[608] - rev_rates[592]);

  //rxn 609
  //sp 37
  sp_rates[36] -= (fwd_rates[609] - rev_rates[593]);
  //sp 38
  sp_rates[37] += (fwd_rates[609] - rev_rates[593]);
  //sp 79
  sp_rates[78] -= (fwd_rates[609] - rev_rates[593]);
  //sp 80
  sp_rates[79] += (fwd_rates[609] - rev_rates[593]);

  //rxn 610
  //sp 81
  sp_rates[80] += 2.0 * (fwd_rates[610] - rev_rates[594]);
  //sp 79
  sp_rates[78] -= 2.0 * (fwd_rates[610] - rev_rates[594]);
  //sp 7
  sp_rates[6] += (fwd_rates[610] - rev_rates[594]);

  //rxn 611
  //sp 81
  sp_rates[80] += (fwd_rates[611] - rev_rates[595]);
  //sp 5
  sp_rates[4] += (fwd_rates[611] - rev_rates[595]);
  //sp 80
  sp_rates[79] -= (fwd_rates[611] - rev_rates[595]);

  //rxn 612
  //sp 81
  sp_rates[80] -= (fwd_rates[612] - rev_rates[596]);
  //sp 26
  sp_rates[25] += (fwd_rates[612] - rev_rates[596]);
  //sp 15
  sp_rates[14] += (fwd_rates[612] - rev_rates[596]);

  //rxn 613
  //sp 81
  sp_rates[80] -= (fwd_rates[613] - rev_rates[597]);
  //sp 82
  sp_rates[81] += (fwd_rates[613] - rev_rates[597]);
  //sp 7
  sp_rates[6] -= (fwd_rates[613] - rev_rates[597]);
  //sp 8
  sp_rates[7] += (fwd_rates[613] - rev_rates[597]);

  //rxn 614
  //sp 85
  sp_rates[84] += (fwd_rates[614] - rev_rates[598]);
  //sp 79
  sp_rates[78] -= (fwd_rates[614] - rev_rates[598]);

  //rxn 615
  //sp 5
  sp_rates[4] += fwd_rates[615];
  //sp 79
  sp_rates[78] -= fwd_rates[615];
  //sp 15
  sp_rates[14] += 2.0 * fwd_rates[615];

  //rxn 616
  //sp 82
  sp_rates[81] += fwd_rates[616];
  //sp 5
  sp_rates[4] += fwd_rates[616];
  //sp 79
  sp_rates[78] -= fwd_rates[616];

  //rxn 617
  //sp 85
  sp_rates[84] -= fwd_rates[617];
  //sp 5
  sp_rates[4] += fwd_rates[617];
  //sp 15
  sp_rates[14] += 2.0 * fwd_rates[617];

  //rxn 618
  //sp 85
  sp_rates[84] -= (fwd_rates[618] - rev_rates[599]);
  //sp 86
  sp_rates[85] = (fwd_rates[618] - rev_rates[599]);
  //sp 7
  sp_rates[6] -= (fwd_rates[618] - rev_rates[599]);

  //rxn 619
  //sp 87
  sp_rates[86] = (fwd_rates[619] - rev_rates[600]);
  //sp 85
  sp_rates[84] -= (fwd_rates[619] - rev_rates[600]);
  //sp 5
  sp_rates[4] += (fwd_rates[619] - rev_rates[600]);
  //sp 7
  sp_rates[6] -= (fwd_rates[619] - rev_rates[600]);

  //rxn 620
  //sp 82
  sp_rates[81] += (fwd_rates[620] - rev_rates[601]);
  //sp 85
  sp_rates[84] -= (fwd_rates[620] - rev_rates[601]);
  //sp 5
  sp_rates[4] += (fwd_rates[620] - rev_rates[601]);

  //rxn 621
  //sp 5
  sp_rates[4] += (fwd_rates[621] - rev_rates[602]);
  //sp 86
  sp_rates[85] -= (fwd_rates[621] - rev_rates[602]);
  //sp 87
  sp_rates[86] += (fwd_rates[621] - rev_rates[602]);

  //rxn 622
  //sp 5
  sp_rates[4] += (fwd_rates[622] - rev_rates[603]);
  //sp 87
  sp_rates[86] -= (fwd_rates[622] - rev_rates[603]);
  //sp 88
  sp_rates[87] = (fwd_rates[622] - rev_rates[603]);

  //rxn 623
  //sp 25
  sp_rates[24] += (fwd_rates[623] - rev_rates[604]);
  //sp 15
  sp_rates[14] += (fwd_rates[623] - rev_rates[604]);
  //sp 88
  sp_rates[87] -= (fwd_rates[623] - rev_rates[604]);

  //rxn 624
  //sp 89
  sp_rates[88] = fwd_rates[624];
  //sp 88
  sp_rates[87] -= fwd_rates[624];

  //rxn 625
  //sp 14
  sp_rates[13] += fwd_rates[625];
  //sp 23
  sp_rates[22] += fwd_rates[625];
  //sp 88
  sp_rates[87] -= fwd_rates[625];

  //rxn 626
  //sp 24
  sp_rates[23] += fwd_rates[626];
  //sp 13
  sp_rates[12] += fwd_rates[626];
  //sp 88
  sp_rates[87] -= fwd_rates[626];

  //rxn 627
  //sp 89
  sp_rates[88] -= (fwd_rates[627] - rev_rates[605]);
  //sp 90
  sp_rates[89] = (fwd_rates[627] - rev_rates[605]);
  //sp 12
  sp_rates[11] += (fwd_rates[627] - rev_rates[605]);

  //rxn 628
  //sp 89
  sp_rates[88] -= (fwd_rates[628] - rev_rates[606]);
  //sp 13
  sp_rates[12] += (fwd_rates[628] - rev_rates[606]);
  //sp 24
  sp_rates[23] += (fwd_rates[628] - rev_rates[606]);

  //rxn 629
  //sp 89
  sp_rates[88] -= fwd_rates[629];
  //sp 14
  sp_rates[13] += fwd_rates[629];
  //sp 23
  sp_rates[22] += fwd_rates[629];

  //rxn 630
  //sp 90
  sp_rates[89] -= (fwd_rates[630] - rev_rates[607]);
  //sp 3
  sp_rates[2] += (fwd_rates[630] - rev_rates[607]);
  //sp 23
  sp_rates[22] += (fwd_rates[630] - rev_rates[607]);

  //rxn 631
  //sp 90
  sp_rates[89] += (fwd_rates[631] - rev_rates[608]);
  //sp 5
  sp_rates[4] -= (fwd_rates[631] - rev_rates[608]);
  //sp 15
  sp_rates[14] -= (fwd_rates[631] - rev_rates[608]);

  //rxn 632
  //sp 90
  sp_rates[89] += (fwd_rates[632] - rev_rates[609]);
  //sp 8
  sp_rates[7] -= (fwd_rates[632] - rev_rates[609]);
  //sp 5
  sp_rates[4] += (fwd_rates[632] - rev_rates[609]);
  //sp 24
  sp_rates[23] -= (fwd_rates[632] - rev_rates[609]);

  //rxn 633
  //sp 11
  sp_rates[10] = (fwd_rates[633] - rev_rates[610]) * pres_mod[74];
  //sp 3
  sp_rates[2] -= (fwd_rates[633] - rev_rates[610]) * pres_mod[74];
  //sp 4
  sp_rates[3] -= (fwd_rates[633] - rev_rates[610]) * pres_mod[74];

  //rxn 634
  //sp 11
  sp_rates[10] += (fwd_rates[634] - rev_rates[611]);
  //sp 12
  sp_rates[11] += (fwd_rates[634] - rev_rates[611]);
  //sp 7
  sp_rates[6] -= (fwd_rates[634] - rev_rates[611]);
  //sp 16
  sp_rates[15] -= (fwd_rates[634] - rev_rates[611]);

  //rxn 635
  //sp 12
  sp_rates[11] += (fwd_rates[635] - rev_rates[612]);
  //sp 11
  sp_rates[10] += (fwd_rates[635] - rev_rates[612]);
  //sp 4
  sp_rates[3] -= (fwd_rates[635] - rev_rates[612]);
  //sp 14
  sp_rates[13] -= (fwd_rates[635] - rev_rates[612]);

  //rxn 636
  //sp 11
  sp_rates[10] -= (fwd_rates[636] - rev_rates[613]);
  //sp 5
  sp_rates[4] += (fwd_rates[636] - rev_rates[613]);

  //rxn 637
  //sp 11
  sp_rates[10] -= (fwd_rates[637] - rev_rates[614]);
  //sp 5
  sp_rates[4] += (fwd_rates[637] - rev_rates[614]);

  //rxn 638
  //sp 11
  sp_rates[10] -= (fwd_rates[638] - rev_rates[615]);
  //sp 5
  sp_rates[4] += (fwd_rates[638] - rev_rates[615]);

  //rxn 639
  //sp 11
  sp_rates[10] -= (fwd_rates[639] - rev_rates[616]);
  //sp 5
  sp_rates[4] += (fwd_rates[639] - rev_rates[616]);

  //rxn 640
  //sp 11
  sp_rates[10] -= (fwd_rates[640] - rev_rates[617]);
  //sp 5
  sp_rates[4] += (fwd_rates[640] - rev_rates[617]);

  //rxn 641
  //sp 11
  sp_rates[10] -= (fwd_rates[641] - rev_rates[618]);
  //sp 5
  sp_rates[4] += (fwd_rates[641] - rev_rates[618]);

  //rxn 642
  //sp 11
  sp_rates[10] -= (fwd_rates[642] - rev_rates[619]);
  //sp 5
  sp_rates[4] += (fwd_rates[642] - rev_rates[619]);

  //rxn 643
  //sp 11
  sp_rates[10] -= (fwd_rates[643] - rev_rates[620]);
  //sp 5
  sp_rates[4] += (fwd_rates[643] - rev_rates[620]);

  //rxn 644
  //sp 11
  sp_rates[10] -= (fwd_rates[644] - rev_rates[621]);
  //sp 5
  sp_rates[4] += (fwd_rates[644] - rev_rates[621]);

  //rxn 645
  //sp 11
  sp_rates[10] -= (fwd_rates[645] - rev_rates[622]);
  //sp 5
  sp_rates[4] += (fwd_rates[645] - rev_rates[622]);

  //rxn 646
  //sp 11
  sp_rates[10] -= (fwd_rates[646] - rev_rates[623]);
  //sp 5
  sp_rates[4] += (fwd_rates[646] - rev_rates[623]);

  //rxn 647
  //sp 91
  sp_rates[90] = (fwd_rates[647] - rev_rates[624]);
  //sp 13
  sp_rates[12] += (fwd_rates[647] - rev_rates[624]);
  //sp 7
  sp_rates[6] -= (fwd_rates[647] - rev_rates[624]);
  //sp 40
  sp_rates[39] -= (fwd_rates[647] - rev_rates[624]);

  //rxn 648
  //sp 12
  sp_rates[11] += (fwd_rates[648] - rev_rates[625]);
  //sp 91
  sp_rates[90] += (fwd_rates[648] - rev_rates[625]);
  //sp 4
  sp_rates[3] -= (fwd_rates[648] - rev_rates[625]);
  //sp 40
  sp_rates[39] -= (fwd_rates[648] - rev_rates[625]);

  //rxn 649
  //sp 41
  sp_rates[40] -= (fwd_rates[649] - rev_rates[626]);
  //sp 91
  sp_rates[90] += (fwd_rates[649] - rev_rates[626]);
  //sp 12
  sp_rates[11] += (fwd_rates[649] - rev_rates[626]);
  //sp 5
  sp_rates[4] -= (fwd_rates[649] - rev_rates[626]);

  //rxn 650
  //sp 91
  sp_rates[90] -= (fwd_rates[650] - rev_rates[627]);
  //sp 16
  sp_rates[15] += (fwd_rates[650] - rev_rates[627]);

  //rxn 651
  //sp 91
  sp_rates[90] -= (fwd_rates[651] - rev_rates[628]);
  //sp 16
  sp_rates[15] += (fwd_rates[651] - rev_rates[628]);

  //rxn 652
  //sp 91
  sp_rates[90] -= (fwd_rates[652] - rev_rates[629]);
  //sp 16
  sp_rates[15] += (fwd_rates[652] - rev_rates[629]);

  //rxn 653
  //sp 91
  sp_rates[90] -= (fwd_rates[653] - rev_rates[630]);
  //sp 16
  sp_rates[15] += (fwd_rates[653] - rev_rates[630]);

  //rxn 654
  //sp 91
  sp_rates[90] -= (fwd_rates[654] - rev_rates[631]);
  //sp 16
  sp_rates[15] += (fwd_rates[654] - rev_rates[631]);

  //rxn 655
  //sp 91
  sp_rates[90] -= (fwd_rates[655] - rev_rates[632]);
  //sp 16
  sp_rates[15] += (fwd_rates[655] - rev_rates[632]);

  //rxn 656
  //sp 91
  sp_rates[90] -= (fwd_rates[656] - rev_rates[633]);
  //sp 16
  sp_rates[15] += (fwd_rates[656] - rev_rates[633]);

  //rxn 657
  //sp 91
  sp_rates[90] -= (fwd_rates[657] - rev_rates[634]);
  //sp 16
  sp_rates[15] += (fwd_rates[657] - rev_rates[634]);

  //rxn 658
  //sp 91
  sp_rates[90] -= (fwd_rates[658] - rev_rates[635]);
  //sp 16
  sp_rates[15] += (fwd_rates[658] - rev_rates[635]);

  //rxn 659
  //sp 91
  sp_rates[90] -= (fwd_rates[659] - rev_rates[636]);
  //sp 16
  sp_rates[15] += (fwd_rates[659] - rev_rates[636]);

  //rxn 660
  //sp 91
  sp_rates[90] -= (fwd_rates[660] - rev_rates[637]);
  //sp 16
  sp_rates[15] += (fwd_rates[660] - rev_rates[637]);

  //rxn 661
  //sp 92
  sp_rates[91] = 2.0 * (fwd_rates[661] - rev_rates[638]) * pres_mod[75];
  //sp 2
  (*dy_N) = -(fwd_rates[661] - rev_rates[638]) * pres_mod[75];

  //rxn 662
  //sp 4
  sp_rates[3] -= (fwd_rates[662] - rev_rates[639]) * pres_mod[76];
  //sp 92
  sp_rates[91] -= (fwd_rates[662] - rev_rates[639]) * pres_mod[76];
  //sp 93
  sp_rates[92] = (fwd_rates[662] - rev_rates[639]) * pres_mod[76];

  //rxn 663
  //sp 4
  sp_rates[3] += (fwd_rates[663] - rev_rates[640]);
  //sp 92
  sp_rates[91] -= (fwd_rates[663] - rev_rates[640]);
  //sp 94
  sp_rates[93] = (fwd_rates[663] - rev_rates[640]);
  //sp 95
  sp_rates[94] = -(fwd_rates[663] - rev_rates[640]);

  //rxn 664
  //sp 4
  sp_rates[3] += (fwd_rates[664] - rev_rates[641]);
  //sp 92
  sp_rates[91] -= (fwd_rates[664] - rev_rates[641]);
  //sp 93
  sp_rates[92] += (fwd_rates[664] - rev_rates[641]);
  //sp 7
  sp_rates[6] -= (fwd_rates[664] - rev_rates[641]);

  //rxn 665
  //sp 3
  sp_rates[2] += (fwd_rates[665] - rev_rates[642]);
  //sp 92
  sp_rates[91] -= (fwd_rates[665] - rev_rates[642]);
  //sp 5
  sp_rates[4] -= (fwd_rates[665] - rev_rates[642]);
  //sp 93
  sp_rates[92] += (fwd_rates[665] - rev_rates[642]);

  //rxn 666
  //sp 2
  (*dy_N) += (fwd_rates[666] - rev_rates[643]);
  //sp 4
  sp_rates[3] += (fwd_rates[666] - rev_rates[643]);
  //sp 92
  sp_rates[91] -= (fwd_rates[666] - rev_rates[643]);
  //sp 93
  sp_rates[92] -= (fwd_rates[666] - rev_rates[643]);

  //rxn 667
  //sp 3
  sp_rates[2] -= (fwd_rates[667] - rev_rates[644]);
  //sp 92
  sp_rates[91] += (fwd_rates[667] - rev_rates[644]);
  //sp 6
  sp_rates[5] += (fwd_rates[667] - rev_rates[644]);
  //sp 96
  sp_rates[95] = -(fwd_rates[667] - rev_rates[644]);

  //rxn 668
  //sp 3
  sp_rates[2] += (fwd_rates[668] - rev_rates[645]);
  //sp 4
  sp_rates[3] -= (fwd_rates[668] - rev_rates[645]);
  //sp 93
  sp_rates[92] += (fwd_rates[668] - rev_rates[645]);
  //sp 96
  sp_rates[95] -= (fwd_rates[668] - rev_rates[645]);

  //rxn 669
  //sp 92
  sp_rates[91] += (fwd_rates[669] - rev_rates[646]);
  //sp 4
  sp_rates[3] -= (fwd_rates[669] - rev_rates[646]);
  //sp 5
  sp_rates[4] += (fwd_rates[669] - rev_rates[646]);
  //sp 96
  sp_rates[95] -= (fwd_rates[669] - rev_rates[646]);

  //rxn 670
  //sp 3
  sp_rates[2] += (fwd_rates[670] - rev_rates[647]);
  //sp 5
  sp_rates[4] -= (fwd_rates[670] - rev_rates[647]);
  //sp 99
  sp_rates[98] = (fwd_rates[670] - rev_rates[647]);
  //sp 96
  sp_rates[95] -= (fwd_rates[670] - rev_rates[647]);

  //rxn 671
  //sp 10
  sp_rates[9] += (fwd_rates[671] - rev_rates[648]);
  //sp 92
  sp_rates[91] += (fwd_rates[671] - rev_rates[648]);
  //sp 5
  sp_rates[4] -= (fwd_rates[671] - rev_rates[648]);
  //sp 96
  sp_rates[95] -= (fwd_rates[671] - rev_rates[648]);

  //rxn 672
  //sp 99
  sp_rates[98] += (fwd_rates[672] - rev_rates[649]);
  //sp 4
  sp_rates[3] += (fwd_rates[672] - rev_rates[649]);
  //sp 7
  sp_rates[6] -= (fwd_rates[672] - rev_rates[649]);
  //sp 96
  sp_rates[95] -= (fwd_rates[672] - rev_rates[649]);

  //rxn 673
  //sp 93
  sp_rates[92] += (fwd_rates[673] - rev_rates[650]);
  //sp 5
  sp_rates[4] += (fwd_rates[673] - rev_rates[650]);
  //sp 7
  sp_rates[6] -= (fwd_rates[673] - rev_rates[650]);
  //sp 96
  sp_rates[95] -= (fwd_rates[673] - rev_rates[650]);

  //rxn 674
  //sp 3
  sp_rates[2] += 2.0 * fwd_rates[674];
  //sp 2
  (*dy_N) += fwd_rates[674];
  //sp 96
  sp_rates[95] -= 2.0 * fwd_rates[674];

  //rxn 675
  //sp 3
  sp_rates[2] += (fwd_rates[675] - rev_rates[651]);
  //sp 92
  sp_rates[91] -= (fwd_rates[675] - rev_rates[651]);
  //sp 2
  (*dy_N) += (fwd_rates[675] - rev_rates[651]);
  //sp 96
  sp_rates[95] -= (fwd_rates[675] - rev_rates[651]);

  //rxn 676
  //sp 3
  sp_rates[2] += (fwd_rates[676] - rev_rates[652]);
  //sp 93
  sp_rates[92] -= (fwd_rates[676] - rev_rates[652]);
  //sp 94
  sp_rates[93] += (fwd_rates[676] - rev_rates[652]);
  //sp 96
  sp_rates[95] -= (fwd_rates[676] - rev_rates[652]);

  //rxn 677
  //sp 2
  (*dy_N) += (fwd_rates[677] - rev_rates[653]);
  //sp 93
  sp_rates[92] -= (fwd_rates[677] - rev_rates[653]);
  //sp 5
  sp_rates[4] += (fwd_rates[677] - rev_rates[653]);
  //sp 96
  sp_rates[95] -= (fwd_rates[677] - rev_rates[653]);

  //rxn 678
  //sp 5
  sp_rates[4] += (fwd_rates[678] - rev_rates[654]);
  //sp 94
  sp_rates[93] += (fwd_rates[678] - rev_rates[654]);
  //sp 95
  sp_rates[94] -= (fwd_rates[678] - rev_rates[654]);
  //sp 96
  sp_rates[95] -= (fwd_rates[678] - rev_rates[654]);

  //rxn 679
  //sp 99
  sp_rates[98] += (fwd_rates[679] - rev_rates[655]);
  //sp 93
  sp_rates[92] += (fwd_rates[679] - rev_rates[655]);
  //sp 95
  sp_rates[94] -= (fwd_rates[679] - rev_rates[655]);
  //sp 96
  sp_rates[95] -= (fwd_rates[679] - rev_rates[655]);

  //rxn 680
  //sp 97
  sp_rates[96] = (fwd_rates[680] - rev_rates[656]);
  //sp 100
  sp_rates[99] = -(fwd_rates[680] - rev_rates[656]);
  //sp 95
  sp_rates[94] += (fwd_rates[680] - rev_rates[656]);
  //sp 96
  sp_rates[95] -= (fwd_rates[680] - rev_rates[656]);

  //rxn 681
  //sp 99
  sp_rates[98] += (fwd_rates[681] - rev_rates[657]);
  //sp 2
  (*dy_N) += (fwd_rates[681] - rev_rates[657]);
  //sp 94
  sp_rates[93] -= (fwd_rates[681] - rev_rates[657]);
  //sp 96
  sp_rates[95] -= (fwd_rates[681] - rev_rates[657]);

  //rxn 682
  //sp 97
  sp_rates[96] -= (fwd_rates[682] - rev_rates[658]);
  //sp 3
  sp_rates[2] -= (fwd_rates[682] - rev_rates[658]);
  //sp 6
  sp_rates[5] += (fwd_rates[682] - rev_rates[658]);
  //sp 96
  sp_rates[95] += (fwd_rates[682] - rev_rates[658]);

  //rxn 683
  //sp 97
  sp_rates[96] -= (fwd_rates[683] - rev_rates[659]);
  //sp 3
  sp_rates[2] += (fwd_rates[683] - rev_rates[659]);
  //sp 4
  sp_rates[3] -= (fwd_rates[683] - rev_rates[659]);
  //sp 99
  sp_rates[98] += (fwd_rates[683] - rev_rates[659]);

  //rxn 684
  //sp 97
  sp_rates[96] -= (fwd_rates[684] - rev_rates[660]);
  //sp 3
  sp_rates[2] += (fwd_rates[684] - rev_rates[660]);
  //sp 4
  sp_rates[3] -= (fwd_rates[684] - rev_rates[660]);
  //sp 99
  sp_rates[98] += (fwd_rates[684] - rev_rates[660]);

  //rxn 685
  //sp 97
  sp_rates[96] -= (fwd_rates[685] - rev_rates[661]);
  //sp 4
  sp_rates[3] -= (fwd_rates[685] - rev_rates[661]);
  //sp 5
  sp_rates[4] += (fwd_rates[685] - rev_rates[661]);
  //sp 96
  sp_rates[95] += (fwd_rates[685] - rev_rates[661]);

  //rxn 686
  //sp 97
  sp_rates[96] -= (fwd_rates[686] - rev_rates[662]);
  //sp 4
  sp_rates[3] -= (fwd_rates[686] - rev_rates[662]);
  //sp 5
  sp_rates[4] += (fwd_rates[686] - rev_rates[662]);
  //sp 96
  sp_rates[95] += (fwd_rates[686] - rev_rates[662]);

  //rxn 687
  //sp 97
  sp_rates[96] -= (fwd_rates[687] - rev_rates[663]);
  //sp 4
  sp_rates[3] -= (fwd_rates[687] - rev_rates[663]);
  //sp 93
  sp_rates[92] += (fwd_rates[687] - rev_rates[663]);
  //sp 6
  sp_rates[5] += (fwd_rates[687] - rev_rates[663]);

  //rxn 688
  //sp 97
  sp_rates[96] -= (fwd_rates[688] - rev_rates[664]);
  //sp 10
  sp_rates[9] += (fwd_rates[688] - rev_rates[664]);
  //sp 5
  sp_rates[4] -= (fwd_rates[688] - rev_rates[664]);
  //sp 96
  sp_rates[95] += (fwd_rates[688] - rev_rates[664]);

  //rxn 689
  //sp 97
  sp_rates[96] -= (fwd_rates[689] - rev_rates[665]);
  //sp 98
  sp_rates[97] = (fwd_rates[689] - rev_rates[665]);
  //sp 7
  sp_rates[6] += (fwd_rates[689] - rev_rates[665]);
  //sp 8
  sp_rates[7] -= (fwd_rates[689] - rev_rates[665]);

  //rxn 690
  //sp 97
  sp_rates[96] -= (fwd_rates[690] - rev_rates[666]);
  //sp 101
  sp_rates[100] = (fwd_rates[690] - rev_rates[666]);
  //sp 5
  sp_rates[4] += (fwd_rates[690] - rev_rates[666]);
  //sp 8
  sp_rates[7] -= (fwd_rates[690] - rev_rates[666]);

  //rxn 691
  //sp 97
  sp_rates[96] -= (fwd_rates[691] - rev_rates[667]);
  //sp 10
  sp_rates[9] += (fwd_rates[691] - rev_rates[667]);
  //sp 99
  sp_rates[98] += (fwd_rates[691] - rev_rates[667]);
  //sp 8
  sp_rates[7] -= (fwd_rates[691] - rev_rates[667]);

  //rxn 692
  //sp 97
  sp_rates[96] -= (fwd_rates[692] - rev_rates[668]);
  //sp 10
  sp_rates[9] += (fwd_rates[692] - rev_rates[668]);
  //sp 99
  sp_rates[98] += (fwd_rates[692] - rev_rates[668]);
  //sp 8
  sp_rates[7] -= (fwd_rates[692] - rev_rates[668]);

  //rxn 693
  //sp 97
  sp_rates[96] -= (fwd_rates[693] - rev_rates[669]);
  //sp 10
  sp_rates[9] += (fwd_rates[693] - rev_rates[669]);
  //sp 123
  sp_rates[122] = (fwd_rates[693] - rev_rates[669]);
  //sp 8
  sp_rates[7] -= (fwd_rates[693] - rev_rates[669]);

  //rxn 694
  //sp 97
  sp_rates[96] -= (fwd_rates[694] - rev_rates[670]);
  //sp 4
  sp_rates[3] += (fwd_rates[694] - rev_rates[670]);
  //sp 101
  sp_rates[100] += (fwd_rates[694] - rev_rates[670]);
  //sp 7
  sp_rates[6] -= (fwd_rates[694] - rev_rates[670]);

  //rxn 695
  //sp 97
  sp_rates[96] -= (fwd_rates[695] - rev_rates[671]);
  //sp 99
  sp_rates[98] += (fwd_rates[695] - rev_rates[671]);
  //sp 5
  sp_rates[4] += (fwd_rates[695] - rev_rates[671]);
  //sp 7
  sp_rates[6] -= (fwd_rates[695] - rev_rates[671]);

  //rxn 696
  //sp 97
  sp_rates[96] -= 2.0 * (fwd_rates[696] - rev_rates[672]);
  //sp 98
  sp_rates[97] += (fwd_rates[696] - rev_rates[672]);
  //sp 96
  sp_rates[95] += (fwd_rates[696] - rev_rates[672]);

  //rxn 697
  //sp 97
  sp_rates[96] -= (fwd_rates[697] - rev_rates[673]);
  //sp 98
  sp_rates[97] += (fwd_rates[697] - rev_rates[673]);
  //sp 92
  sp_rates[91] += (fwd_rates[697] - rev_rates[673]);
  //sp 96
  sp_rates[95] -= (fwd_rates[697] - rev_rates[673]);

  //rxn 698
  //sp 97
  sp_rates[96] -= fwd_rates[698];
  //sp 3
  sp_rates[2] += 2.0 * fwd_rates[698];
  //sp 92
  sp_rates[91] -= fwd_rates[698];
  //sp 2
  (*dy_N) += fwd_rates[698];

  //rxn 699
  //sp 97
  sp_rates[96] -= (fwd_rates[699] - rev_rates[674]);
  //sp 10
  sp_rates[9] += (fwd_rates[699] - rev_rates[674]);
  //sp 93
  sp_rates[92] -= (fwd_rates[699] - rev_rates[674]);
  //sp 2
  (*dy_N) += (fwd_rates[699] - rev_rates[674]);

  //rxn 700
  //sp 97
  sp_rates[96] -= (fwd_rates[700] - rev_rates[675]);
  //sp 5
  sp_rates[4] += (fwd_rates[700] - rev_rates[675]);
  //sp 93
  sp_rates[92] -= (fwd_rates[700] - rev_rates[675]);
  //sp 102
  sp_rates[101] = (fwd_rates[700] - rev_rates[675]);

  //rxn 701
  //sp 97
  sp_rates[96] -= (fwd_rates[701] - rev_rates[676]);
  //sp 10
  sp_rates[9] += (fwd_rates[701] - rev_rates[676]);
  //sp 94
  sp_rates[93] += (fwd_rates[701] - rev_rates[676]);
  //sp 95
  sp_rates[94] -= (fwd_rates[701] - rev_rates[676]);

  //rxn 702
  //sp 97
  sp_rates[96] -= (fwd_rates[702] - rev_rates[677]);
  //sp 101
  sp_rates[100] += (fwd_rates[702] - rev_rates[677]);
  //sp 95
  sp_rates[94] -= (fwd_rates[702] - rev_rates[677]);
  //sp 93
  sp_rates[92] += (fwd_rates[702] - rev_rates[677]);

  //rxn 703
  //sp 97
  sp_rates[96] -= (fwd_rates[703] - rev_rates[678]);
  //sp 98
  sp_rates[97] += (fwd_rates[703] - rev_rates[678]);
  //sp 99
  sp_rates[98] -= (fwd_rates[703] - rev_rates[678]);
  //sp 93
  sp_rates[92] += (fwd_rates[703] - rev_rates[678]);

  //rxn 704
  //sp 97
  sp_rates[96] -= (fwd_rates[704] - rev_rates[679]);
  //sp 98
  sp_rates[97] += (fwd_rates[704] - rev_rates[679]);
  //sp 100
  sp_rates[99] -= (fwd_rates[704] - rev_rates[679]);
  //sp 95
  sp_rates[94] += (fwd_rates[704] - rev_rates[679]);

  //rxn 705
  //sp 97
  sp_rates[96] += (fwd_rates[705] - rev_rates[680]);
  //sp 92
  sp_rates[91] += (fwd_rates[705] - rev_rates[680]);
  //sp 96
  sp_rates[95] -= 2.0 * (fwd_rates[705] - rev_rates[680]);

  //rxn 706
  //sp 3
  sp_rates[2] += (fwd_rates[706] - rev_rates[681]) * pres_mod[77];
  //sp 102
  sp_rates[101] += (fwd_rates[706] - rev_rates[681]) * pres_mod[77];
  //sp 103
  sp_rates[102] = -(fwd_rates[706] - rev_rates[681]) * pres_mod[77];

  //rxn 707
  //sp 3
  sp_rates[2] -= (fwd_rates[707] - rev_rates[682]);
  //sp 102
  sp_rates[101] += (fwd_rates[707] - rev_rates[682]);
  //sp 6
  sp_rates[5] += (fwd_rates[707] - rev_rates[682]);
  //sp 103
  sp_rates[102] -= (fwd_rates[707] - rev_rates[682]);

  //rxn 708
  //sp 97
  sp_rates[96] += (fwd_rates[708] - rev_rates[683]);
  //sp 4
  sp_rates[3] -= (fwd_rates[708] - rev_rates[683]);
  //sp 93
  sp_rates[92] += (fwd_rates[708] - rev_rates[683]);
  //sp 103
  sp_rates[102] -= (fwd_rates[708] - rev_rates[683]);

  //rxn 709
  //sp 4
  sp_rates[3] -= (fwd_rates[709] - rev_rates[684]);
  //sp 5
  sp_rates[4] += (fwd_rates[709] - rev_rates[684]);
  //sp 102
  sp_rates[101] += (fwd_rates[709] - rev_rates[684]);
  //sp 103
  sp_rates[102] -= (fwd_rates[709] - rev_rates[684]);

  //rxn 710
  //sp 10
  sp_rates[9] += (fwd_rates[710] - rev_rates[685]);
  //sp 5
  sp_rates[4] -= (fwd_rates[710] - rev_rates[685]);
  //sp 102
  sp_rates[101] += (fwd_rates[710] - rev_rates[685]);
  //sp 103
  sp_rates[102] -= (fwd_rates[710] - rev_rates[685]);

  //rxn 711
  //sp 97
  sp_rates[96] += (fwd_rates[711] - rev_rates[686]);
  //sp 93
  sp_rates[92] -= (fwd_rates[711] - rev_rates[686]);
  //sp 94
  sp_rates[93] += (fwd_rates[711] - rev_rates[686]);
  //sp 103
  sp_rates[102] -= (fwd_rates[711] - rev_rates[686]);

  //rxn 712
  //sp 97
  sp_rates[96] += (fwd_rates[712] - rev_rates[687]);
  //sp 102
  sp_rates[101] += (fwd_rates[712] - rev_rates[687]);
  //sp 103
  sp_rates[102] -= (fwd_rates[712] - rev_rates[687]);
  //sp 96
  sp_rates[95] -= (fwd_rates[712] - rev_rates[687]);

  //rxn 713
  //sp 97
  sp_rates[96] -= (fwd_rates[713] - rev_rates[688]);
  //sp 98
  sp_rates[97] += (fwd_rates[713] - rev_rates[688]);
  //sp 102
  sp_rates[101] += (fwd_rates[713] - rev_rates[688]);
  //sp 103
  sp_rates[102] -= (fwd_rates[713] - rev_rates[688]);

  //rxn 714
  //sp 97
  sp_rates[96] -= (fwd_rates[714] - rev_rates[689]);
  //sp 3
  sp_rates[2] += (fwd_rates[714] - rev_rates[689]);
  //sp 103
  sp_rates[102] += (fwd_rates[714] - rev_rates[689]);
  //sp 96
  sp_rates[95] -= (fwd_rates[714] - rev_rates[689]);

  //rxn 715
  //sp 97
  sp_rates[96] -= 2.0 * (fwd_rates[715] - rev_rates[690]);
  //sp 6
  sp_rates[5] += (fwd_rates[715] - rev_rates[690]);
  //sp 103
  sp_rates[102] += (fwd_rates[715] - rev_rates[690]);

  //rxn 716
  //sp 9
  sp_rates[8] += (fwd_rates[716] - rev_rates[691]);
  //sp 102
  sp_rates[101] += (fwd_rates[716] - rev_rates[691]);
  //sp 103
  sp_rates[102] -= (fwd_rates[716] - rev_rates[691]);
  //sp 8
  sp_rates[7] -= (fwd_rates[716] - rev_rates[691]);

  //rxn 717
  //sp 103
  sp_rates[102] -= (fwd_rates[717] - rev_rates[692]) * pres_mod[78];
  //sp 96
  sp_rates[95] += 2.0 * (fwd_rates[717] - rev_rates[692]) * pres_mod[78];

  //rxn 718
  //sp 97
  sp_rates[96] -= 2.0 * (fwd_rates[718] - rev_rates[693]) * pres_mod[79];
  //sp 105
  sp_rates[104] = (fwd_rates[718] - rev_rates[693]) * pres_mod[79];

  //rxn 719
  //sp 105
  sp_rates[104] -= (fwd_rates[719] - rev_rates[694]) * pres_mod[80];
  //sp 3
  sp_rates[2] += (fwd_rates[719] - rev_rates[694]) * pres_mod[80];
  //sp 104
  sp_rates[103] = (fwd_rates[719] - rev_rates[694]) * pres_mod[80];

  //rxn 720
  //sp 105
  sp_rates[104] -= (fwd_rates[720] - rev_rates[695]);
  //sp 3
  sp_rates[2] -= (fwd_rates[720] - rev_rates[695]);
  //sp 6
  sp_rates[5] += (fwd_rates[720] - rev_rates[695]);
  //sp 104
  sp_rates[103] += (fwd_rates[720] - rev_rates[695]);

  //rxn 721
  //sp 105
  sp_rates[104] -= (fwd_rates[721] - rev_rates[696]);
  //sp 10
  sp_rates[9] += (fwd_rates[721] - rev_rates[696]);
  //sp 4
  sp_rates[3] -= (fwd_rates[721] - rev_rates[696]);
  //sp 103
  sp_rates[102] += (fwd_rates[721] - rev_rates[696]);

  //rxn 722
  //sp 105
  sp_rates[104] -= (fwd_rates[722] - rev_rates[697]);
  //sp 4
  sp_rates[3] -= (fwd_rates[722] - rev_rates[697]);
  //sp 5
  sp_rates[4] += (fwd_rates[722] - rev_rates[697]);
  //sp 104
  sp_rates[103] += (fwd_rates[722] - rev_rates[697]);

  //rxn 723
  //sp 105
  sp_rates[104] -= (fwd_rates[723] - rev_rates[698]);
  //sp 10
  sp_rates[9] += (fwd_rates[723] - rev_rates[698]);
  //sp 5
  sp_rates[4] -= (fwd_rates[723] - rev_rates[698]);
  //sp 104
  sp_rates[103] += (fwd_rates[723] - rev_rates[698]);

  //rxn 724
  //sp 105
  sp_rates[104] -= (fwd_rates[724] - rev_rates[699]);
  //sp 98
  sp_rates[97] += (fwd_rates[724] - rev_rates[699]);
  //sp 97
  sp_rates[96] -= (fwd_rates[724] - rev_rates[699]);
  //sp 104
  sp_rates[103] += (fwd_rates[724] - rev_rates[699]);

  //rxn 725
  //sp 105
  sp_rates[104] -= (fwd_rates[725] - rev_rates[700]);
  //sp 98
  sp_rates[97] += (fwd_rates[725] - rev_rates[700]);
  //sp 3
  sp_rates[2] -= (fwd_rates[725] - rev_rates[700]);
  //sp 97
  sp_rates[96] += (fwd_rates[725] - rev_rates[700]);

  //rxn 726
  //sp 105
  sp_rates[104] -= (fwd_rates[726] - rev_rates[701]);
  //sp 103
  sp_rates[102] -= (fwd_rates[726] - rev_rates[701]);
  //sp 104
  sp_rates[103] += 2.0 * (fwd_rates[726] - rev_rates[701]);

  //rxn 727
  //sp 105
  sp_rates[104] -= (fwd_rates[727] - rev_rates[702]);
  //sp 9
  sp_rates[8] += (fwd_rates[727] - rev_rates[702]);
  //sp 104
  sp_rates[103] += (fwd_rates[727] - rev_rates[702]);
  //sp 8
  sp_rates[7] -= (fwd_rates[727] - rev_rates[702]);

  //rxn 728
  //sp 105
  sp_rates[104] -= (fwd_rates[728] - rev_rates[703]);
  //sp 104
  sp_rates[103] += (fwd_rates[728] - rev_rates[703]);
  //sp 97
  sp_rates[96] += (fwd_rates[728] - rev_rates[703]);
  //sp 96
  sp_rates[95] -= (fwd_rates[728] - rev_rates[703]);

  //rxn 729
  //sp 3
  sp_rates[2] += (fwd_rates[729] - rev_rates[704]);
  //sp 103
  sp_rates[102] += (fwd_rates[729] - rev_rates[704]);
  //sp 104
  sp_rates[103] -= (fwd_rates[729] - rev_rates[704]);

  //rxn 730
  //sp 97
  sp_rates[96] += 2.0 * (fwd_rates[730] - rev_rates[705]);
  //sp 3
  sp_rates[2] -= (fwd_rates[730] - rev_rates[705]);
  //sp 104
  sp_rates[103] -= (fwd_rates[730] - rev_rates[705]);

  //rxn 731
  //sp 3
  sp_rates[2] -= (fwd_rates[731] - rev_rates[706]);
  //sp 6
  sp_rates[5] += (fwd_rates[731] - rev_rates[706]);
  //sp 103
  sp_rates[102] += (fwd_rates[731] - rev_rates[706]);
  //sp 104
  sp_rates[103] -= (fwd_rates[731] - rev_rates[706]);

  //rxn 732
  //sp 98
  sp_rates[97] += (fwd_rates[732] - rev_rates[707]);
  //sp 3
  sp_rates[2] -= (fwd_rates[732] - rev_rates[707]);
  //sp 96
  sp_rates[95] += (fwd_rates[732] - rev_rates[707]);
  //sp 104
  sp_rates[103] -= (fwd_rates[732] - rev_rates[707]);

  //rxn 733
  //sp 4
  sp_rates[3] -= (fwd_rates[733] - rev_rates[708]);
  //sp 5
  sp_rates[4] += (fwd_rates[733] - rev_rates[708]);
  //sp 103
  sp_rates[102] += (fwd_rates[733] - rev_rates[708]);
  //sp 104
  sp_rates[103] -= (fwd_rates[733] - rev_rates[708]);

  //rxn 734
  //sp 97
  sp_rates[96] += (fwd_rates[734] - rev_rates[709]);
  //sp 99
  sp_rates[98] += (fwd_rates[734] - rev_rates[709]);
  //sp 4
  sp_rates[3] -= (fwd_rates[734] - rev_rates[709]);
  //sp 104
  sp_rates[103] -= (fwd_rates[734] - rev_rates[709]);

  //rxn 735
  //sp 97
  sp_rates[96] += fwd_rates[735];
  //sp 3
  sp_rates[2] += fwd_rates[735];
  //sp 4
  sp_rates[3] -= fwd_rates[735];
  //sp 104
  sp_rates[103] -= fwd_rates[735];
  //sp 93
  sp_rates[92] += fwd_rates[735];

  //rxn 736
  //sp 10
  sp_rates[9] += (fwd_rates[736] - rev_rates[710]);
  //sp 5
  sp_rates[4] -= (fwd_rates[736] - rev_rates[710]);
  //sp 103
  sp_rates[102] += (fwd_rates[736] - rev_rates[710]);
  //sp 104
  sp_rates[103] -= (fwd_rates[736] - rev_rates[710]);

  //rxn 737
  //sp 106
  sp_rates[105] = (fwd_rates[737] - rev_rates[711]);
  //sp 5
  sp_rates[4] -= (fwd_rates[737] - rev_rates[711]);
  //sp 10
  sp_rates[9] += (fwd_rates[737] - rev_rates[711]);
  //sp 104
  sp_rates[103] -= (fwd_rates[737] - rev_rates[711]);

  //rxn 738
  //sp 98
  sp_rates[97] += (fwd_rates[738] - rev_rates[712]);
  //sp 99
  sp_rates[98] += (fwd_rates[738] - rev_rates[712]);
  //sp 5
  sp_rates[4] -= (fwd_rates[738] - rev_rates[712]);
  //sp 104
  sp_rates[103] -= (fwd_rates[738] - rev_rates[712]);

  //rxn 739
  //sp 9
  sp_rates[8] += (fwd_rates[739] - rev_rates[713]);
  //sp 104
  sp_rates[103] -= (fwd_rates[739] - rev_rates[713]);
  //sp 103
  sp_rates[102] += (fwd_rates[739] - rev_rates[713]);
  //sp 8
  sp_rates[7] -= (fwd_rates[739] - rev_rates[713]);

  //rxn 740
  //sp 105
  sp_rates[104] += (fwd_rates[740] - rev_rates[714]);
  //sp 104
  sp_rates[103] -= (fwd_rates[740] - rev_rates[714]);
  //sp 7
  sp_rates[6] += (fwd_rates[740] - rev_rates[714]);
  //sp 8
  sp_rates[7] -= (fwd_rates[740] - rev_rates[714]);

  //rxn 741
  //sp 97
  sp_rates[96] -= (fwd_rates[741] - rev_rates[715]);
  //sp 98
  sp_rates[97] += (fwd_rates[741] - rev_rates[715]);
  //sp 103
  sp_rates[102] += (fwd_rates[741] - rev_rates[715]);
  //sp 104
  sp_rates[103] -= (fwd_rates[741] - rev_rates[715]);

  //rxn 742
  //sp 97
  sp_rates[96] -= (fwd_rates[742] - rev_rates[716]);
  //sp 106
  sp_rates[105] += (fwd_rates[742] - rev_rates[716]);
  //sp 98
  sp_rates[97] += (fwd_rates[742] - rev_rates[716]);
  //sp 104
  sp_rates[103] -= (fwd_rates[742] - rev_rates[716]);

  //rxn 743
  //sp 97
  sp_rates[96] += (fwd_rates[743] - rev_rates[717]);
  //sp 103
  sp_rates[102] += (fwd_rates[743] - rev_rates[717]);
  //sp 96
  sp_rates[95] -= (fwd_rates[743] - rev_rates[717]);
  //sp 104
  sp_rates[103] -= (fwd_rates[743] - rev_rates[717]);

  //rxn 744
  //sp 97
  sp_rates[96] += (fwd_rates[744] - rev_rates[718]) * pres_mod[81];
  //sp 96
  sp_rates[95] += (fwd_rates[744] - rev_rates[718]) * pres_mod[81];
  //sp 104
  sp_rates[103] -= (fwd_rates[744] - rev_rates[718]) * pres_mod[81];

  //rxn 745
  //sp 105
  sp_rates[104] += (fwd_rates[745] - rev_rates[719]);
  //sp 102
  sp_rates[101] += (fwd_rates[745] - rev_rates[719]);
  //sp 103
  sp_rates[102] -= (fwd_rates[745] - rev_rates[719]);
  //sp 104
  sp_rates[103] -= (fwd_rates[745] - rev_rates[719]);

  //rxn 746
  //sp 98
  sp_rates[97] += 2.0 * (fwd_rates[746] - rev_rates[720]);
  //sp 2
  (*dy_N) += (fwd_rates[746] - rev_rates[720]);
  //sp 104
  sp_rates[103] -= 2.0 * (fwd_rates[746] - rev_rates[720]);

  //rxn 747
  //sp 102
  sp_rates[101] += (fwd_rates[747] - rev_rates[721]);
  //sp 103
  sp_rates[102] -= 2.0 * (fwd_rates[747] - rev_rates[721]);
  //sp 104
  sp_rates[103] += (fwd_rates[747] - rev_rates[721]);

  //rxn 748
  //sp 97
  sp_rates[96] -= 2.0 * (fwd_rates[748] - rev_rates[722]);
  //sp 106
  sp_rates[105] += (fwd_rates[748] - rev_rates[722]);
  //sp 6
  sp_rates[5] += (fwd_rates[748] - rev_rates[722]);

  //rxn 749
  //sp 106
  sp_rates[105] -= (fwd_rates[749] - rev_rates[723]);
  //sp 3
  sp_rates[2] += (fwd_rates[749] - rev_rates[723]);
  //sp 102
  sp_rates[101] += (fwd_rates[749] - rev_rates[723]);

  //rxn 750
  //sp 106
  sp_rates[105] -= (fwd_rates[750] - rev_rates[724]);
  //sp 3
  sp_rates[2] -= (fwd_rates[750] - rev_rates[724]);
  //sp 102
  sp_rates[101] += (fwd_rates[750] - rev_rates[724]);
  //sp 6
  sp_rates[5] += (fwd_rates[750] - rev_rates[724]);

  //rxn 751
  //sp 106
  sp_rates[105] -= (fwd_rates[751] - rev_rates[725]);
  //sp 103
  sp_rates[102] += (fwd_rates[751] - rev_rates[725]);

  //rxn 752
  //sp 106
  sp_rates[105] -= (fwd_rates[752] - rev_rates[726]);
  //sp 4
  sp_rates[3] -= (fwd_rates[752] - rev_rates[726]);
  //sp 5
  sp_rates[4] += (fwd_rates[752] - rev_rates[726]);
  //sp 102
  sp_rates[101] += (fwd_rates[752] - rev_rates[726]);

  //rxn 753
  //sp 97
  sp_rates[96] += (fwd_rates[753] - rev_rates[727]);
  //sp 106
  sp_rates[105] -= (fwd_rates[753] - rev_rates[727]);
  //sp 4
  sp_rates[3] -= (fwd_rates[753] - rev_rates[727]);
  //sp 93
  sp_rates[92] += (fwd_rates[753] - rev_rates[727]);

  //rxn 754
  //sp 106
  sp_rates[105] -= (fwd_rates[754] - rev_rates[728]);
  //sp 5
  sp_rates[4] -= (fwd_rates[754] - rev_rates[728]);
  //sp 102
  sp_rates[101] += (fwd_rates[754] - rev_rates[728]);
  //sp 10
  sp_rates[9] += (fwd_rates[754] - rev_rates[728]);

  //rxn 755
  //sp 97
  sp_rates[96] += fwd_rates[755];
  //sp 3
  sp_rates[2] += fwd_rates[755];
  //sp 5
  sp_rates[4] -= fwd_rates[755];
  //sp 106
  sp_rates[105] -= fwd_rates[755];
  //sp 93
  sp_rates[92] += fwd_rates[755];

  //rxn 756
  //sp 97
  sp_rates[96] += fwd_rates[756];
  //sp 5
  sp_rates[4] += fwd_rates[756];
  //sp 8
  sp_rates[7] -= fwd_rates[756];
  //sp 106
  sp_rates[105] -= fwd_rates[756];
  //sp 93
  sp_rates[92] += fwd_rates[756];

  //rxn 757
  //sp 9
  sp_rates[8] += (fwd_rates[757] - rev_rates[729]);
  //sp 106
  sp_rates[105] -= (fwd_rates[757] - rev_rates[729]);
  //sp 102
  sp_rates[101] += (fwd_rates[757] - rev_rates[729]);
  //sp 8
  sp_rates[7] -= (fwd_rates[757] - rev_rates[729]);

  //rxn 758
  //sp 97
  sp_rates[96] += (fwd_rates[758] - rev_rates[730]);
  //sp 106
  sp_rates[105] -= (fwd_rates[758] - rev_rates[730]);
  //sp 7
  sp_rates[6] -= (fwd_rates[758] - rev_rates[730]);
  //sp 95
  sp_rates[94] += (fwd_rates[758] - rev_rates[730]);

  //rxn 759
  //sp 97
  sp_rates[96] -= (fwd_rates[759] - rev_rates[731]);
  //sp 106
  sp_rates[105] -= (fwd_rates[759] - rev_rates[731]);
  //sp 102
  sp_rates[101] += (fwd_rates[759] - rev_rates[731]);
  //sp 98
  sp_rates[97] += (fwd_rates[759] - rev_rates[731]);

  //rxn 760
  //sp 107
  sp_rates[106] = (fwd_rates[760] - rev_rates[732]) * pres_mod[82];
  //sp 101
  sp_rates[100] -= (fwd_rates[760] - rev_rates[732]) * pres_mod[82];

  //rxn 761
  //sp 107
  sp_rates[106] -= (fwd_rates[761] - rev_rates[733]) * pres_mod[83];
  //sp 3
  sp_rates[2] += (fwd_rates[761] - rev_rates[733]) * pres_mod[83];
  //sp 99
  sp_rates[98] += (fwd_rates[761] - rev_rates[733]) * pres_mod[83];

  //rxn 762
  //sp 97
  sp_rates[96] += (fwd_rates[762] - rev_rates[734]);
  //sp 3
  sp_rates[2] -= (fwd_rates[762] - rev_rates[734]);
  //sp 5
  sp_rates[4] += (fwd_rates[762] - rev_rates[734]);
  //sp 107
  sp_rates[106] -= (fwd_rates[762] - rev_rates[734]);

  //rxn 763
  //sp 3
  sp_rates[2] -= (fwd_rates[763] - rev_rates[735]);
  //sp 6
  sp_rates[5] += (fwd_rates[763] - rev_rates[735]);
  //sp 107
  sp_rates[106] -= (fwd_rates[763] - rev_rates[735]);
  //sp 99
  sp_rates[98] += (fwd_rates[763] - rev_rates[735]);

  //rxn 764
  //sp 107
  sp_rates[106] -= (fwd_rates[764] - rev_rates[736]);
  //sp 4
  sp_rates[3] -= (fwd_rates[764] - rev_rates[736]);
  //sp 5
  sp_rates[4] += (fwd_rates[764] - rev_rates[736]);
  //sp 99
  sp_rates[98] += (fwd_rates[764] - rev_rates[736]);

  //rxn 765
  //sp 107
  sp_rates[106] -= (fwd_rates[765] - rev_rates[737]);
  //sp 4
  sp_rates[3] -= (fwd_rates[765] - rev_rates[737]);
  //sp 5
  sp_rates[4] += (fwd_rates[765] - rev_rates[737]);
  //sp 99
  sp_rates[98] += (fwd_rates[765] - rev_rates[737]);

  //rxn 766
  //sp 10
  sp_rates[9] += (fwd_rates[766] - rev_rates[738]);
  //sp 107
  sp_rates[106] -= (fwd_rates[766] - rev_rates[738]);
  //sp 5
  sp_rates[4] -= (fwd_rates[766] - rev_rates[738]);
  //sp 99
  sp_rates[98] += (fwd_rates[766] - rev_rates[738]);

  //rxn 767
  //sp 9
  sp_rates[8] += (fwd_rates[767] - rev_rates[739]);
  //sp 107
  sp_rates[106] -= (fwd_rates[767] - rev_rates[739]);
  //sp 99
  sp_rates[98] += (fwd_rates[767] - rev_rates[739]);
  //sp 8
  sp_rates[7] -= (fwd_rates[767] - rev_rates[739]);

  //rxn 768
  //sp 107
  sp_rates[106] -= (fwd_rates[768] - rev_rates[740]);
  //sp 8
  sp_rates[7] += (fwd_rates[768] - rev_rates[740]);
  //sp 7
  sp_rates[6] -= (fwd_rates[768] - rev_rates[740]);
  //sp 99
  sp_rates[98] += (fwd_rates[768] - rev_rates[740]);

  //rxn 769
  //sp 97
  sp_rates[96] -= (fwd_rates[769] - rev_rates[741]);
  //sp 107
  sp_rates[106] -= (fwd_rates[769] - rev_rates[741]);
  //sp 5
  sp_rates[4] += (fwd_rates[769] - rev_rates[741]);
  //sp 104
  sp_rates[103] += (fwd_rates[769] - rev_rates[741]);

  //rxn 770
  //sp 97
  sp_rates[96] -= (fwd_rates[770] - rev_rates[742]);
  //sp 106
  sp_rates[105] += (fwd_rates[770] - rev_rates[742]);
  //sp 107
  sp_rates[106] -= (fwd_rates[770] - rev_rates[742]);
  //sp 10
  sp_rates[9] += (fwd_rates[770] - rev_rates[742]);

  //rxn 771
  //sp 97
  sp_rates[96] -= (fwd_rates[771] - rev_rates[743]);
  //sp 98
  sp_rates[97] += (fwd_rates[771] - rev_rates[743]);
  //sp 107
  sp_rates[106] -= (fwd_rates[771] - rev_rates[743]);
  //sp 99
  sp_rates[98] += (fwd_rates[771] - rev_rates[743]);

  //rxn 772
  //sp 107
  sp_rates[106] -= (fwd_rates[772] - rev_rates[744]);
  //sp 100
  sp_rates[99] += (fwd_rates[772] - rev_rates[744]);
  //sp 95
  sp_rates[94] -= (fwd_rates[772] - rev_rates[744]);
  //sp 99
  sp_rates[98] += (fwd_rates[772] - rev_rates[744]);

  //rxn 773
  //sp 97
  sp_rates[96] += (fwd_rates[773] - rev_rates[745]) * pres_mod[84];
  //sp 108
  sp_rates[107] = -(fwd_rates[773] - rev_rates[745]) * pres_mod[84];
  //sp 5
  sp_rates[4] += (fwd_rates[773] - rev_rates[745]) * pres_mod[84];

  //rxn 774
  //sp 3
  sp_rates[2] -= (fwd_rates[774] - rev_rates[746]);
  //sp 108
  sp_rates[107] -= (fwd_rates[774] - rev_rates[746]);
  //sp 6
  sp_rates[5] += (fwd_rates[774] - rev_rates[746]);
  //sp 107
  sp_rates[106] += (fwd_rates[774] - rev_rates[746]);

  //rxn 775
  //sp 3
  sp_rates[2] -= (fwd_rates[775] - rev_rates[747]);
  //sp 108
  sp_rates[107] -= (fwd_rates[775] - rev_rates[747]);
  //sp 101
  sp_rates[100] += (fwd_rates[775] - rev_rates[747]);
  //sp 6
  sp_rates[5] += (fwd_rates[775] - rev_rates[747]);

  //rxn 776
  //sp 4
  sp_rates[3] -= (fwd_rates[776] - rev_rates[748]);
  //sp 107
  sp_rates[106] += (fwd_rates[776] - rev_rates[748]);
  //sp 108
  sp_rates[107] -= (fwd_rates[776] - rev_rates[748]);
  //sp 5
  sp_rates[4] += (fwd_rates[776] - rev_rates[748]);

  //rxn 777
  //sp 4
  sp_rates[3] -= (fwd_rates[777] - rev_rates[749]);
  //sp 108
  sp_rates[107] -= (fwd_rates[777] - rev_rates[749]);
  //sp 101
  sp_rates[100] += (fwd_rates[777] - rev_rates[749]);
  //sp 5
  sp_rates[4] += (fwd_rates[777] - rev_rates[749]);

  //rxn 778
  //sp 10
  sp_rates[9] += (fwd_rates[778] - rev_rates[750]);
  //sp 107
  sp_rates[106] += (fwd_rates[778] - rev_rates[750]);
  //sp 108
  sp_rates[107] -= (fwd_rates[778] - rev_rates[750]);
  //sp 5
  sp_rates[4] -= (fwd_rates[778] - rev_rates[750]);

  //rxn 779
  //sp 101
  sp_rates[100] += (fwd_rates[779] - rev_rates[751]);
  //sp 10
  sp_rates[9] += (fwd_rates[779] - rev_rates[751]);
  //sp 108
  sp_rates[107] -= (fwd_rates[779] - rev_rates[751]);
  //sp 5
  sp_rates[4] -= (fwd_rates[779] - rev_rates[751]);

  //rxn 780
  //sp 97
  sp_rates[96] -= (fwd_rates[780] - rev_rates[752]);
  //sp 98
  sp_rates[97] += (fwd_rates[780] - rev_rates[752]);
  //sp 107
  sp_rates[106] += (fwd_rates[780] - rev_rates[752]);
  //sp 108
  sp_rates[107] -= (fwd_rates[780] - rev_rates[752]);

  //rxn 781
  //sp 97
  sp_rates[96] -= (fwd_rates[781] - rev_rates[753]);
  //sp 98
  sp_rates[97] += (fwd_rates[781] - rev_rates[753]);
  //sp 108
  sp_rates[107] -= (fwd_rates[781] - rev_rates[753]);
  //sp 101
  sp_rates[100] += (fwd_rates[781] - rev_rates[753]);

  //rxn 782
  //sp 97
  sp_rates[96] += (fwd_rates[782] - rev_rates[754]);
  //sp 107
  sp_rates[106] += (fwd_rates[782] - rev_rates[754]);
  //sp 108
  sp_rates[107] -= (fwd_rates[782] - rev_rates[754]);
  //sp 96
  sp_rates[95] -= (fwd_rates[782] - rev_rates[754]);

  //rxn 783
  //sp 97
  sp_rates[96] += (fwd_rates[783] - rev_rates[755]);
  //sp 108
  sp_rates[107] -= (fwd_rates[783] - rev_rates[755]);
  //sp 101
  sp_rates[100] += (fwd_rates[783] - rev_rates[755]);
  //sp 96
  sp_rates[95] -= (fwd_rates[783] - rev_rates[755]);

  //rxn 784
  //sp 9
  sp_rates[8] += (fwd_rates[784] - rev_rates[756]);
  //sp 107
  sp_rates[106] += (fwd_rates[784] - rev_rates[756]);
  //sp 108
  sp_rates[107] -= (fwd_rates[784] - rev_rates[756]);
  //sp 8
  sp_rates[7] -= (fwd_rates[784] - rev_rates[756]);

  //rxn 785
  //sp 9
  sp_rates[8] += (fwd_rates[785] - rev_rates[757]);
  //sp 108
  sp_rates[107] -= (fwd_rates[785] - rev_rates[757]);
  //sp 101
  sp_rates[100] += (fwd_rates[785] - rev_rates[757]);
  //sp 8
  sp_rates[7] -= (fwd_rates[785] - rev_rates[757]);

  //rxn 786
  //sp 99
  sp_rates[98] -= (fwd_rates[786] - rev_rates[758]);
  //sp 108
  sp_rates[107] += (fwd_rates[786] - rev_rates[758]);
  //sp 93
  sp_rates[92] += (fwd_rates[786] - rev_rates[758]);
  //sp 107
  sp_rates[106] -= (fwd_rates[786] - rev_rates[758]);

  //rxn 787
  //sp 97
  sp_rates[96] -= (fwd_rates[787] - rev_rates[759]) * pres_mod[85];
  //sp 98
  sp_rates[97] += (fwd_rates[787] - rev_rates[759]) * pres_mod[85];
  //sp 3
  sp_rates[2] -= (fwd_rates[787] - rev_rates[759]) * pres_mod[85];

  //rxn 788
  //sp 98
  sp_rates[97] -= (fwd_rates[788] - rev_rates[760]) * pres_mod[86];
  //sp 6
  sp_rates[5] += (fwd_rates[788] - rev_rates[760]) * pres_mod[86];
  //sp 96
  sp_rates[95] += (fwd_rates[788] - rev_rates[760]) * pres_mod[86];

  //rxn 789
  //sp 97
  sp_rates[96] += (fwd_rates[789] - rev_rates[761]);
  //sp 98
  sp_rates[97] -= (fwd_rates[789] - rev_rates[761]);
  //sp 3
  sp_rates[2] -= (fwd_rates[789] - rev_rates[761]);
  //sp 6
  sp_rates[5] += (fwd_rates[789] - rev_rates[761]);

  //rxn 790
  //sp 97
  sp_rates[96] += (fwd_rates[790] - rev_rates[762]);
  //sp 98
  sp_rates[97] -= (fwd_rates[790] - rev_rates[762]);
  //sp 4
  sp_rates[3] -= (fwd_rates[790] - rev_rates[762]);
  //sp 5
  sp_rates[4] += (fwd_rates[790] - rev_rates[762]);

  //rxn 791
  //sp 97
  sp_rates[96] += (fwd_rates[791] - rev_rates[763]);
  //sp 98
  sp_rates[97] -= (fwd_rates[791] - rev_rates[763]);
  //sp 5
  sp_rates[4] -= (fwd_rates[791] - rev_rates[763]);
  //sp 10
  sp_rates[9] += (fwd_rates[791] - rev_rates[763]);

  //rxn 792
  //sp 9
  sp_rates[8] += (fwd_rates[792] - rev_rates[764]);
  //sp 98
  sp_rates[97] -= (fwd_rates[792] - rev_rates[764]);
  //sp 97
  sp_rates[96] += (fwd_rates[792] - rev_rates[764]);
  //sp 8
  sp_rates[7] -= (fwd_rates[792] - rev_rates[764]);

  //rxn 793
  //sp 97
  sp_rates[96] -= (fwd_rates[793] - rev_rates[765]);
  //sp 98
  sp_rates[97] -= (fwd_rates[793] - rev_rates[765]);
  //sp 6
  sp_rates[5] += (fwd_rates[793] - rev_rates[765]);
  //sp 104
  sp_rates[103] += (fwd_rates[793] - rev_rates[765]);

  //rxn 794
  //sp 3
  sp_rates[2] += (fwd_rates[794] - rev_rates[766]);
  //sp 2
  (*dy_N) += (fwd_rates[794] - rev_rates[766]);
  //sp 102
  sp_rates[101] -= (fwd_rates[794] - rev_rates[766]);

  //rxn 795
  //sp 3
  sp_rates[2] -= (fwd_rates[795] - rev_rates[767]);
  //sp 6
  sp_rates[5] += (fwd_rates[795] - rev_rates[767]);
  //sp 2
  (*dy_N) += (fwd_rates[795] - rev_rates[767]);
  //sp 102
  sp_rates[101] -= (fwd_rates[795] - rev_rates[767]);

  //rxn 796
  //sp 3
  sp_rates[2] += (fwd_rates[796] - rev_rates[768]);
  //sp 4
  sp_rates[3] -= (fwd_rates[796] - rev_rates[768]);
  //sp 94
  sp_rates[93] += (fwd_rates[796] - rev_rates[768]);
  //sp 102
  sp_rates[101] -= (fwd_rates[796] - rev_rates[768]);

  //rxn 797
  //sp 4
  sp_rates[3] -= (fwd_rates[797] - rev_rates[769]);
  //sp 93
  sp_rates[92] += (fwd_rates[797] - rev_rates[769]);
  //sp 102
  sp_rates[101] -= (fwd_rates[797] - rev_rates[769]);
  //sp 96
  sp_rates[95] += (fwd_rates[797] - rev_rates[769]);

  //rxn 798
  //sp 5
  sp_rates[4] += (fwd_rates[798] - rev_rates[770]);
  //sp 4
  sp_rates[3] -= (fwd_rates[798] - rev_rates[770]);
  //sp 2
  (*dy_N) += (fwd_rates[798] - rev_rates[770]);
  //sp 102
  sp_rates[101] -= (fwd_rates[798] - rev_rates[770]);

  //rxn 799
  //sp 2
  (*dy_N) += (fwd_rates[799] - rev_rates[771]);
  //sp 10
  sp_rates[9] += (fwd_rates[799] - rev_rates[771]);
  //sp 5
  sp_rates[4] -= (fwd_rates[799] - rev_rates[771]);
  //sp 102
  sp_rates[101] -= (fwd_rates[799] - rev_rates[771]);

  //rxn 800
  //sp 2
  (*dy_N) += (fwd_rates[800] - rev_rates[772]);
  //sp 102
  sp_rates[101] -= (fwd_rates[800] - rev_rates[772]);
  //sp 7
  sp_rates[6] -= (fwd_rates[800] - rev_rates[772]);
  //sp 8
  sp_rates[7] += (fwd_rates[800] - rev_rates[772]);

  //rxn 801
  //sp 2
  (*dy_N) += (fwd_rates[801] - rev_rates[773]);
  //sp 99
  sp_rates[98] += (fwd_rates[801] - rev_rates[773]);
  //sp 93
  sp_rates[92] -= (fwd_rates[801] - rev_rates[773]);
  //sp 102
  sp_rates[101] -= (fwd_rates[801] - rev_rates[773]);

  //rxn 802
  //sp 97
  sp_rates[96] += (fwd_rates[802] - rev_rates[774]);
  //sp 2
  (*dy_N) += (fwd_rates[802] - rev_rates[774]);
  //sp 102
  sp_rates[101] -= (fwd_rates[802] - rev_rates[774]);
  //sp 96
  sp_rates[95] -= (fwd_rates[802] - rev_rates[774]);

  //rxn 803
  //sp 97
  sp_rates[96] -= (fwd_rates[803] - rev_rates[775]);
  //sp 98
  sp_rates[97] += (fwd_rates[803] - rev_rates[775]);
  //sp 2
  (*dy_N) += (fwd_rates[803] - rev_rates[775]);
  //sp 102
  sp_rates[101] -= (fwd_rates[803] - rev_rates[775]);

  //rxn 804
  //sp 2
  (*dy_N) += (fwd_rates[804] - rev_rates[776]);
  //sp 102
  sp_rates[101] -= 2.0 * (fwd_rates[804] - rev_rates[776]);
  //sp 103
  sp_rates[102] += (fwd_rates[804] - rev_rates[776]);

  //rxn 805
  //sp 3
  sp_rates[2] -= (fwd_rates[805] - rev_rates[777]) * pres_mod[87];
  //sp 93
  sp_rates[92] -= (fwd_rates[805] - rev_rates[777]) * pres_mod[87];
  //sp 99
  sp_rates[98] += (fwd_rates[805] - rev_rates[777]) * pres_mod[87];

  //rxn 806
  //sp 99
  sp_rates[98] -= (fwd_rates[806] - rev_rates[778]);
  //sp 4
  sp_rates[3] -= (fwd_rates[806] - rev_rates[778]);
  //sp 93
  sp_rates[92] += (fwd_rates[806] - rev_rates[778]);
  //sp 5
  sp_rates[4] += (fwd_rates[806] - rev_rates[778]);

  //rxn 807
  //sp 3
  sp_rates[2] -= (fwd_rates[807] - rev_rates[779]);
  //sp 93
  sp_rates[92] += (fwd_rates[807] - rev_rates[779]);
  //sp 6
  sp_rates[5] += (fwd_rates[807] - rev_rates[779]);
  //sp 99
  sp_rates[98] -= (fwd_rates[807] - rev_rates[779]);

  //rxn 808
  //sp 10
  sp_rates[9] += (fwd_rates[808] - rev_rates[780]);
  //sp 99
  sp_rates[98] -= (fwd_rates[808] - rev_rates[780]);
  //sp 5
  sp_rates[4] -= (fwd_rates[808] - rev_rates[780]);
  //sp 93
  sp_rates[92] += (fwd_rates[808] - rev_rates[780]);

  //rxn 809
  //sp 99
  sp_rates[98] -= (fwd_rates[809] - rev_rates[781]);
  //sp 93
  sp_rates[92] += (fwd_rates[809] - rev_rates[781]);
  //sp 7
  sp_rates[6] -= (fwd_rates[809] - rev_rates[781]);
  //sp 8
  sp_rates[7] += (fwd_rates[809] - rev_rates[781]);

  //rxn 810
  //sp 99
  sp_rates[98] -= (fwd_rates[810] - rev_rates[782]);
  //sp 92
  sp_rates[91] -= (fwd_rates[810] - rev_rates[782]);
  //sp 93
  sp_rates[92] += (fwd_rates[810] - rev_rates[782]);
  //sp 96
  sp_rates[95] += (fwd_rates[810] - rev_rates[782]);

  //rxn 811
  //sp 3
  sp_rates[2] += (fwd_rates[811] - rev_rates[783]);
  //sp 99
  sp_rates[98] -= (fwd_rates[811] - rev_rates[783]);
  //sp 92
  sp_rates[91] -= (fwd_rates[811] - rev_rates[783]);
  //sp 94
  sp_rates[93] += (fwd_rates[811] - rev_rates[783]);

  //rxn 812
  //sp 99
  sp_rates[98] -= (fwd_rates[812] - rev_rates[784]);
  //sp 5
  sp_rates[4] += (fwd_rates[812] - rev_rates[784]);
  //sp 93
  sp_rates[92] -= (fwd_rates[812] - rev_rates[784]);
  //sp 94
  sp_rates[93] += (fwd_rates[812] - rev_rates[784]);

  //rxn 813
  //sp 99
  sp_rates[98] -= (fwd_rates[813] - rev_rates[785]);
  //sp 100
  sp_rates[99] += (fwd_rates[813] - rev_rates[785]);
  //sp 93
  sp_rates[92] += (fwd_rates[813] - rev_rates[785]);
  //sp 95
  sp_rates[94] -= (fwd_rates[813] - rev_rates[785]);

  //rxn 814
  //sp 10
  sp_rates[9] += (fwd_rates[814] - rev_rates[786]);
  //sp 99
  sp_rates[98] -= 2.0 * (fwd_rates[814] - rev_rates[786]);
  //sp 94
  sp_rates[93] += (fwd_rates[814] - rev_rates[786]);

  //rxn 815
  //sp 97
  sp_rates[96] += (fwd_rates[815] - rev_rates[787]);
  //sp 99
  sp_rates[98] -= (fwd_rates[815] - rev_rates[787]);
  //sp 93
  sp_rates[92] += (fwd_rates[815] - rev_rates[787]);
  //sp 96
  sp_rates[95] -= (fwd_rates[815] - rev_rates[787]);

  //rxn 816
  //sp 122
  sp_rates[121] = -(fwd_rates[816] - rev_rates[788]);
  //sp 123
  sp_rates[122] += (fwd_rates[816] - rev_rates[788]);
  //sp 12
  sp_rates[11] += (fwd_rates[816] - rev_rates[788]);
  //sp 5
  sp_rates[4] -= (fwd_rates[816] - rev_rates[788]);

  //rxn 817
  //sp 123
  sp_rates[122] -= (fwd_rates[817] - rev_rates[789]) * pres_mod[88];
  //sp 3
  sp_rates[2] += (fwd_rates[817] - rev_rates[789]) * pres_mod[88];
  //sp 93
  sp_rates[92] += (fwd_rates[817] - rev_rates[789]) * pres_mod[88];

  //rxn 818
  //sp 123
  sp_rates[122] -= (fwd_rates[818] - rev_rates[790]);
  //sp 99
  sp_rates[98] += (fwd_rates[818] - rev_rates[790]);

  //rxn 819
  //sp 3
  sp_rates[2] -= (fwd_rates[819] - rev_rates[791]);
  //sp 5
  sp_rates[4] += (fwd_rates[819] - rev_rates[791]);
  //sp 123
  sp_rates[122] -= (fwd_rates[819] - rev_rates[791]);
  //sp 96
  sp_rates[95] += (fwd_rates[819] - rev_rates[791]);

  //rxn 820
  //sp 123
  sp_rates[122] -= (fwd_rates[820] - rev_rates[792]);
  //sp 4
  sp_rates[3] -= (fwd_rates[820] - rev_rates[792]);
  //sp 93
  sp_rates[92] += (fwd_rates[820] - rev_rates[792]);
  //sp 5
  sp_rates[4] += (fwd_rates[820] - rev_rates[792]);

  //rxn 821
  //sp 100
  sp_rates[99] += (fwd_rates[821] - rev_rates[793]);
  //sp 123
  sp_rates[122] -= (fwd_rates[821] - rev_rates[793]);
  //sp 3
  sp_rates[2] += (fwd_rates[821] - rev_rates[793]);
  //sp 5
  sp_rates[4] -= (fwd_rates[821] - rev_rates[793]);

  //rxn 822
  //sp 95
  sp_rates[94] += (fwd_rates[822] - rev_rates[794]);
  //sp 123
  sp_rates[122] -= (fwd_rates[822] - rev_rates[794]);
  //sp 5
  sp_rates[4] += (fwd_rates[822] - rev_rates[794]);
  //sp 7
  sp_rates[6] -= (fwd_rates[822] - rev_rates[794]);

  //rxn 823
  //sp 100
  sp_rates[99] += (fwd_rates[823] - rev_rates[795]) * pres_mod[89];
  //sp 93
  sp_rates[92] -= (fwd_rates[823] - rev_rates[795]) * pres_mod[89];
  //sp 5
  sp_rates[4] -= (fwd_rates[823] - rev_rates[795]) * pres_mod[89];

  //rxn 824
  //sp 93
  sp_rates[92] -= (fwd_rates[824] - rev_rates[796]);
  //sp 5
  sp_rates[4] += (fwd_rates[824] - rev_rates[796]);
  //sp 95
  sp_rates[94] += (fwd_rates[824] - rev_rates[796]);
  //sp 8
  sp_rates[7] -= (fwd_rates[824] - rev_rates[796]);

  //rxn 825
  //sp 3
  sp_rates[2] += (fwd_rates[825] - rev_rates[797]);
  //sp 100
  sp_rates[99] += (fwd_rates[825] - rev_rates[797]);
  //sp 6
  sp_rates[5] -= (fwd_rates[825] - rev_rates[797]);
  //sp 95
  sp_rates[94] -= (fwd_rates[825] - rev_rates[797]);

  //rxn 826
  //sp 3
  sp_rates[2] += (fwd_rates[826] - rev_rates[798]);
  //sp 109
  sp_rates[108] = (fwd_rates[826] - rev_rates[798]);
  //sp 6
  sp_rates[5] -= (fwd_rates[826] - rev_rates[798]);
  //sp 95
  sp_rates[94] -= (fwd_rates[826] - rev_rates[798]);

  //rxn 827
  //sp 3
  sp_rates[2] -= (fwd_rates[827] - rev_rates[799]);
  //sp 93
  sp_rates[92] += (fwd_rates[827] - rev_rates[799]);
  //sp 5
  sp_rates[4] += (fwd_rates[827] - rev_rates[799]);
  //sp 95
  sp_rates[94] -= (fwd_rates[827] - rev_rates[799]);

  //rxn 828
  //sp 4
  sp_rates[3] -= (fwd_rates[828] - rev_rates[800]);
  //sp 93
  sp_rates[92] += (fwd_rates[828] - rev_rates[800]);
  //sp 95
  sp_rates[94] -= (fwd_rates[828] - rev_rates[800]);
  //sp 7
  sp_rates[6] += (fwd_rates[828] - rev_rates[800]);

  //rxn 829
  //sp 4
  sp_rates[3] -= (fwd_rates[829] - rev_rates[801]) * pres_mod[90];
  //sp 93
  sp_rates[92] -= (fwd_rates[829] - rev_rates[801]) * pres_mod[90];
  //sp 95
  sp_rates[94] += (fwd_rates[829] - rev_rates[801]) * pres_mod[90];

  //rxn 830
  //sp 111
  sp_rates[110] = (fwd_rates[830] - rev_rates[802]) * pres_mod[91];
  //sp 5
  sp_rates[4] -= (fwd_rates[830] - rev_rates[802]) * pres_mod[91];
  //sp 95
  sp_rates[94] -= (fwd_rates[830] - rev_rates[802]) * pres_mod[91];

  //rxn 831
  //sp 93
  sp_rates[92] += 2.0 * (fwd_rates[831] - rev_rates[803]);
  //sp 95
  sp_rates[94] -= 2.0 * (fwd_rates[831] - rev_rates[803]);
  //sp 7
  sp_rates[6] += (fwd_rates[831] - rev_rates[803]);

  //rxn 832
  //sp 100
  sp_rates[99] += (fwd_rates[832] - rev_rates[804]);
  //sp 7
  sp_rates[6] += (fwd_rates[832] - rev_rates[804]);
  //sp 95
  sp_rates[94] -= (fwd_rates[832] - rev_rates[804]);
  //sp 8
  sp_rates[7] -= (fwd_rates[832] - rev_rates[804]);

  //rxn 833
  //sp 109
  sp_rates[108] += (fwd_rates[833] - rev_rates[805]);
  //sp 7
  sp_rates[6] += (fwd_rates[833] - rev_rates[805]);
  //sp 95
  sp_rates[94] -= (fwd_rates[833] - rev_rates[805]);
  //sp 8
  sp_rates[7] -= (fwd_rates[833] - rev_rates[805]);

  //rxn 834
  //sp 4
  sp_rates[3] += (fwd_rates[834] - rev_rates[806]) * pres_mod[92];
  //sp 2
  (*dy_N) += (fwd_rates[834] - rev_rates[806]) * pres_mod[92];
  //sp 94
  sp_rates[93] -= (fwd_rates[834] - rev_rates[806]) * pres_mod[92];

  //rxn 835
  //sp 3
  sp_rates[2] -= (fwd_rates[835] - rev_rates[807]);
  //sp 5
  sp_rates[4] += (fwd_rates[835] - rev_rates[807]);
  //sp 2
  (*dy_N) += (fwd_rates[835] - rev_rates[807]);
  //sp 94
  sp_rates[93] -= (fwd_rates[835] - rev_rates[807]);

  //rxn 836
  //sp 3
  sp_rates[2] -= (fwd_rates[836] - rev_rates[808]);
  //sp 5
  sp_rates[4] += (fwd_rates[836] - rev_rates[808]);
  //sp 2
  (*dy_N) += (fwd_rates[836] - rev_rates[808]);
  //sp 94
  sp_rates[93] -= (fwd_rates[836] - rev_rates[808]);

  //rxn 837
  //sp 4
  sp_rates[3] -= (fwd_rates[837] - rev_rates[809]);
  //sp 93
  sp_rates[92] += 2.0 * (fwd_rates[837] - rev_rates[809]);
  //sp 94
  sp_rates[93] -= (fwd_rates[837] - rev_rates[809]);

  //rxn 838
  //sp 4
  sp_rates[3] -= (fwd_rates[838] - rev_rates[810]);
  //sp 2
  (*dy_N) += (fwd_rates[838] - rev_rates[810]);
  //sp 94
  sp_rates[93] -= (fwd_rates[838] - rev_rates[810]);
  //sp 7
  sp_rates[6] += (fwd_rates[838] - rev_rates[810]);

  //rxn 839
  //sp 2
  (*dy_N) += (fwd_rates[839] - rev_rates[811]);
  //sp 5
  sp_rates[4] -= (fwd_rates[839] - rev_rates[811]);
  //sp 94
  sp_rates[93] -= (fwd_rates[839] - rev_rates[811]);
  //sp 8
  sp_rates[7] += (fwd_rates[839] - rev_rates[811]);

  //rxn 840
  //sp 2
  (*dy_N) += (fwd_rates[840] - rev_rates[812]);
  //sp 93
  sp_rates[92] -= (fwd_rates[840] - rev_rates[812]);
  //sp 94
  sp_rates[93] -= (fwd_rates[840] - rev_rates[812]);
  //sp 95
  sp_rates[94] += (fwd_rates[840] - rev_rates[812]);

  //rxn 841
  //sp 92
  sp_rates[91] -= (fwd_rates[841] - rev_rates[813]);
  //sp 2
  (*dy_N) += (fwd_rates[841] - rev_rates[813]);
  //sp 94
  sp_rates[93] -= (fwd_rates[841] - rev_rates[813]);
  //sp 93
  sp_rates[92] += (fwd_rates[841] - rev_rates[813]);

  //rxn 842
  //sp 3
  sp_rates[2] -= (fwd_rates[842] - rev_rates[814]);
  //sp 11
  sp_rates[10] += (fwd_rates[842] - rev_rates[814]);
  //sp 2
  (*dy_N) += (fwd_rates[842] - rev_rates[814]);
  //sp 94
  sp_rates[93] -= (fwd_rates[842] - rev_rates[814]);

  //rxn 843
  //sp 3
  sp_rates[2] += (fwd_rates[843] - rev_rates[815]) * pres_mod[93];
  //sp 101
  sp_rates[100] -= (fwd_rates[843] - rev_rates[815]) * pres_mod[93];
  //sp 99
  sp_rates[98] += (fwd_rates[843] - rev_rates[815]) * pres_mod[93];

  //rxn 844
  //sp 3
  sp_rates[2] -= (fwd_rates[844] - rev_rates[816]);
  //sp 101
  sp_rates[100] -= (fwd_rates[844] - rev_rates[816]);
  //sp 6
  sp_rates[5] += (fwd_rates[844] - rev_rates[816]);
  //sp 99
  sp_rates[98] += (fwd_rates[844] - rev_rates[816]);

  //rxn 845
  //sp 97
  sp_rates[96] += (fwd_rates[845] - rev_rates[817]);
  //sp 3
  sp_rates[2] -= (fwd_rates[845] - rev_rates[817]);
  //sp 101
  sp_rates[100] -= (fwd_rates[845] - rev_rates[817]);
  //sp 5
  sp_rates[4] += (fwd_rates[845] - rev_rates[817]);

  //rxn 846
  //sp 99
  sp_rates[98] += (fwd_rates[846] - rev_rates[818]);
  //sp 4
  sp_rates[3] -= (fwd_rates[846] - rev_rates[818]);
  //sp 101
  sp_rates[100] -= (fwd_rates[846] - rev_rates[818]);
  //sp 5
  sp_rates[4] += (fwd_rates[846] - rev_rates[818]);

  //rxn 847
  //sp 10
  sp_rates[9] += (fwd_rates[847] - rev_rates[819]);
  //sp 99
  sp_rates[98] += (fwd_rates[847] - rev_rates[819]);
  //sp 101
  sp_rates[100] -= (fwd_rates[847] - rev_rates[819]);
  //sp 5
  sp_rates[4] -= (fwd_rates[847] - rev_rates[819]);

  //rxn 848
  //sp 99
  sp_rates[98] += 2.0 * (fwd_rates[848] - rev_rates[820]);
  //sp 101
  sp_rates[100] -= (fwd_rates[848] - rev_rates[820]);
  //sp 93
  sp_rates[92] -= (fwd_rates[848] - rev_rates[820]);

  //rxn 849
  //sp 97
  sp_rates[96] -= (fwd_rates[849] - rev_rates[821]);
  //sp 98
  sp_rates[97] += (fwd_rates[849] - rev_rates[821]);
  //sp 99
  sp_rates[98] += (fwd_rates[849] - rev_rates[821]);
  //sp 101
  sp_rates[100] -= (fwd_rates[849] - rev_rates[821]);

  //rxn 850
  //sp 9
  sp_rates[8] += (fwd_rates[850] - rev_rates[822]);
  //sp 99
  sp_rates[98] += (fwd_rates[850] - rev_rates[822]);
  //sp 101
  sp_rates[100] -= (fwd_rates[850] - rev_rates[822]);
  //sp 8
  sp_rates[7] -= (fwd_rates[850] - rev_rates[822]);

  //rxn 851
  //sp 99
  sp_rates[98] += (fwd_rates[851] - rev_rates[823]);
  //sp 101
  sp_rates[100] -= (fwd_rates[851] - rev_rates[823]);
  //sp 7
  sp_rates[6] -= (fwd_rates[851] - rev_rates[823]);
  //sp 8
  sp_rates[7] += (fwd_rates[851] - rev_rates[823]);

  //rxn 852
  //sp 99
  sp_rates[98] += (fwd_rates[852] - rev_rates[824]);
  //sp 100
  sp_rates[99] += (fwd_rates[852] - rev_rates[824]);
  //sp 101
  sp_rates[100] -= (fwd_rates[852] - rev_rates[824]);
  //sp 95
  sp_rates[94] -= (fwd_rates[852] - rev_rates[824]);

  //rxn 853
  //sp 4
  sp_rates[3] -= (fwd_rates[853] - rev_rates[825]);
  //sp 100
  sp_rates[99] -= (fwd_rates[853] - rev_rates[825]);
  //sp 5
  sp_rates[4] += (fwd_rates[853] - rev_rates[825]);
  //sp 95
  sp_rates[94] += (fwd_rates[853] - rev_rates[825]);

  //rxn 854
  //sp 3
  sp_rates[2] -= (fwd_rates[854] - rev_rates[826]);
  //sp 100
  sp_rates[99] -= (fwd_rates[854] - rev_rates[826]);
  //sp 5
  sp_rates[4] += (fwd_rates[854] - rev_rates[826]);
  //sp 99
  sp_rates[98] += (fwd_rates[854] - rev_rates[826]);

  //rxn 855
  //sp 10
  sp_rates[9] += (fwd_rates[855] - rev_rates[827]);
  //sp 3
  sp_rates[2] -= (fwd_rates[855] - rev_rates[827]);
  //sp 100
  sp_rates[99] -= (fwd_rates[855] - rev_rates[827]);
  //sp 93
  sp_rates[92] += (fwd_rates[855] - rev_rates[827]);

  //rxn 856
  //sp 10
  sp_rates[9] += (fwd_rates[856] - rev_rates[828]);
  //sp 100
  sp_rates[99] -= (fwd_rates[856] - rev_rates[828]);
  //sp 5
  sp_rates[4] -= (fwd_rates[856] - rev_rates[828]);
  //sp 95
  sp_rates[94] += (fwd_rates[856] - rev_rates[828]);

  //rxn 857
  //sp 111
  sp_rates[110] += (fwd_rates[857] - rev_rates[829]);
  //sp 100
  sp_rates[99] -= (fwd_rates[857] - rev_rates[829]);
  //sp 93
  sp_rates[92] += (fwd_rates[857] - rev_rates[829]);
  //sp 95
  sp_rates[94] -= (fwd_rates[857] - rev_rates[829]);

  //rxn 858
  //sp 10
  sp_rates[9] += (fwd_rates[858] - rev_rates[830]);
  //sp 100
  sp_rates[99] -= 2.0 * (fwd_rates[858] - rev_rates[830]);
  //sp 93
  sp_rates[92] += (fwd_rates[858] - rev_rates[830]);
  //sp 95
  sp_rates[94] += (fwd_rates[858] - rev_rates[830]);

  //rxn 859
  //sp 100
  sp_rates[99] += (fwd_rates[859] - rev_rates[831]) * pres_mod[94];
  //sp 109
  sp_rates[108] -= (fwd_rates[859] - rev_rates[831]) * pres_mod[94];

  //rxn 860
  //sp 4
  sp_rates[3] -= (fwd_rates[860] - rev_rates[832]);
  //sp 109
  sp_rates[108] -= (fwd_rates[860] - rev_rates[832]);
  //sp 5
  sp_rates[4] += (fwd_rates[860] - rev_rates[832]);
  //sp 95
  sp_rates[94] += (fwd_rates[860] - rev_rates[832]);

  //rxn 861
  //sp 10
  sp_rates[9] += (fwd_rates[861] - rev_rates[833]);
  //sp 109
  sp_rates[108] -= (fwd_rates[861] - rev_rates[833]);
  //sp 5
  sp_rates[4] -= (fwd_rates[861] - rev_rates[833]);
  //sp 95
  sp_rates[94] += (fwd_rates[861] - rev_rates[833]);

  //rxn 862
  //sp 4
  sp_rates[3] -= (fwd_rates[862] - rev_rates[834]) * pres_mod[95];
  //sp 110
  sp_rates[109] = (fwd_rates[862] - rev_rates[834]) * pres_mod[95];
  //sp 95
  sp_rates[94] -= (fwd_rates[862] - rev_rates[834]) * pres_mod[95];

  //rxn 863
  //sp 93
  sp_rates[92] += (fwd_rates[863] - rev_rates[835]);
  //sp 110
  sp_rates[109] += (fwd_rates[863] - rev_rates[835]);
  //sp 95
  sp_rates[94] -= 2.0 * (fwd_rates[863] - rev_rates[835]);

  //rxn 864
  //sp 3
  sp_rates[2] -= (fwd_rates[864] - rev_rates[836]);
  //sp 5
  sp_rates[4] += (fwd_rates[864] - rev_rates[836]);
  //sp 110
  sp_rates[109] -= (fwd_rates[864] - rev_rates[836]);
  //sp 95
  sp_rates[94] += (fwd_rates[864] - rev_rates[836]);

  //rxn 865
  //sp 4
  sp_rates[3] -= (fwd_rates[865] - rev_rates[837]);
  //sp 110
  sp_rates[109] -= (fwd_rates[865] - rev_rates[837]);
  //sp 95
  sp_rates[94] += (fwd_rates[865] - rev_rates[837]);
  //sp 7
  sp_rates[6] += (fwd_rates[865] - rev_rates[837]);

  //rxn 866
  //sp 5
  sp_rates[4] -= (fwd_rates[866] - rev_rates[838]);
  //sp 110
  sp_rates[109] -= (fwd_rates[866] - rev_rates[838]);
  //sp 95
  sp_rates[94] += (fwd_rates[866] - rev_rates[838]);
  //sp 8
  sp_rates[7] += (fwd_rates[866] - rev_rates[838]);

  //rxn 867
  //sp 5
  sp_rates[4] += (fwd_rates[867] - rev_rates[839]);
  //sp 7
  sp_rates[6] += (fwd_rates[867] - rev_rates[839]);
  //sp 8
  sp_rates[7] -= (fwd_rates[867] - rev_rates[839]);
  //sp 110
  sp_rates[109] -= (fwd_rates[867] - rev_rates[839]);
  //sp 95
  sp_rates[94] += (fwd_rates[867] - rev_rates[839]);

  //rxn 868
  //sp 93
  sp_rates[92] += (fwd_rates[868] - rev_rates[840]) * pres_mod[96];
  //sp 110
  sp_rates[109] -= (fwd_rates[868] - rev_rates[840]) * pres_mod[96];
  //sp 7
  sp_rates[6] += (fwd_rates[868] - rev_rates[840]) * pres_mod[96];

  //rxn 869
  //sp 110
  sp_rates[109] -= 2.0 * (fwd_rates[869] - rev_rates[841]);
  //sp 95
  sp_rates[94] += 2.0 * (fwd_rates[869] - rev_rates[841]);
  //sp 7
  sp_rates[6] += (fwd_rates[869] - rev_rates[841]);

  //rxn 870
  //sp 7
  sp_rates[6] += (fwd_rates[870] - rev_rates[842]);
  //sp 110
  sp_rates[109] -= (fwd_rates[870] - rev_rates[842]);
  //sp 111
  sp_rates[110] += (fwd_rates[870] - rev_rates[842]);
  //sp 8
  sp_rates[7] -= (fwd_rates[870] - rev_rates[842]);

  //rxn 871
  //sp 93
  sp_rates[92] -= (fwd_rates[871] - rev_rates[843]) * pres_mod[97];
  //sp 111
  sp_rates[110] += (fwd_rates[871] - rev_rates[843]) * pres_mod[97];
  //sp 8
  sp_rates[7] -= (fwd_rates[871] - rev_rates[843]) * pres_mod[97];

  //rxn 872
  //sp 3
  sp_rates[2] -= (fwd_rates[872] - rev_rates[844]);
  //sp 110
  sp_rates[109] += (fwd_rates[872] - rev_rates[844]);
  //sp 6
  sp_rates[5] += (fwd_rates[872] - rev_rates[844]);
  //sp 111
  sp_rates[110] -= (fwd_rates[872] - rev_rates[844]);

  //rxn 873
  //sp 10
  sp_rates[9] += (fwd_rates[873] - rev_rates[845]);
  //sp 3
  sp_rates[2] -= (fwd_rates[873] - rev_rates[845]);
  //sp 111
  sp_rates[110] -= (fwd_rates[873] - rev_rates[845]);
  //sp 95
  sp_rates[94] += (fwd_rates[873] - rev_rates[845]);

  //rxn 874
  //sp 3
  sp_rates[2] -= (fwd_rates[874] - rev_rates[846]);
  //sp 100
  sp_rates[99] += (fwd_rates[874] - rev_rates[846]);
  //sp 5
  sp_rates[4] += (fwd_rates[874] - rev_rates[846]);
  //sp 111
  sp_rates[110] -= (fwd_rates[874] - rev_rates[846]);

  //rxn 875
  //sp 10
  sp_rates[9] += (fwd_rates[875] - rev_rates[847]);
  //sp 5
  sp_rates[4] -= (fwd_rates[875] - rev_rates[847]);
  //sp 110
  sp_rates[109] += (fwd_rates[875] - rev_rates[847]);
  //sp 111
  sp_rates[110] -= (fwd_rates[875] - rev_rates[847]);

  //rxn 876
  //sp 97
  sp_rates[96] -= (fwd_rates[876] - rev_rates[848]);
  //sp 98
  sp_rates[97] += (fwd_rates[876] - rev_rates[848]);
  //sp 110
  sp_rates[109] += (fwd_rates[876] - rev_rates[848]);
  //sp 111
  sp_rates[110] -= (fwd_rates[876] - rev_rates[848]);

  //rxn 877
  //sp 97
  sp_rates[96] -= (fwd_rates[877] - rev_rates[849]);
  //sp 98
  sp_rates[97] += (fwd_rates[877] - rev_rates[849]);
  //sp 110
  sp_rates[109] += (fwd_rates[877] - rev_rates[849]);
  //sp 111
  sp_rates[110] -= (fwd_rates[877] - rev_rates[849]);

  //rxn 878
  //sp 2
  (*dy_N) += (fwd_rates[878] - rev_rates[850]);
  //sp 12
  sp_rates[11] -= (fwd_rates[878] - rev_rates[850]);
  //sp 13
  sp_rates[12] += (fwd_rates[878] - rev_rates[850]);
  //sp 94
  sp_rates[93] -= (fwd_rates[878] - rev_rates[850]);

  //rxn 879
  //sp 12
  sp_rates[11] -= (fwd_rates[879] - rev_rates[851]);
  //sp 13
  sp_rates[12] += (fwd_rates[879] - rev_rates[851]);
  //sp 95
  sp_rates[94] -= (fwd_rates[879] - rev_rates[851]);
  //sp 93
  sp_rates[92] += (fwd_rates[879] - rev_rates[851]);

  //rxn 880
  //sp 99
  sp_rates[98] += (fwd_rates[880] - rev_rates[852]);
  //sp 12
  sp_rates[11] += (fwd_rates[880] - rev_rates[852]);
  //sp 93
  sp_rates[92] -= (fwd_rates[880] - rev_rates[852]);
  //sp 14
  sp_rates[13] -= (fwd_rates[880] - rev_rates[852]);

  //rxn 881
  //sp 5
  sp_rates[4] += (fwd_rates[881] - rev_rates[853]);
  //sp 12
  sp_rates[11] += (fwd_rates[881] - rev_rates[853]);
  //sp 14
  sp_rates[13] -= (fwd_rates[881] - rev_rates[853]);
  //sp 93
  sp_rates[92] += (fwd_rates[881] - rev_rates[853]);
  //sp 95
  sp_rates[94] -= (fwd_rates[881] - rev_rates[853]);

  //rxn 882
  //sp 3
  sp_rates[2] += (fwd_rates[882] - rev_rates[854]);
  //sp 13
  sp_rates[12] += (fwd_rates[882] - rev_rates[854]);
  //sp 14
  sp_rates[13] -= (fwd_rates[882] - rev_rates[854]);
  //sp 93
  sp_rates[92] += (fwd_rates[882] - rev_rates[854]);
  //sp 95
  sp_rates[94] -= (fwd_rates[882] - rev_rates[854]);

  //rxn 883
  //sp 100
  sp_rates[99] += (fwd_rates[883] - rev_rates[855]);
  //sp 12
  sp_rates[11] += (fwd_rates[883] - rev_rates[855]);
  //sp 14
  sp_rates[13] -= (fwd_rates[883] - rev_rates[855]);
  //sp 95
  sp_rates[94] -= (fwd_rates[883] - rev_rates[855]);

  //rxn 884
  //sp 99
  sp_rates[98] -= (fwd_rates[884] - rev_rates[856]);
  //sp 93
  sp_rates[92] += (fwd_rates[884] - rev_rates[856]);
  //sp 14
  sp_rates[13] -= (fwd_rates[884] - rev_rates[856]);
  //sp 15
  sp_rates[14] += (fwd_rates[884] - rev_rates[856]);

  //rxn 885
  //sp 12
  sp_rates[11] += (fwd_rates[885] - rev_rates[857]);
  //sp 92
  sp_rates[91] -= (fwd_rates[885] - rev_rates[857]);
  //sp 13
  sp_rates[12] -= (fwd_rates[885] - rev_rates[857]);
  //sp 93
  sp_rates[92] += (fwd_rates[885] - rev_rates[857]);

  //rxn 886
  //sp 95
  sp_rates[94] -= (fwd_rates[886] - rev_rates[858]);
  //sp 100
  sp_rates[99] += (fwd_rates[886] - rev_rates[858]);
  //sp 14
  sp_rates[13] += (fwd_rates[886] - rev_rates[858]);
  //sp 15
  sp_rates[14] -= (fwd_rates[886] - rev_rates[858]);

  //rxn 887
  //sp 95
  sp_rates[94] -= (fwd_rates[887] - rev_rates[859]);
  //sp 109
  sp_rates[108] += (fwd_rates[887] - rev_rates[859]);
  //sp 14
  sp_rates[13] += (fwd_rates[887] - rev_rates[859]);
  //sp 15
  sp_rates[14] -= (fwd_rates[887] - rev_rates[859]);

  //rxn 888
  //sp 14
  sp_rates[13] += (fwd_rates[888] - rev_rates[860]);
  //sp 15
  sp_rates[14] += (fwd_rates[888] - rev_rates[860]);
  //sp 48
  sp_rates[47] -= (fwd_rates[888] - rev_rates[860]);
  //sp 93
  sp_rates[92] += (fwd_rates[888] - rev_rates[860]);
  //sp 95
  sp_rates[94] -= (fwd_rates[888] - rev_rates[860]);

  //rxn 889
  //sp 3
  sp_rates[2] += (fwd_rates[889] - rev_rates[861]);
  //sp 2
  (*dy_N) -= (fwd_rates[889] - rev_rates[861]);
  //sp 115
  sp_rates[114] = (fwd_rates[889] - rev_rates[861]);
  //sp 16
  sp_rates[15] -= (fwd_rates[889] - rev_rates[861]);

  //rxn 890
  //sp 12
  sp_rates[11] += (fwd_rates[890] - rev_rates[862]);
  //sp 93
  sp_rates[92] -= (fwd_rates[890] - rev_rates[862]);
  //sp 96
  sp_rates[95] += (fwd_rates[890] - rev_rates[862]);
  //sp 16
  sp_rates[15] -= (fwd_rates[890] - rev_rates[862]);

  //rxn 891
  //sp 92
  sp_rates[91] += (fwd_rates[891] - rev_rates[863]);
  //sp 93
  sp_rates[92] -= (fwd_rates[891] - rev_rates[863]);
  //sp 14
  sp_rates[13] += (fwd_rates[891] - rev_rates[863]);
  //sp 16
  sp_rates[15] -= (fwd_rates[891] - rev_rates[863]);

  //rxn 892
  //sp 122
  sp_rates[121] += (fwd_rates[892] - rev_rates[864]);
  //sp 3
  sp_rates[2] += (fwd_rates[892] - rev_rates[864]);
  //sp 93
  sp_rates[92] -= (fwd_rates[892] - rev_rates[864]);
  //sp 16
  sp_rates[15] -= (fwd_rates[892] - rev_rates[864]);

  //rxn 893
  //sp 113
  sp_rates[112] = (fwd_rates[893] - rev_rates[865]);
  //sp 4
  sp_rates[3] += (fwd_rates[893] - rev_rates[865]);
  //sp 93
  sp_rates[92] -= (fwd_rates[893] - rev_rates[865]);
  //sp 16
  sp_rates[15] -= (fwd_rates[893] - rev_rates[865]);

  //rxn 894
  //sp 5
  sp_rates[4] += (fwd_rates[894] - rev_rates[866]);
  //sp 93
  sp_rates[92] -= (fwd_rates[894] - rev_rates[866]);
  //sp 112
  sp_rates[111] = (fwd_rates[894] - rev_rates[866]);
  //sp 16
  sp_rates[15] -= (fwd_rates[894] - rev_rates[866]);

  //rxn 895
  //sp 93
  sp_rates[92] += (fwd_rates[895] - rev_rates[867]);
  //sp 14
  sp_rates[13] += (fwd_rates[895] - rev_rates[867]);
  //sp 95
  sp_rates[94] -= (fwd_rates[895] - rev_rates[867]);
  //sp 16
  sp_rates[15] -= (fwd_rates[895] - rev_rates[867]);

  //rxn 896
  //sp 113
  sp_rates[112] += (fwd_rates[896] - rev_rates[868]);
  //sp 93
  sp_rates[92] += (fwd_rates[896] - rev_rates[868]);
  //sp 94
  sp_rates[93] -= (fwd_rates[896] - rev_rates[868]);
  //sp 16
  sp_rates[15] -= (fwd_rates[896] - rev_rates[868]);

  //rxn 897
  //sp 3
  sp_rates[2] += (fwd_rates[897] - rev_rates[869]);
  //sp 92
  sp_rates[91] -= (fwd_rates[897] - rev_rates[869]);
  //sp 112
  sp_rates[111] += (fwd_rates[897] - rev_rates[869]);
  //sp 16
  sp_rates[15] -= (fwd_rates[897] - rev_rates[869]);

  //rxn 898
  //sp 98
  sp_rates[97] -= (fwd_rates[898] - rev_rates[870]);
  //sp 3
  sp_rates[2] += 2.0 * (fwd_rates[898] - rev_rates[870]);
  //sp 118
  sp_rates[117] = (fwd_rates[898] - rev_rates[870]);
  //sp 16
  sp_rates[15] -= (fwd_rates[898] - rev_rates[870]);

  //rxn 899
  //sp 97
  sp_rates[96] -= (fwd_rates[899] - rev_rates[871]);
  //sp 3
  sp_rates[2] += (fwd_rates[899] - rev_rates[871]);
  //sp 118
  sp_rates[117] += (fwd_rates[899] - rev_rates[871]);
  //sp 16
  sp_rates[15] -= (fwd_rates[899] - rev_rates[871]);

  //rxn 900
  //sp 113
  sp_rates[112] += (fwd_rates[900] - rev_rates[872]);
  //sp 3
  sp_rates[2] += (fwd_rates[900] - rev_rates[872]);
  //sp 96
  sp_rates[95] -= (fwd_rates[900] - rev_rates[872]);
  //sp 16
  sp_rates[15] -= (fwd_rates[900] - rev_rates[872]);

  //rxn 901
  //sp 17
  sp_rates[16] -= (fwd_rates[901] - rev_rates[873]);
  //sp 4
  sp_rates[3] += (fwd_rates[901] - rev_rates[873]);
  //sp 93
  sp_rates[92] -= (fwd_rates[901] - rev_rates[873]);
  //sp 112
  sp_rates[111] += (fwd_rates[901] - rev_rates[873]);

  //rxn 902
  //sp 17
  sp_rates[16] -= (fwd_rates[902] - rev_rates[874]);
  //sp 92
  sp_rates[91] += (fwd_rates[902] - rev_rates[874]);
  //sp 12
  sp_rates[11] += (fwd_rates[902] - rev_rates[874]);
  //sp 93
  sp_rates[92] -= (fwd_rates[902] - rev_rates[874]);

  //rxn 903
  //sp 17
  sp_rates[16] -= (fwd_rates[903] - rev_rates[875]);
  //sp 93
  sp_rates[92] += (fwd_rates[903] - rev_rates[875]);
  //sp 94
  sp_rates[93] -= (fwd_rates[903] - rev_rates[875]);
  //sp 112
  sp_rates[111] += (fwd_rates[903] - rev_rates[875]);

  //rxn 904
  //sp 113
  sp_rates[112] += (fwd_rates[904] - rev_rates[876]);
  //sp 19
  sp_rates[18] -= (fwd_rates[904] - rev_rates[876]);
  //sp 93
  sp_rates[92] -= (fwd_rates[904] - rev_rates[876]);
  //sp 5
  sp_rates[4] += (fwd_rates[904] - rev_rates[876]);

  //rxn 905
  //sp 18
  sp_rates[17] += (fwd_rates[905] - rev_rates[877]);
  //sp 19
  sp_rates[18] -= (fwd_rates[905] - rev_rates[877]);

  //rxn 906
  //sp 19
  sp_rates[18] -= (fwd_rates[906] - rev_rates[878]);
  //sp 2
  (*dy_N) += (fwd_rates[906] - rev_rates[878]);
  //sp 94
  sp_rates[93] -= (fwd_rates[906] - rev_rates[878]);
  //sp 15
  sp_rates[14] += (fwd_rates[906] - rev_rates[878]);

  //rxn 907
  //sp 113
  sp_rates[112] -= (fwd_rates[907] - rev_rates[879]);
  //sp 19
  sp_rates[18] -= (fwd_rates[907] - rev_rates[879]);
  //sp 21
  sp_rates[20] += (fwd_rates[907] - rev_rates[879]);
  //sp 112
  sp_rates[111] += (fwd_rates[907] - rev_rates[879]);

  //rxn 908
  //sp 113
  sp_rates[112] += (fwd_rates[908] - rev_rates[880]);
  //sp 18
  sp_rates[17] -= (fwd_rates[908] - rev_rates[880]);
  //sp 3
  sp_rates[2] += (fwd_rates[908] - rev_rates[880]);
  //sp 92
  sp_rates[91] -= (fwd_rates[908] - rev_rates[880]);

  //rxn 909
  //sp 18
  sp_rates[17] -= (fwd_rates[909] - rev_rates[881]);
  //sp 3
  sp_rates[2] += (fwd_rates[909] - rev_rates[881]);
  //sp 93
  sp_rates[92] -= (fwd_rates[909] - rev_rates[881]);
  //sp 114
  sp_rates[113] = (fwd_rates[909] - rev_rates[881]);

  //rxn 910
  //sp 113
  sp_rates[112] += (fwd_rates[910] - rev_rates[882]);
  //sp 18
  sp_rates[17] -= (fwd_rates[910] - rev_rates[882]);
  //sp 93
  sp_rates[92] -= (fwd_rates[910] - rev_rates[882]);
  //sp 5
  sp_rates[4] += (fwd_rates[910] - rev_rates[882]);

  //rxn 911
  //sp 113
  sp_rates[112] += (fwd_rates[911] - rev_rates[883]);
  //sp 18
  sp_rates[17] -= (fwd_rates[911] - rev_rates[883]);
  //sp 2
  (*dy_N) -= (fwd_rates[911] - rev_rates[883]);
  //sp 96
  sp_rates[95] += (fwd_rates[911] - rev_rates[883]);

  //rxn 912
  //sp 18
  sp_rates[17] -= (fwd_rates[912] - rev_rates[884]);
  //sp 93
  sp_rates[92] += (fwd_rates[912] - rev_rates[884]);
  //sp 95
  sp_rates[94] -= (fwd_rates[912] - rev_rates[884]);
  //sp 15
  sp_rates[14] += (fwd_rates[912] - rev_rates[884]);

  //rxn 913
  //sp 113
  sp_rates[112] += (fwd_rates[913] - rev_rates[885]);
  //sp 21
  sp_rates[20] += (fwd_rates[913] - rev_rates[885]);
  //sp 22
  sp_rates[21] -= (fwd_rates[913] - rev_rates[885]);
  //sp 112
  sp_rates[111] -= (fwd_rates[913] - rev_rates[885]);

  //rxn 914
  //sp 113
  sp_rates[112] += (fwd_rates[914] - rev_rates[886]);
  //sp 21
  sp_rates[20] += (fwd_rates[914] - rev_rates[886]);
  //sp 22
  sp_rates[21] -= (fwd_rates[914] - rev_rates[886]);
  //sp 112
  sp_rates[111] -= (fwd_rates[914] - rev_rates[886]);

  //rxn 915
  //sp 97
  sp_rates[96] += (fwd_rates[915] - rev_rates[887]);
  //sp 21
  sp_rates[20] += (fwd_rates[915] - rev_rates[887]);
  //sp 22
  sp_rates[21] -= (fwd_rates[915] - rev_rates[887]);
  //sp 96
  sp_rates[95] -= (fwd_rates[915] - rev_rates[887]);

  //rxn 916
  //sp 97
  sp_rates[96] -= (fwd_rates[916] - rev_rates[888]);
  //sp 98
  sp_rates[97] += (fwd_rates[916] - rev_rates[888]);
  //sp 21
  sp_rates[20] += (fwd_rates[916] - rev_rates[888]);
  //sp 22
  sp_rates[21] -= (fwd_rates[916] - rev_rates[888]);

  //rxn 917
  //sp 100
  sp_rates[99] += (fwd_rates[917] - rev_rates[889]);
  //sp 21
  sp_rates[20] += (fwd_rates[917] - rev_rates[889]);
  //sp 22
  sp_rates[21] -= (fwd_rates[917] - rev_rates[889]);
  //sp 95
  sp_rates[94] -= (fwd_rates[917] - rev_rates[889]);

  //rxn 918
  //sp 100
  sp_rates[99] += (fwd_rates[918] - rev_rates[890]);
  //sp 21
  sp_rates[20] += (fwd_rates[918] - rev_rates[890]);
  //sp 22
  sp_rates[21] -= (fwd_rates[918] - rev_rates[890]);
  //sp 95
  sp_rates[94] -= (fwd_rates[918] - rev_rates[890]);

  //rxn 919
  //sp 109
  sp_rates[108] += (fwd_rates[919] - rev_rates[891]);
  //sp 21
  sp_rates[20] += (fwd_rates[919] - rev_rates[891]);
  //sp 22
  sp_rates[21] -= (fwd_rates[919] - rev_rates[891]);
  //sp 95
  sp_rates[94] -= (fwd_rates[919] - rev_rates[891]);

  //rxn 920
  //sp 3
  sp_rates[2] += (fwd_rates[920] - rev_rates[892]);
  //sp 92
  sp_rates[91] -= (fwd_rates[920] - rev_rates[892]);
  //sp 21
  sp_rates[20] -= (fwd_rates[920] - rev_rates[892]);
  //sp 118
  sp_rates[117] += (fwd_rates[920] - rev_rates[892]);

  //rxn 921
  //sp 113
  sp_rates[112] += (fwd_rates[921] - rev_rates[893]);
  //sp 10
  sp_rates[9] += (fwd_rates[921] - rev_rates[893]);
  //sp 21
  sp_rates[20] -= (fwd_rates[921] - rev_rates[893]);
  //sp 93
  sp_rates[92] -= (fwd_rates[921] - rev_rates[893]);

  //rxn 922
  //sp 5
  sp_rates[4] += (fwd_rates[922] - rev_rates[894]);
  //sp 21
  sp_rates[20] -= (fwd_rates[922] - rev_rates[894]);
  //sp 118
  sp_rates[117] += (fwd_rates[922] - rev_rates[894]);
  //sp 93
  sp_rates[92] -= (fwd_rates[922] - rev_rates[894]);

  //rxn 923
  //sp 26
  sp_rates[25] += (fwd_rates[923] - rev_rates[895]);
  //sp 21
  sp_rates[20] -= (fwd_rates[923] - rev_rates[895]);
  //sp 95
  sp_rates[94] -= (fwd_rates[923] - rev_rates[895]);
  //sp 93
  sp_rates[92] += (fwd_rates[923] - rev_rates[895]);

  //rxn 924
  //sp 99
  sp_rates[98] -= (fwd_rates[924] - rev_rates[896]);
  //sp 21
  sp_rates[20] -= (fwd_rates[924] - rev_rates[896]);
  //sp 22
  sp_rates[21] += (fwd_rates[924] - rev_rates[896]);
  //sp 93
  sp_rates[92] += (fwd_rates[924] - rev_rates[896]);

  //rxn 925
  //sp 100
  sp_rates[99] += (fwd_rates[925] - rev_rates[897]);
  //sp 30
  sp_rates[29] -= (fwd_rates[925] - rev_rates[897]);
  //sp 95
  sp_rates[94] -= (fwd_rates[925] - rev_rates[897]);
  //sp 24
  sp_rates[23] += (fwd_rates[925] - rev_rates[897]);

  //rxn 926
  //sp 109
  sp_rates[108] += (fwd_rates[926] - rev_rates[898]);
  //sp 30
  sp_rates[29] -= (fwd_rates[926] - rev_rates[898]);
  //sp 95
  sp_rates[94] -= (fwd_rates[926] - rev_rates[898]);
  //sp 24
  sp_rates[23] += (fwd_rates[926] - rev_rates[898]);

  //rxn 927
  //sp 26
  sp_rates[25] -= (fwd_rates[927] - rev_rates[899]);
  //sp 99
  sp_rates[98] += (fwd_rates[927] - rev_rates[899]);
  //sp 93
  sp_rates[92] -= (fwd_rates[927] - rev_rates[899]);
  //sp 15
  sp_rates[14] += (fwd_rates[927] - rev_rates[899]);

  //rxn 928
  //sp 26
  sp_rates[25] -= (fwd_rates[928] - rev_rates[900]);
  //sp 100
  sp_rates[99] += (fwd_rates[928] - rev_rates[900]);
  //sp 95
  sp_rates[94] -= (fwd_rates[928] - rev_rates[900]);
  //sp 15
  sp_rates[14] += (fwd_rates[928] - rev_rates[900]);

  //rxn 929
  //sp 26
  sp_rates[25] -= (fwd_rates[929] - rev_rates[901]);
  //sp 99
  sp_rates[98] -= (fwd_rates[929] - rev_rates[901]);
  //sp 93
  sp_rates[92] += (fwd_rates[929] - rev_rates[901]);
  //sp 30
  sp_rates[29] += (fwd_rates[929] - rev_rates[901]);

  //rxn 930
  //sp 26
  sp_rates[25] += (fwd_rates[930] - rev_rates[902]);
  //sp 93
  sp_rates[92] -= (fwd_rates[930] - rev_rates[902]);
  //sp 31
  sp_rates[30] -= (fwd_rates[930] - rev_rates[902]);
  //sp 95
  sp_rates[94] += (fwd_rates[930] - rev_rates[902]);

  //rxn 931
  //sp 12
  sp_rates[11] += (fwd_rates[931] - rev_rates[903]);
  //sp 43
  sp_rates[42] -= (fwd_rates[931] - rev_rates[903]);
  //sp 92
  sp_rates[91] -= (fwd_rates[931] - rev_rates[903]);
  //sp 113
  sp_rates[112] += (fwd_rates[931] - rev_rates[903]);

  //rxn 932
  //sp 114
  sp_rates[113] += (fwd_rates[932] - rev_rates[904]);
  //sp 43
  sp_rates[42] -= (fwd_rates[932] - rev_rates[904]);
  //sp 12
  sp_rates[11] += (fwd_rates[932] - rev_rates[904]);
  //sp 93
  sp_rates[92] -= (fwd_rates[932] - rev_rates[904]);

  //rxn 933
  //sp 113
  sp_rates[112] += (fwd_rates[933] - rev_rates[905]);
  //sp 43
  sp_rates[42] -= (fwd_rates[933] - rev_rates[905]);
  //sp 93
  sp_rates[92] -= (fwd_rates[933] - rev_rates[905]);
  //sp 13
  sp_rates[12] += (fwd_rates[933] - rev_rates[905]);

  //rxn 934
  //sp 114
  sp_rates[113] += (fwd_rates[934] - rev_rates[906]);
  //sp 43
  sp_rates[42] -= (fwd_rates[934] - rev_rates[906]);
  //sp 13
  sp_rates[12] += (fwd_rates[934] - rev_rates[906]);
  //sp 95
  sp_rates[94] -= (fwd_rates[934] - rev_rates[906]);

  //rxn 935
  //sp 5
  sp_rates[4] += fwd_rates[935];
  //sp 43
  sp_rates[42] -= fwd_rates[935];
  //sp 12
  sp_rates[11] += fwd_rates[935];
  //sp 122
  sp_rates[121] += fwd_rates[935];
  //sp 95
  sp_rates[94] -= fwd_rates[935];

  //rxn 936
  //sp 121
  sp_rates[120] = (fwd_rates[936] - rev_rates[907]);
  //sp 43
  sp_rates[42] -= (fwd_rates[936] - rev_rates[907]);
  //sp 13
  sp_rates[12] += (fwd_rates[936] - rev_rates[907]);
  //sp 95
  sp_rates[94] -= (fwd_rates[936] - rev_rates[907]);

  //rxn 937
  //sp 4
  sp_rates[3] += (fwd_rates[937] - rev_rates[908]);
  //sp 43
  sp_rates[42] -= (fwd_rates[937] - rev_rates[908]);
  //sp 13
  sp_rates[12] += (fwd_rates[937] - rev_rates[908]);
  //sp 113
  sp_rates[112] += (fwd_rates[937] - rev_rates[908]);
  //sp 95
  sp_rates[94] -= (fwd_rates[937] - rev_rates[908]);

  //rxn 938
  //sp 113
  sp_rates[112] += (fwd_rates[938] - rev_rates[909]);
  //sp 114
  sp_rates[113] -= (fwd_rates[938] - rev_rates[909]);
  //sp 4
  sp_rates[3] += (fwd_rates[938] - rev_rates[909]);

  //rxn 939
  //sp 113
  sp_rates[112] += (fwd_rates[939] - rev_rates[910]);
  //sp 114
  sp_rates[113] -= (fwd_rates[939] - rev_rates[910]);
  //sp 3
  sp_rates[2] -= (fwd_rates[939] - rev_rates[910]);
  //sp 5
  sp_rates[4] += (fwd_rates[939] - rev_rates[910]);

  //rxn 940
  //sp 114
  sp_rates[113] -= (fwd_rates[940] - rev_rates[911]);
  //sp 4
  sp_rates[3] -= (fwd_rates[940] - rev_rates[911]);
  //sp 93
  sp_rates[92] += (fwd_rates[940] - rev_rates[911]);
  //sp 14
  sp_rates[13] += (fwd_rates[940] - rev_rates[911]);

  //rxn 941
  //sp 114
  sp_rates[113] -= (fwd_rates[941] - rev_rates[912]);
  //sp 4
  sp_rates[3] -= (fwd_rates[941] - rev_rates[912]);
  //sp 5
  sp_rates[4] += (fwd_rates[941] - rev_rates[912]);
  //sp 122
  sp_rates[121] += (fwd_rates[941] - rev_rates[912]);

  //rxn 942
  //sp 101
  sp_rates[100] += (fwd_rates[942] - rev_rates[913]);
  //sp 114
  sp_rates[113] -= (fwd_rates[942] - rev_rates[913]);
  //sp 12
  sp_rates[11] += (fwd_rates[942] - rev_rates[913]);
  //sp 5
  sp_rates[4] -= (fwd_rates[942] - rev_rates[913]);

  //rxn 943
  //sp 114
  sp_rates[113] -= (fwd_rates[943] - rev_rates[914]);
  //sp 99
  sp_rates[98] += (fwd_rates[943] - rev_rates[914]);
  //sp 5
  sp_rates[4] -= (fwd_rates[943] - rev_rates[914]);
  //sp 14
  sp_rates[13] += (fwd_rates[943] - rev_rates[914]);

  //rxn 944
  //sp 113
  sp_rates[112] += (fwd_rates[944] - rev_rates[915]);
  //sp 114
  sp_rates[113] -= (fwd_rates[944] - rev_rates[915]);
  //sp 122
  sp_rates[121] += (fwd_rates[944] - rev_rates[915]);
  //sp 112
  sp_rates[111] -= (fwd_rates[944] - rev_rates[915]);

  //rxn 945
  //sp 12
  sp_rates[11] += fwd_rates[945];
  //sp 113
  sp_rates[112] += fwd_rates[945];
  //sp 114
  sp_rates[113] -= fwd_rates[945];
  //sp 122
  sp_rates[121] -= fwd_rates[945];
  //sp 93
  sp_rates[92] += fwd_rates[945];

  //rxn 946
  //sp 114
  sp_rates[113] -= (fwd_rates[946] - rev_rates[916]);
  //sp 124
  sp_rates[123] = (fwd_rates[946] - rev_rates[916]);
  //sp 93
  sp_rates[92] += (fwd_rates[946] - rev_rates[916]);
  //sp 96
  sp_rates[95] -= (fwd_rates[946] - rev_rates[916]);

  //rxn 947
  //sp 122
  sp_rates[121] += (fwd_rates[947] - rev_rates[917]);
  //sp 4
  sp_rates[3] -= (fwd_rates[947] - rev_rates[917]);
  //sp 5
  sp_rates[4] += (fwd_rates[947] - rev_rates[917]);
  //sp 120
  sp_rates[119] = -(fwd_rates[947] - rev_rates[917]);

  //rxn 948
  //sp 121
  sp_rates[120] += (fwd_rates[948] - rev_rates[918]);
  //sp 120
  sp_rates[119] -= (fwd_rates[948] - rev_rates[918]);

  //rxn 949
  //sp 97
  sp_rates[96] += (fwd_rates[949] - rev_rates[919]);
  //sp 3
  sp_rates[2] -= (fwd_rates[949] - rev_rates[919]);
  //sp 12
  sp_rates[11] += (fwd_rates[949] - rev_rates[919]);
  //sp 120
  sp_rates[119] -= (fwd_rates[949] - rev_rates[919]);

  //rxn 950
  //sp 122
  sp_rates[121] += (fwd_rates[950] - rev_rates[920]);
  //sp 3
  sp_rates[2] -= (fwd_rates[950] - rev_rates[920]);
  //sp 6
  sp_rates[5] += (fwd_rates[950] - rev_rates[920]);
  //sp 120
  sp_rates[119] -= (fwd_rates[950] - rev_rates[920]);

  //rxn 951
  //sp 10
  sp_rates[9] += (fwd_rates[951] - rev_rates[921]);
  //sp 5
  sp_rates[4] -= (fwd_rates[951] - rev_rates[921]);
  //sp 122
  sp_rates[121] += (fwd_rates[951] - rev_rates[921]);
  //sp 120
  sp_rates[119] -= (fwd_rates[951] - rev_rates[921]);

  //rxn 952
  //sp 97
  sp_rates[96] -= (fwd_rates[952] - rev_rates[922]);
  //sp 122
  sp_rates[121] += (fwd_rates[952] - rev_rates[922]);
  //sp 98
  sp_rates[97] += (fwd_rates[952] - rev_rates[922]);
  //sp 120
  sp_rates[119] -= (fwd_rates[952] - rev_rates[922]);

  //rxn 953
  //sp 121
  sp_rates[120] -= (fwd_rates[953] - rev_rates[923]) * pres_mod[98];
  //sp 12
  sp_rates[11] += (fwd_rates[953] - rev_rates[923]) * pres_mod[98];
  //sp 96
  sp_rates[95] += (fwd_rates[953] - rev_rates[923]) * pres_mod[98];

  //rxn 954
  //sp 121
  sp_rates[120] -= (fwd_rates[954] - rev_rates[924]);
  //sp 122
  sp_rates[121] += (fwd_rates[954] - rev_rates[924]);
  //sp 4
  sp_rates[3] -= (fwd_rates[954] - rev_rates[924]);
  //sp 5
  sp_rates[4] += (fwd_rates[954] - rev_rates[924]);

  //rxn 955
  //sp 121
  sp_rates[120] -= (fwd_rates[955] - rev_rates[925]);
  //sp 4
  sp_rates[3] -= (fwd_rates[955] - rev_rates[925]);
  //sp 13
  sp_rates[12] += (fwd_rates[955] - rev_rates[925]);
  //sp 96
  sp_rates[95] += (fwd_rates[955] - rev_rates[925]);

  //rxn 956
  //sp 121
  sp_rates[120] -= (fwd_rates[956] - rev_rates[926]);
  //sp 12
  sp_rates[11] += (fwd_rates[956] - rev_rates[926]);
  //sp 99
  sp_rates[98] += (fwd_rates[956] - rev_rates[926]);
  //sp 4
  sp_rates[3] -= (fwd_rates[956] - rev_rates[926]);

  //rxn 957
  //sp 121
  sp_rates[120] -= (fwd_rates[957] - rev_rates[927]);
  //sp 3
  sp_rates[2] -= (fwd_rates[957] - rev_rates[927]);
  //sp 12
  sp_rates[11] += (fwd_rates[957] - rev_rates[927]);
  //sp 97
  sp_rates[96] += (fwd_rates[957] - rev_rates[927]);

  //rxn 958
  //sp 121
  sp_rates[120] -= (fwd_rates[958] - rev_rates[928]);
  //sp 122
  sp_rates[121] += (fwd_rates[958] - rev_rates[928]);
  //sp 3
  sp_rates[2] -= (fwd_rates[958] - rev_rates[928]);
  //sp 6
  sp_rates[5] += (fwd_rates[958] - rev_rates[928]);

  //rxn 959
  //sp 121
  sp_rates[120] -= (fwd_rates[959] - rev_rates[929]);
  //sp 10
  sp_rates[9] += (fwd_rates[959] - rev_rates[929]);
  //sp 5
  sp_rates[4] -= (fwd_rates[959] - rev_rates[929]);
  //sp 122
  sp_rates[121] += (fwd_rates[959] - rev_rates[929]);

  //rxn 960
  //sp 121
  sp_rates[120] -= (fwd_rates[960] - rev_rates[930]);
  //sp 97
  sp_rates[96] += (fwd_rates[960] - rev_rates[930]);
  //sp 5
  sp_rates[4] -= (fwd_rates[960] - rev_rates[930]);
  //sp 13
  sp_rates[12] += (fwd_rates[960] - rev_rates[930]);

  //rxn 961
  //sp 121
  sp_rates[120] -= (fwd_rates[961] - rev_rates[931]);
  //sp 99
  sp_rates[98] += (fwd_rates[961] - rev_rates[931]);
  //sp 13
  sp_rates[12] += (fwd_rates[961] - rev_rates[931]);
  //sp 7
  sp_rates[6] -= (fwd_rates[961] - rev_rates[931]);

  //rxn 962
  //sp 121
  sp_rates[120] -= (fwd_rates[962] - rev_rates[932]);
  //sp 9
  sp_rates[8] += (fwd_rates[962] - rev_rates[932]);
  //sp 122
  sp_rates[121] += (fwd_rates[962] - rev_rates[932]);
  //sp 8
  sp_rates[7] -= (fwd_rates[962] - rev_rates[932]);

  //rxn 963
  //sp 121
  sp_rates[120] -= (fwd_rates[963] - rev_rates[933]);
  //sp 122
  sp_rates[121] += (fwd_rates[963] - rev_rates[933]);
  //sp 97
  sp_rates[96] += (fwd_rates[963] - rev_rates[933]);
  //sp 96
  sp_rates[95] -= (fwd_rates[963] - rev_rates[933]);

  //rxn 964
  //sp 121
  sp_rates[120] -= (fwd_rates[964] - rev_rates[934]);
  //sp 122
  sp_rates[121] += (fwd_rates[964] - rev_rates[934]);
  //sp 97
  sp_rates[96] -= (fwd_rates[964] - rev_rates[934]);
  //sp 98
  sp_rates[97] += (fwd_rates[964] - rev_rates[934]);

  //rxn 965
  //sp 121
  sp_rates[120] -= (fwd_rates[965] - rev_rates[935]);
  //sp 122
  sp_rates[121] += (fwd_rates[965] - rev_rates[935]);
  //sp 113
  sp_rates[112] += (fwd_rates[965] - rev_rates[935]);
  //sp 112
  sp_rates[111] -= (fwd_rates[965] - rev_rates[935]);

  //rxn 966
  //sp 113
  sp_rates[112] -= (fwd_rates[966] - rev_rates[936]) * pres_mod[99];
  //sp 116
  sp_rates[115] = (fwd_rates[966] - rev_rates[936]) * pres_mod[99];

  //rxn 967
  //sp 113
  sp_rates[112] += (fwd_rates[967] - rev_rates[937]);
  //sp 116
  sp_rates[115] -= (fwd_rates[967] - rev_rates[937]);

  //rxn 968
  //sp 4
  sp_rates[3] -= (fwd_rates[968] - rev_rates[938]);
  //sp 12
  sp_rates[11] += (fwd_rates[968] - rev_rates[938]);
  //sp 116
  sp_rates[115] -= (fwd_rates[968] - rev_rates[938]);
  //sp 96
  sp_rates[95] += (fwd_rates[968] - rev_rates[938]);

  //rxn 969
  //sp 121
  sp_rates[120] += (fwd_rates[969] - rev_rates[939]);
  //sp 3
  sp_rates[2] += (fwd_rates[969] - rev_rates[939]);
  //sp 116
  sp_rates[115] -= (fwd_rates[969] - rev_rates[939]);
  //sp 5
  sp_rates[4] -= (fwd_rates[969] - rev_rates[939]);

  //rxn 970
  //sp 3
  sp_rates[2] += (fwd_rates[970] - rev_rates[940]);
  //sp 116
  sp_rates[115] -= (fwd_rates[970] - rev_rates[940]);
  //sp 119
  sp_rates[118] = (fwd_rates[970] - rev_rates[940]);
  //sp 112
  sp_rates[111] -= (fwd_rates[970] - rev_rates[940]);

  //rxn 971
  //sp 113
  sp_rates[112] -= (fwd_rates[971] - rev_rates[941]);
  //sp 10
  sp_rates[9] += (fwd_rates[971] - rev_rates[941]);
  //sp 5
  sp_rates[4] -= (fwd_rates[971] - rev_rates[941]);
  //sp 112
  sp_rates[111] += (fwd_rates[971] - rev_rates[941]);

  //rxn 972
  //sp 113
  sp_rates[112] -= (fwd_rates[972] - rev_rates[942]);
  //sp 3
  sp_rates[2] += (fwd_rates[972] - rev_rates[942]);
  //sp 5
  sp_rates[4] -= (fwd_rates[972] - rev_rates[942]);
  //sp 120
  sp_rates[119] += (fwd_rates[972] - rev_rates[942]);

  //rxn 973
  //sp 113
  sp_rates[112] -= (fwd_rates[973] - rev_rates[943]);
  //sp 3
  sp_rates[2] += (fwd_rates[973] - rev_rates[943]);
  //sp 5
  sp_rates[4] -= (fwd_rates[973] - rev_rates[943]);
  //sp 121
  sp_rates[120] += (fwd_rates[973] - rev_rates[943]);

  //rxn 974
  //sp 113
  sp_rates[112] -= (fwd_rates[974] - rev_rates[944]);
  //sp 12
  sp_rates[11] += (fwd_rates[974] - rev_rates[944]);
  //sp 5
  sp_rates[4] -= (fwd_rates[974] - rev_rates[944]);
  //sp 97
  sp_rates[96] += (fwd_rates[974] - rev_rates[944]);

  //rxn 975
  //sp 113
  sp_rates[112] -= (fwd_rates[975] - rev_rates[945]);
  //sp 8
  sp_rates[7] += (fwd_rates[975] - rev_rates[945]);
  //sp 7
  sp_rates[6] -= (fwd_rates[975] - rev_rates[945]);
  //sp 112
  sp_rates[111] += (fwd_rates[975] - rev_rates[945]);

  //rxn 976
  //sp 113
  sp_rates[112] -= (fwd_rates[976] - rev_rates[946]);
  //sp 122
  sp_rates[121] += (fwd_rates[976] - rev_rates[946]);
  //sp 3
  sp_rates[2] += (fwd_rates[976] - rev_rates[946]);
  //sp 4
  sp_rates[3] -= (fwd_rates[976] - rev_rates[946]);

  //rxn 977
  //sp 113
  sp_rates[112] -= (fwd_rates[977] - rev_rates[947]);
  //sp 12
  sp_rates[11] += (fwd_rates[977] - rev_rates[947]);
  //sp 4
  sp_rates[3] -= (fwd_rates[977] - rev_rates[947]);
  //sp 96
  sp_rates[95] += (fwd_rates[977] - rev_rates[947]);

  //rxn 978
  //sp 113
  sp_rates[112] -= (fwd_rates[978] - rev_rates[948]) * pres_mod[100];
  //sp 3
  sp_rates[2] += (fwd_rates[978] - rev_rates[948]) * pres_mod[100];
  //sp 112
  sp_rates[111] += (fwd_rates[978] - rev_rates[948]) * pres_mod[100];

  //rxn 979
  //sp 113
  sp_rates[112] -= (fwd_rates[979] - rev_rates[949]) * pres_mod[101];
  //sp 3
  sp_rates[2] += (fwd_rates[979] - rev_rates[949]) * pres_mod[101];
  //sp 112
  sp_rates[111] += (fwd_rates[979] - rev_rates[949]) * pres_mod[101];

  //rxn 980
  //sp 17
  sp_rates[16] += (fwd_rates[980] - rev_rates[950]);
  //sp 92
  sp_rates[91] -= (fwd_rates[980] - rev_rates[950]);
  //sp 2
  (*dy_N) += (fwd_rates[980] - rev_rates[950]);
  //sp 112
  sp_rates[111] -= (fwd_rates[980] - rev_rates[950]);

  //rxn 981
  //sp 12
  sp_rates[11] += (fwd_rates[981] - rev_rates[951]);
  //sp 92
  sp_rates[91] += (fwd_rates[981] - rev_rates[951]);
  //sp 4
  sp_rates[3] -= (fwd_rates[981] - rev_rates[951]);
  //sp 112
  sp_rates[111] -= (fwd_rates[981] - rev_rates[951]);

  //rxn 982
  //sp 12
  sp_rates[11] += (fwd_rates[982] - rev_rates[952]);
  //sp 5
  sp_rates[4] -= (fwd_rates[982] - rev_rates[952]);
  //sp 96
  sp_rates[95] += (fwd_rates[982] - rev_rates[952]);
  //sp 112
  sp_rates[111] -= (fwd_rates[982] - rev_rates[952]);

  //rxn 983
  //sp 121
  sp_rates[120] += (fwd_rates[983] - rev_rates[953]);
  //sp 5
  sp_rates[4] -= (fwd_rates[983] - rev_rates[953]);
  //sp 112
  sp_rates[111] -= (fwd_rates[983] - rev_rates[953]);

  //rxn 984
  //sp 122
  sp_rates[121] += (fwd_rates[984] - rev_rates[954]);
  //sp 3
  sp_rates[2] += (fwd_rates[984] - rev_rates[954]);
  //sp 5
  sp_rates[4] -= (fwd_rates[984] - rev_rates[954]);
  //sp 112
  sp_rates[111] -= (fwd_rates[984] - rev_rates[954]);

  //rxn 985
  //sp 113
  sp_rates[112] += (fwd_rates[985] - rev_rates[955]);
  //sp 4
  sp_rates[3] += (fwd_rates[985] - rev_rates[955]);
  //sp 5
  sp_rates[4] -= (fwd_rates[985] - rev_rates[955]);
  //sp 112
  sp_rates[111] -= (fwd_rates[985] - rev_rates[955]);

  //rxn 986
  //sp 122
  sp_rates[121] += (fwd_rates[986] - rev_rates[956]);
  //sp 92
  sp_rates[91] += (fwd_rates[986] - rev_rates[956]);
  //sp 93
  sp_rates[92] -= (fwd_rates[986] - rev_rates[956]);
  //sp 112
  sp_rates[111] -= (fwd_rates[986] - rev_rates[956]);

  //rxn 987
  //sp 2
  (*dy_N) += (fwd_rates[987] - rev_rates[957]);
  //sp 12
  sp_rates[11] += (fwd_rates[987] - rev_rates[957]);
  //sp 93
  sp_rates[92] -= (fwd_rates[987] - rev_rates[957]);
  //sp 112
  sp_rates[111] -= (fwd_rates[987] - rev_rates[957]);

  //rxn 988
  //sp 113
  sp_rates[112] += (fwd_rates[988] - rev_rates[958]);
  //sp 3
  sp_rates[2] += (fwd_rates[988] - rev_rates[958]);
  //sp 6
  sp_rates[5] -= (fwd_rates[988] - rev_rates[958]);
  //sp 112
  sp_rates[111] -= (fwd_rates[988] - rev_rates[958]);

  //rxn 989
  //sp 12
  sp_rates[11] += (fwd_rates[989] - rev_rates[959]);
  //sp 93
  sp_rates[92] += (fwd_rates[989] - rev_rates[959]);
  //sp 7
  sp_rates[6] -= (fwd_rates[989] - rev_rates[959]);
  //sp 112
  sp_rates[111] -= (fwd_rates[989] - rev_rates[959]);

  //rxn 990
  //sp 113
  sp_rates[112] += (fwd_rates[990] - rev_rates[960]);
  //sp 99
  sp_rates[98] -= (fwd_rates[990] - rev_rates[960]);
  //sp 93
  sp_rates[92] += (fwd_rates[990] - rev_rates[960]);
  //sp 112
  sp_rates[111] -= (fwd_rates[990] - rev_rates[960]);

  //rxn 991
  //sp 113
  sp_rates[112] -= (fwd_rates[991] - rev_rates[961]);
  //sp 3
  sp_rates[2] += (fwd_rates[991] - rev_rates[961]);
  //sp 119
  sp_rates[118] += (fwd_rates[991] - rev_rates[961]);
  //sp 112
  sp_rates[111] -= (fwd_rates[991] - rev_rates[961]);

  //rxn 992
  //sp 122
  sp_rates[121] += (fwd_rates[992] - rev_rates[962]);
  //sp 2
  (*dy_N) += (fwd_rates[992] - rev_rates[962]);
  //sp 94
  sp_rates[93] -= (fwd_rates[992] - rev_rates[962]);
  //sp 112
  sp_rates[111] -= (fwd_rates[992] - rev_rates[962]);

  //rxn 993
  //sp 115
  sp_rates[114] += (fwd_rates[993] - rev_rates[963]);
  //sp 93
  sp_rates[92] += (fwd_rates[993] - rev_rates[963]);
  //sp 94
  sp_rates[93] -= (fwd_rates[993] - rev_rates[963]);
  //sp 112
  sp_rates[111] -= (fwd_rates[993] - rev_rates[963]);

  //rxn 994
  //sp 122
  sp_rates[121] += (fwd_rates[994] - rev_rates[964]);
  //sp 12
  sp_rates[11] += (fwd_rates[994] - rev_rates[964]);
  //sp 13
  sp_rates[12] -= (fwd_rates[994] - rev_rates[964]);
  //sp 112
  sp_rates[111] -= (fwd_rates[994] - rev_rates[964]);

  //rxn 995
  //sp 122
  sp_rates[121] += (fwd_rates[995] - rev_rates[965]);
  //sp 93
  sp_rates[92] += (fwd_rates[995] - rev_rates[965]);
  //sp 95
  sp_rates[94] -= (fwd_rates[995] - rev_rates[965]);
  //sp 112
  sp_rates[111] -= (fwd_rates[995] - rev_rates[965]);

  //rxn 996
  //sp 12
  sp_rates[11] += (fwd_rates[996] - rev_rates[966]);
  //sp 94
  sp_rates[93] += (fwd_rates[996] - rev_rates[966]);
  //sp 95
  sp_rates[94] -= (fwd_rates[996] - rev_rates[966]);
  //sp 112
  sp_rates[111] -= (fwd_rates[996] - rev_rates[966]);

  //rxn 997
  //sp 2
  (*dy_N) += (fwd_rates[997] - rev_rates[967]);
  //sp 13
  sp_rates[12] += (fwd_rates[997] - rev_rates[967]);
  //sp 95
  sp_rates[94] -= (fwd_rates[997] - rev_rates[967]);
  //sp 112
  sp_rates[111] -= (fwd_rates[997] - rev_rates[967]);

  //rxn 998
  //sp 113
  sp_rates[112] += (fwd_rates[998] - rev_rates[968]);
  //sp 14
  sp_rates[13] += (fwd_rates[998] - rev_rates[968]);
  //sp 15
  sp_rates[14] -= (fwd_rates[998] - rev_rates[968]);
  //sp 112
  sp_rates[111] -= (fwd_rates[998] - rev_rates[968]);

  //rxn 999
  //sp 113
  sp_rates[112] += (fwd_rates[999] - rev_rates[969]);
  //sp 100
  sp_rates[99] -= (fwd_rates[999] - rev_rates[969]);
  //sp 95
  sp_rates[94] += (fwd_rates[999] - rev_rates[969]);
  //sp 112
  sp_rates[111] -= (fwd_rates[999] - rev_rates[969]);

  //rxn 1000
  //sp 119
  sp_rates[118] -= (fwd_rates[1000] - rev_rates[970]) * pres_mod[102];
  //sp 112
  sp_rates[111] += 2.0 * (fwd_rates[1000] - rev_rates[970]) * pres_mod[102];

  //rxn 1001
  //sp 122
  sp_rates[121] += (fwd_rates[1001] - rev_rates[971]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1001] - rev_rates[971]);
  //sp 119
  sp_rates[118] -= (fwd_rates[1001] - rev_rates[971]);
  //sp 112
  sp_rates[111] += (fwd_rates[1001] - rev_rates[971]);

  //rxn 1002
  //sp 5
  sp_rates[4] -= (fwd_rates[1002] - rev_rates[972]);
  //sp 120
  sp_rates[119] += (fwd_rates[1002] - rev_rates[972]);
  //sp 119
  sp_rates[118] -= (fwd_rates[1002] - rev_rates[972]);
  //sp 112
  sp_rates[111] += (fwd_rates[1002] - rev_rates[972]);

  //rxn 1003
  //sp 12
  sp_rates[11] += (fwd_rates[1003] - rev_rates[973]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1003] - rev_rates[973]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1003] - rev_rates[973]);
  //sp 93
  sp_rates[92] += (fwd_rates[1003] - rev_rates[973]);

  //rxn 1004
  //sp 122
  sp_rates[121] -= (fwd_rates[1004] - rev_rates[974]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1004] - rev_rates[974]);
  //sp 7
  sp_rates[6] += (fwd_rates[1004] - rev_rates[974]);
  //sp 112
  sp_rates[111] += (fwd_rates[1004] - rev_rates[974]);

  //rxn 1005
  //sp 12
  sp_rates[11] += (fwd_rates[1005] - rev_rates[975]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1005] - rev_rates[975]);
  //sp 92
  sp_rates[91] -= (fwd_rates[1005] - rev_rates[975]);
  //sp 2
  (*dy_N) += (fwd_rates[1005] - rev_rates[975]);

  //rxn 1006
  //sp 122
  sp_rates[121] -= (fwd_rates[1006] - rev_rates[976]);
  //sp 3
  sp_rates[2] -= (fwd_rates[1006] - rev_rates[976]);
  //sp 12
  sp_rates[11] += (fwd_rates[1006] - rev_rates[976]);
  //sp 96
  sp_rates[95] += (fwd_rates[1006] - rev_rates[976]);

  //rxn 1007
  //sp 122
  sp_rates[121] -= (fwd_rates[1007] - rev_rates[977]);
  //sp 3
  sp_rates[2] -= (fwd_rates[1007] - rev_rates[977]);
  //sp 12
  sp_rates[11] += (fwd_rates[1007] - rev_rates[977]);
  //sp 96
  sp_rates[95] += (fwd_rates[1007] - rev_rates[977]);

  //rxn 1008
  //sp 3
  sp_rates[2] += (fwd_rates[1008] - rev_rates[978]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1008] - rev_rates[978]);
  //sp 12
  sp_rates[11] += (fwd_rates[1008] - rev_rates[978]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1008] - rev_rates[978]);
  //sp 93
  sp_rates[92] += (fwd_rates[1008] - rev_rates[978]);

  //rxn 1009
  //sp 122
  sp_rates[121] -= (fwd_rates[1009] - rev_rates[979]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1009] - rev_rates[979]);
  //sp 14
  sp_rates[13] += (fwd_rates[1009] - rev_rates[979]);
  //sp 93
  sp_rates[92] += (fwd_rates[1009] - rev_rates[979]);

  //rxn 1010
  //sp 122
  sp_rates[121] -= (fwd_rates[1010] - rev_rates[980]);
  //sp 13
  sp_rates[12] += (fwd_rates[1010] - rev_rates[980]);
  //sp 7
  sp_rates[6] -= (fwd_rates[1010] - rev_rates[980]);
  //sp 93
  sp_rates[92] += (fwd_rates[1010] - rev_rates[980]);

  //rxn 1011
  //sp 122
  sp_rates[121] -= (fwd_rates[1011] - rev_rates[981]);
  //sp 115
  sp_rates[114] += (fwd_rates[1011] - rev_rates[981]);
  //sp 12
  sp_rates[11] += (fwd_rates[1011] - rev_rates[981]);
  //sp 112
  sp_rates[111] -= (fwd_rates[1011] - rev_rates[981]);

  //rxn 1012
  //sp 122
  sp_rates[121] -= (fwd_rates[1012] - rev_rates[982]);
  //sp 12
  sp_rates[11] += (fwd_rates[1012] - rev_rates[982]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1012] - rev_rates[982]);
  //sp 94
  sp_rates[93] += (fwd_rates[1012] - rev_rates[982]);

  //rxn 1013
  //sp 2
  (*dy_N) += (fwd_rates[1013] - rev_rates[983]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1013] - rev_rates[983]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1013] - rev_rates[983]);
  //sp 13
  sp_rates[12] += (fwd_rates[1013] - rev_rates[983]);

  //rxn 1014
  //sp 121
  sp_rates[120] += (fwd_rates[1014] - rev_rates[984]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1014] - rev_rates[984]);
  //sp 7
  sp_rates[6] += (fwd_rates[1014] - rev_rates[984]);
  //sp 8
  sp_rates[7] -= (fwd_rates[1014] - rev_rates[984]);

  //rxn 1015
  //sp 122
  sp_rates[121] -= (fwd_rates[1015] - rev_rates[985]);
  //sp 12
  sp_rates[11] += (fwd_rates[1015] - rev_rates[985]);
  //sp 93
  sp_rates[92] += 2.0 * (fwd_rates[1015] - rev_rates[985]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1015] - rev_rates[985]);

  //rxn 1016
  //sp 122
  sp_rates[121] -= (fwd_rates[1016] - rev_rates[986]);
  //sp 13
  sp_rates[12] += (fwd_rates[1016] - rev_rates[986]);
  //sp 94
  sp_rates[93] += (fwd_rates[1016] - rev_rates[986]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1016] - rev_rates[986]);

  //rxn 1017
  //sp 12
  sp_rates[11] += (fwd_rates[1017] - rev_rates[987]);
  //sp 93
  sp_rates[92] += (fwd_rates[1017] - rev_rates[987]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1017] - rev_rates[987]);
  //sp 2
  (*dy_N) += (fwd_rates[1017] - rev_rates[987]);
  //sp 94
  sp_rates[93] -= (fwd_rates[1017] - rev_rates[987]);

  //rxn 1018
  //sp 121
  sp_rates[120] += (fwd_rates[1018] - rev_rates[988]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1018] - rev_rates[988]);
  //sp 12
  sp_rates[11] += (fwd_rates[1018] - rev_rates[988]);
  //sp 14
  sp_rates[13] -= (fwd_rates[1018] - rev_rates[988]);

  //rxn 1019
  //sp 122
  sp_rates[121] -= 2.0 * (fwd_rates[1019] - rev_rates[989]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[1019] - rev_rates[989]);
  //sp 2
  (*dy_N) += (fwd_rates[1019] - rev_rates[989]);

  //rxn 1020
  //sp 121
  sp_rates[120] += (fwd_rates[1020] - rev_rates[990]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1020] - rev_rates[990]);
  //sp 99
  sp_rates[98] -= (fwd_rates[1020] - rev_rates[990]);
  //sp 93
  sp_rates[92] += (fwd_rates[1020] - rev_rates[990]);

  //rxn 1021
  //sp 121
  sp_rates[120] += (fwd_rates[1021] - rev_rates[991]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1021] - rev_rates[991]);
  //sp 14
  sp_rates[13] += (fwd_rates[1021] - rev_rates[991]);
  //sp 15
  sp_rates[14] -= (fwd_rates[1021] - rev_rates[991]);

  //rxn 1022
  //sp 121
  sp_rates[120] += (fwd_rates[1022] - rev_rates[992]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1022] - rev_rates[992]);
  //sp 100
  sp_rates[99] -= (fwd_rates[1022] - rev_rates[992]);
  //sp 95
  sp_rates[94] += (fwd_rates[1022] - rev_rates[992]);

  //rxn 1023
  //sp 122
  sp_rates[121] -= (fwd_rates[1023] - rev_rates[993]) * pres_mod[103];
  //sp 92
  sp_rates[91] += (fwd_rates[1023] - rev_rates[993]) * pres_mod[103];
  //sp 12
  sp_rates[11] += (fwd_rates[1023] - rev_rates[993]) * pres_mod[103];

  //rxn 1024
  //sp 121
  sp_rates[120] += (fwd_rates[1024] - rev_rates[994]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1024] - rev_rates[994]);
  //sp 21
  sp_rates[20] += (fwd_rates[1024] - rev_rates[994]);
  //sp 22
  sp_rates[21] -= (fwd_rates[1024] - rev_rates[994]);

  //rxn 1025
  //sp 121
  sp_rates[120] += (fwd_rates[1025] - rev_rates[995]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1025] - rev_rates[995]);
  //sp 29
  sp_rates[28] -= (fwd_rates[1025] - rev_rates[995]);
  //sp 34
  sp_rates[33] += (fwd_rates[1025] - rev_rates[995]);

  //rxn 1026
  //sp 17
  sp_rates[16] += (fwd_rates[1026] - rev_rates[996]) * pres_mod[104];
  //sp 115
  sp_rates[114] -= (fwd_rates[1026] - rev_rates[996]) * pres_mod[104];
  //sp 2
  (*dy_N) += (fwd_rates[1026] - rev_rates[996]) * pres_mod[104];

  //rxn 1027
  //sp 113
  sp_rates[112] += (fwd_rates[1027] - rev_rates[997]);
  //sp 3
  sp_rates[2] -= (fwd_rates[1027] - rev_rates[997]);
  //sp 92
  sp_rates[91] += (fwd_rates[1027] - rev_rates[997]);
  //sp 115
  sp_rates[114] -= (fwd_rates[1027] - rev_rates[997]);

  //rxn 1028
  //sp 115
  sp_rates[114] -= (fwd_rates[1028] - rev_rates[998]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1028] - rev_rates[998]);
  //sp 93
  sp_rates[92] += (fwd_rates[1028] - rev_rates[998]);
  //sp 112
  sp_rates[111] += (fwd_rates[1028] - rev_rates[998]);

  //rxn 1029
  //sp 17
  sp_rates[16] -= (fwd_rates[1029] - rev_rates[999]);
  //sp 115
  sp_rates[114] -= (fwd_rates[1029] - rev_rates[999]);
  //sp 112
  sp_rates[111] += 2.0 * (fwd_rates[1029] - rev_rates[999]);

  //rxn 1030
  //sp 115
  sp_rates[114] -= (fwd_rates[1030] - rev_rates[1000]);
  //sp 92
  sp_rates[91] -= (fwd_rates[1030] - rev_rates[1000]);
  //sp 2
  (*dy_N) += (fwd_rates[1030] - rev_rates[1000]);
  //sp 112
  sp_rates[111] += (fwd_rates[1030] - rev_rates[1000]);

  //rxn 1031
  //sp 113
  sp_rates[112] += (fwd_rates[1031] - rev_rates[1001]);
  //sp 115
  sp_rates[114] -= (fwd_rates[1031] - rev_rates[1001]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1031] - rev_rates[1001]);
  //sp 93
  sp_rates[92] += (fwd_rates[1031] - rev_rates[1001]);

  //rxn 1032
  //sp 122
  sp_rates[121] += (fwd_rates[1032] - rev_rates[1002]);
  //sp 115
  sp_rates[114] -= (fwd_rates[1032] - rev_rates[1002]);
  //sp 93
  sp_rates[92] += (fwd_rates[1032] - rev_rates[1002]);
  //sp 7
  sp_rates[6] -= (fwd_rates[1032] - rev_rates[1002]);

  //rxn 1033
  //sp 115
  sp_rates[114] -= (fwd_rates[1033] - rev_rates[1003]);
  //sp 3
  sp_rates[2] += (fwd_rates[1033] - rev_rates[1003]);
  //sp 117
  sp_rates[116] = (fwd_rates[1033] - rev_rates[1003]);
  //sp 6
  sp_rates[5] -= (fwd_rates[1033] - rev_rates[1003]);

  //rxn 1034
  //sp 115
  sp_rates[114] -= (fwd_rates[1034] - rev_rates[1004]);
  //sp 92
  sp_rates[91] += (fwd_rates[1034] - rev_rates[1004]);
  //sp 119
  sp_rates[118] += (fwd_rates[1034] - rev_rates[1004]);
  //sp 112
  sp_rates[111] -= (fwd_rates[1034] - rev_rates[1004]);

  //rxn 1035
  //sp 115
  sp_rates[114] -= 2.0 * (fwd_rates[1035] - rev_rates[1005]);
  //sp 2
  (*dy_N) += (fwd_rates[1035] - rev_rates[1005]);
  //sp 112
  sp_rates[111] += 2.0 * (fwd_rates[1035] - rev_rates[1005]);

  //rxn 1036
  //sp 113
  sp_rates[112] += (fwd_rates[1036] - rev_rates[1006]);
  //sp 115
  sp_rates[114] -= (fwd_rates[1036] - rev_rates[1006]);
  //sp 112
  sp_rates[111] += (fwd_rates[1036] - rev_rates[1006]);
  //sp 16
  sp_rates[15] -= (fwd_rates[1036] - rev_rates[1006]);

  //rxn 1037
  //sp 18
  sp_rates[17] -= (fwd_rates[1037] - rev_rates[1007]);
  //sp 115
  sp_rates[114] -= (fwd_rates[1037] - rev_rates[1007]);
  //sp 118
  sp_rates[117] += (fwd_rates[1037] - rev_rates[1007]);
  //sp 112
  sp_rates[111] += (fwd_rates[1037] - rev_rates[1007]);

  //rxn 1038
  //sp 3
  sp_rates[2] -= (fwd_rates[1038] - rev_rates[1008]);
  //sp 117
  sp_rates[116] += (fwd_rates[1038] - rev_rates[1008]);
  //sp 115
  sp_rates[114] -= (fwd_rates[1038] - rev_rates[1008]);

  //rxn 1039
  //sp 116
  sp_rates[115] += (fwd_rates[1039] - rev_rates[1009]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1039] - rev_rates[1009]);
  //sp 117
  sp_rates[116] -= (fwd_rates[1039] - rev_rates[1009]);
  //sp 93
  sp_rates[92] += (fwd_rates[1039] - rev_rates[1009]);

  //rxn 1040
  //sp 122
  sp_rates[121] += (fwd_rates[1040] - rev_rates[1010]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1040] - rev_rates[1010]);
  //sp 117
  sp_rates[116] -= (fwd_rates[1040] - rev_rates[1010]);
  //sp 96
  sp_rates[95] += (fwd_rates[1040] - rev_rates[1010]);

  //rxn 1041
  //sp 99
  sp_rates[98] += (fwd_rates[1041] - rev_rates[1011]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1041] - rev_rates[1011]);
  //sp 117
  sp_rates[116] -= (fwd_rates[1041] - rev_rates[1011]);
  //sp 112
  sp_rates[111] += (fwd_rates[1041] - rev_rates[1011]);

  //rxn 1042
  //sp 115
  sp_rates[114] += (fwd_rates[1042] - rev_rates[1012]);
  //sp 117
  sp_rates[116] -= (fwd_rates[1042] - rev_rates[1012]);
  //sp 7
  sp_rates[6] -= (fwd_rates[1042] - rev_rates[1012]);
  //sp 8
  sp_rates[7] += (fwd_rates[1042] - rev_rates[1012]);

  //rxn 1043
  //sp 10
  sp_rates[9] += (fwd_rates[1043] - rev_rates[1013]);
  //sp 115
  sp_rates[114] += (fwd_rates[1043] - rev_rates[1013]);
  //sp 117
  sp_rates[116] -= (fwd_rates[1043] - rev_rates[1013]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1043] - rev_rates[1013]);

  //rxn 1044
  //sp 18
  sp_rates[17] += (fwd_rates[1044] - rev_rates[1014]);
  //sp 92
  sp_rates[91] -= (fwd_rates[1044] - rev_rates[1014]);
  //sp 2
  (*dy_N) += (fwd_rates[1044] - rev_rates[1014]);
  //sp 118
  sp_rates[117] -= (fwd_rates[1044] - rev_rates[1014]);

  //rxn 1045
  //sp 113
  sp_rates[112] += (fwd_rates[1045] - rev_rates[1015]) * pres_mod[105];
  //sp 3
  sp_rates[2] += (fwd_rates[1045] - rev_rates[1015]) * pres_mod[105];
  //sp 118
  sp_rates[117] -= (fwd_rates[1045] - rev_rates[1015]) * pres_mod[105];

  //rxn 1046
  //sp 113
  sp_rates[112] += (fwd_rates[1046] - rev_rates[1016]);
  //sp 3
  sp_rates[2] -= (fwd_rates[1046] - rev_rates[1016]);
  //sp 6
  sp_rates[5] += (fwd_rates[1046] - rev_rates[1016]);
  //sp 118
  sp_rates[117] -= (fwd_rates[1046] - rev_rates[1016]);

  //rxn 1047
  //sp 113
  sp_rates[112] += (fwd_rates[1047] - rev_rates[1017]);
  //sp 4
  sp_rates[3] -= (fwd_rates[1047] - rev_rates[1017]);
  //sp 5
  sp_rates[4] += (fwd_rates[1047] - rev_rates[1017]);
  //sp 118
  sp_rates[117] -= (fwd_rates[1047] - rev_rates[1017]);

  //rxn 1048
  //sp 113
  sp_rates[112] += (fwd_rates[1048] - rev_rates[1018]);
  //sp 10
  sp_rates[9] += (fwd_rates[1048] - rev_rates[1018]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1048] - rev_rates[1018]);
  //sp 118
  sp_rates[117] -= (fwd_rates[1048] - rev_rates[1018]);

  //rxn 1049
  //sp 113
  sp_rates[112] += (fwd_rates[1049] - rev_rates[1019]);
  //sp 10
  sp_rates[9] += (fwd_rates[1049] - rev_rates[1019]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1049] - rev_rates[1019]);
  //sp 118
  sp_rates[117] -= (fwd_rates[1049] - rev_rates[1019]);

  //rxn 1050
  //sp 93
  sp_rates[92] += (fwd_rates[1050] - rev_rates[1020]);
  //sp 118
  sp_rates[117] -= (fwd_rates[1050] - rev_rates[1020]);
  //sp 7
  sp_rates[6] -= (fwd_rates[1050] - rev_rates[1020]);
  //sp 15
  sp_rates[14] += (fwd_rates[1050] - rev_rates[1020]);

  //rxn 1051
  //sp 113
  sp_rates[112] += (fwd_rates[1051] - rev_rates[1021]);
  //sp 97
  sp_rates[96] += (fwd_rates[1051] - rev_rates[1021]);
  //sp 118
  sp_rates[117] -= (fwd_rates[1051] - rev_rates[1021]);
  //sp 96
  sp_rates[95] -= (fwd_rates[1051] - rev_rates[1021]);

  //rxn 1052
  //sp 97
  sp_rates[96] -= (fwd_rates[1052] - rev_rates[1022]);
  //sp 98
  sp_rates[97] += (fwd_rates[1052] - rev_rates[1022]);
  //sp 113
  sp_rates[112] += (fwd_rates[1052] - rev_rates[1022]);
  //sp 118
  sp_rates[117] -= (fwd_rates[1052] - rev_rates[1022]);

  //rxn 1053
  //sp 105
  sp_rates[104] -= (fwd_rates[1053] - rev_rates[1023]);
  //sp 99
  sp_rates[98] += (fwd_rates[1053] - rev_rates[1023]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1053] - rev_rates[1023]);
  //sp 104
  sp_rates[103] += (fwd_rates[1053] - rev_rates[1023]);

  //rxn 1054
  //sp 105
  sp_rates[104] -= (fwd_rates[1054] - rev_rates[1024]);
  //sp 100
  sp_rates[99] += (fwd_rates[1054] - rev_rates[1024]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1054] - rev_rates[1024]);
  //sp 104
  sp_rates[103] += (fwd_rates[1054] - rev_rates[1024]);

  //rxn 1055
  //sp 105
  sp_rates[104] -= (fwd_rates[1055] - rev_rates[1025]);
  //sp 109
  sp_rates[108] += (fwd_rates[1055] - rev_rates[1025]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1055] - rev_rates[1025]);
  //sp 104
  sp_rates[103] += (fwd_rates[1055] - rev_rates[1025]);

  //rxn 1056
  //sp 105
  sp_rates[104] -= (fwd_rates[1056] - rev_rates[1026]);
  //sp 106
  sp_rates[105] += (fwd_rates[1056] - rev_rates[1026]);
  //sp 6
  sp_rates[5] += (fwd_rates[1056] - rev_rates[1026]);

  //rxn 1057
  //sp 103
  sp_rates[102] += (fwd_rates[1057] - rev_rates[1027]);
  //sp 100
  sp_rates[99] += (fwd_rates[1057] - rev_rates[1027]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1057] - rev_rates[1027]);
  //sp 104
  sp_rates[103] -= (fwd_rates[1057] - rev_rates[1027]);

  //rxn 1058
  //sp 95
  sp_rates[94] -= (fwd_rates[1058] - rev_rates[1028]);
  //sp 100
  sp_rates[99] += (fwd_rates[1058] - rev_rates[1028]);
  //sp 102
  sp_rates[101] += (fwd_rates[1058] - rev_rates[1028]);
  //sp 103
  sp_rates[102] -= (fwd_rates[1058] - rev_rates[1028]);

  //rxn 1059
  //sp 100
  sp_rates[99] += (fwd_rates[1059] - rev_rates[1029]);
  //sp 2
  (*dy_N) += (fwd_rates[1059] - rev_rates[1029]);
  //sp 102
  sp_rates[101] -= (fwd_rates[1059] - rev_rates[1029]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1059] - rev_rates[1029]);

  //rxn 1060
  //sp 113
  sp_rates[112] += (fwd_rates[1060] - rev_rates[1030]);
  //sp 12
  sp_rates[11] += (fwd_rates[1060] - rev_rates[1030]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1060] - rev_rates[1030]);
  //sp 40
  sp_rates[39] -= (fwd_rates[1060] - rev_rates[1030]);

  //rxn 1061
  //sp 14
  sp_rates[13] += (fwd_rates[1061] - rev_rates[1031]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1061] - rev_rates[1031]);
  //sp 112
  sp_rates[111] += (fwd_rates[1061] - rev_rates[1031]);
  //sp 40
  sp_rates[39] -= (fwd_rates[1061] - rev_rates[1031]);

  //rxn 1062
  //sp 43
  sp_rates[42] += (fwd_rates[1062] - rev_rates[1032]);
  //sp 93
  sp_rates[92] += (fwd_rates[1062] - rev_rates[1032]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1062] - rev_rates[1032]);
  //sp 40
  sp_rates[39] -= (fwd_rates[1062] - rev_rates[1032]);

  //rxn 1063
  //sp 97
  sp_rates[96] += (fwd_rates[1063] - rev_rates[1033]);
  //sp 98
  sp_rates[97] -= (fwd_rates[1063] - rev_rates[1033]);
  //sp 20
  sp_rates[19] += (fwd_rates[1063] - rev_rates[1033]);
  //sp 40
  sp_rates[39] -= (fwd_rates[1063] - rev_rates[1033]);

  //rxn 1064
  //sp 113
  sp_rates[112] -= (fwd_rates[1064] - rev_rates[1034]);
  //sp 20
  sp_rates[19] += (fwd_rates[1064] - rev_rates[1034]);
  //sp 112
  sp_rates[111] += (fwd_rates[1064] - rev_rates[1034]);
  //sp 40
  sp_rates[39] -= (fwd_rates[1064] - rev_rates[1034]);

  //rxn 1065
  //sp 113
  sp_rates[112] += (fwd_rates[1065] - rev_rates[1035]);
  //sp 92
  sp_rates[91] -= (fwd_rates[1065] - rev_rates[1035]);
  //sp 20
  sp_rates[19] -= (fwd_rates[1065] - rev_rates[1035]);
  //sp 16
  sp_rates[15] += (fwd_rates[1065] - rev_rates[1035]);

  //rxn 1066
  //sp 113
  sp_rates[112] += (fwd_rates[1066] - rev_rates[1036]);
  //sp 122
  sp_rates[121] -= (fwd_rates[1066] - rev_rates[1036]);
  //sp 43
  sp_rates[42] += (fwd_rates[1066] - rev_rates[1036]);
  //sp 20
  sp_rates[19] -= (fwd_rates[1066] - rev_rates[1036]);

  //rxn 1067
  //sp 33
  sp_rates[32] -= (fwd_rates[1067] - rev_rates[1037]);
  //sp 18
  sp_rates[17] += (fwd_rates[1067] - rev_rates[1037]);
  //sp 92
  sp_rates[91] -= (fwd_rates[1067] - rev_rates[1037]);
  //sp 113
  sp_rates[112] += (fwd_rates[1067] - rev_rates[1037]);

  //rxn 1068
  //sp 43
  sp_rates[42] += (fwd_rates[1068] - rev_rates[1038]);
  //sp 2
  (*dy_N) += (fwd_rates[1068] - rev_rates[1038]);
  //sp 94
  sp_rates[93] -= (fwd_rates[1068] - rev_rates[1038]);
  //sp 40
  sp_rates[39] -= (fwd_rates[1068] - rev_rates[1038]);

  //rxn 1069
  //sp 33
  sp_rates[32] -= (fwd_rates[1069] - rev_rates[1039]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1069] - rev_rates[1039]);
  //sp 113
  sp_rates[112] += (fwd_rates[1069] - rev_rates[1039]);
  //sp 15
  sp_rates[14] += (fwd_rates[1069] - rev_rates[1039]);

  //rxn 1070
  //sp 33
  sp_rates[32] -= (fwd_rates[1070] - rev_rates[1040]);
  //sp 93
  sp_rates[92] += (fwd_rates[1070] - rev_rates[1040]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1070] - rev_rates[1040]);
  //sp 48
  sp_rates[47] += (fwd_rates[1070] - rev_rates[1040]);

  //rxn 1071
  //sp 20
  sp_rates[19] += (fwd_rates[1071] - rev_rates[1041]);
  //sp 92
  sp_rates[91] -= (fwd_rates[1071] - rev_rates[1041]);
  //sp 45
  sp_rates[44] -= (fwd_rates[1071] - rev_rates[1041]);
  //sp 113
  sp_rates[112] += (fwd_rates[1071] - rev_rates[1041]);

  //rxn 1072
  //sp 121
  sp_rates[120] += (fwd_rates[1072] - rev_rates[1042]);
  //sp 10
  sp_rates[9] += (fwd_rates[1072] - rev_rates[1042]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1072] - rev_rates[1042]);
  //sp 24
  sp_rates[23] -= (fwd_rates[1072] - rev_rates[1042]);

  //rxn 1073
  //sp 100
  sp_rates[99] += (fwd_rates[1073] - rev_rates[1043]);
  //sp 15
  sp_rates[14] += (fwd_rates[1073] - rev_rates[1043]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1073] - rev_rates[1043]);
  //sp 24
  sp_rates[23] -= (fwd_rates[1073] - rev_rates[1043]);

  //rxn 1074
  //sp 99
  sp_rates[98] -= (fwd_rates[1074] - rev_rates[1044]);
  //sp 93
  sp_rates[92] += (fwd_rates[1074] - rev_rates[1044]);
  //sp 30
  sp_rates[29] += (fwd_rates[1074] - rev_rates[1044]);
  //sp 24
  sp_rates[23] -= (fwd_rates[1074] - rev_rates[1044]);

  //rxn 1075
  //sp 99
  sp_rates[98] += (fwd_rates[1075] - rev_rates[1045]);
  //sp 77
  sp_rates[76] -= (fwd_rates[1075] - rev_rates[1045]);
  //sp 78
  sp_rates[77] += (fwd_rates[1075] - rev_rates[1045]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1075] - rev_rates[1045]);

  //rxn 1076
  //sp 100
  sp_rates[99] += (fwd_rates[1076] - rev_rates[1046]);
  //sp 77
  sp_rates[76] -= (fwd_rates[1076] - rev_rates[1046]);
  //sp 78
  sp_rates[77] += (fwd_rates[1076] - rev_rates[1046]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1076] - rev_rates[1046]);

  //rxn 1077
  //sp 81
  sp_rates[80] += (fwd_rates[1077] - rev_rates[1047]);
  //sp 93
  sp_rates[92] += (fwd_rates[1077] - rev_rates[1047]);
  //sp 78
  sp_rates[77] -= (fwd_rates[1077] - rev_rates[1047]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1077] - rev_rates[1047]);

  //rxn 1078
  //sp 81
  sp_rates[80] -= (fwd_rates[1078] - rev_rates[1048]);
  //sp 82
  sp_rates[81] += (fwd_rates[1078] - rev_rates[1048]);
  //sp 100
  sp_rates[99] += (fwd_rates[1078] - rev_rates[1048]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1078] - rev_rates[1048]);

  //rxn 1079
  //sp 81
  sp_rates[80] -= (fwd_rates[1079] - rev_rates[1049]);
  //sp 82
  sp_rates[81] += (fwd_rates[1079] - rev_rates[1049]);
  //sp 99
  sp_rates[98] += (fwd_rates[1079] - rev_rates[1049]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1079] - rev_rates[1049]);

  //rxn 1080
  //sp 81
  sp_rates[80] += (fwd_rates[1080] - rev_rates[1050]);
  //sp 95
  sp_rates[94] += (fwd_rates[1080] - rev_rates[1050]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1080] - rev_rates[1050]);
  //sp 79
  sp_rates[78] -= (fwd_rates[1080] - rev_rates[1050]);

  //rxn 1081
  //sp 95
  sp_rates[94] += (fwd_rates[1081] - rev_rates[1051]);
  //sp 100
  sp_rates[99] -= (fwd_rates[1081] - rev_rates[1051]);
  //sp 79
  sp_rates[78] -= (fwd_rates[1081] - rev_rates[1051]);
  //sp 80
  sp_rates[79] += (fwd_rates[1081] - rev_rates[1051]);

  //rxn 1082
  //sp 58
  sp_rates[57] += (fwd_rates[1082] - rev_rates[1052]);
  //sp 59
  sp_rates[58] -= (fwd_rates[1082] - rev_rates[1052]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1082] - rev_rates[1052]);
  //sp 95
  sp_rates[94] += (fwd_rates[1082] - rev_rates[1052]);

  //rxn 1083
  //sp 65
  sp_rates[64] += (fwd_rates[1083] - rev_rates[1053]);
  //sp 99
  sp_rates[98] += (fwd_rates[1083] - rev_rates[1053]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1083] - rev_rates[1053]);
  //sp 47
  sp_rates[46] -= (fwd_rates[1083] - rev_rates[1053]);

  //rxn 1084
  //sp 47
  sp_rates[46] += (fwd_rates[1084] - rev_rates[1054]);
  //sp 100
  sp_rates[99] += (fwd_rates[1084] - rev_rates[1054]);
  //sp 61
  sp_rates[60] -= (fwd_rates[1084] - rev_rates[1054]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1084] - rev_rates[1054]);

  //rxn 1085
  //sp 47
  sp_rates[46] += (fwd_rates[1085] - rev_rates[1055]);
  //sp 109
  sp_rates[108] += (fwd_rates[1085] - rev_rates[1055]);
  //sp 61
  sp_rates[60] -= (fwd_rates[1085] - rev_rates[1055]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1085] - rev_rates[1055]);

  //rxn 1086
  //sp 109
  sp_rates[108] += (fwd_rates[1086] - rev_rates[1056]);
  //sp 62
  sp_rates[61] = -(fwd_rates[1086] - rev_rates[1056]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1086] - rev_rates[1056]);
  //sp 64
  sp_rates[63] = (fwd_rates[1086] - rev_rates[1056]);

  //rxn 1087
  //sp 63
  sp_rates[62] = (fwd_rates[1087] - rev_rates[1057]);
  //sp 109
  sp_rates[108] += (fwd_rates[1087] - rev_rates[1057]);
  //sp 62
  sp_rates[61] -= (fwd_rates[1087] - rev_rates[1057]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1087] - rev_rates[1057]);

  //rxn 1088
  //sp 41
  sp_rates[40] -= (fwd_rates[1088] - rev_rates[1058]);
  //sp 42
  sp_rates[41] += (fwd_rates[1088] - rev_rates[1058]);
  //sp 92
  sp_rates[91] += (fwd_rates[1088] - rev_rates[1058]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1088] - rev_rates[1058]);

  //rxn 1089
  //sp 41
  sp_rates[40] -= (fwd_rates[1089] - rev_rates[1059]);
  //sp 2
  (*dy_N) -= (fwd_rates[1089] - rev_rates[1059]);
  //sp 112
  sp_rates[111] += 2.0 * (fwd_rates[1089] - rev_rates[1059]);

  //rxn 1090
  //sp 42
  sp_rates[41] -= (fwd_rates[1090] - rev_rates[1060]);
  //sp 12
  sp_rates[11] += (fwd_rates[1090] - rev_rates[1060]);
  //sp 93
  sp_rates[92] -= (fwd_rates[1090] - rev_rates[1060]);
  //sp 122
  sp_rates[121] += (fwd_rates[1090] - rev_rates[1060]);

  //rxn 1091
  //sp 42
  sp_rates[41] -= (fwd_rates[1091] - rev_rates[1061]);
  //sp 13
  sp_rates[12] += (fwd_rates[1091] - rev_rates[1061]);
  //sp 95
  sp_rates[94] -= (fwd_rates[1091] - rev_rates[1061]);
  //sp 122
  sp_rates[121] += (fwd_rates[1091] - rev_rates[1061]);

  //rxn 1092
  //sp 113
  sp_rates[112] += (fwd_rates[1092] - rev_rates[1062]);
  //sp 3
  sp_rates[2] += (fwd_rates[1092] - rev_rates[1062]);
  //sp 124
  sp_rates[123] -= (fwd_rates[1092] - rev_rates[1062]);

  //rxn 1093
  //sp 124
  sp_rates[123] -= (fwd_rates[1093] - rev_rates[1063]);
  //sp 118
  sp_rates[117] += (fwd_rates[1093] - rev_rates[1063]);

  //rxn 1094
  //sp 113
  sp_rates[112] += (fwd_rates[1094] - rev_rates[1064]);
  //sp 3
  sp_rates[2] -= (fwd_rates[1094] - rev_rates[1064]);
  //sp 124
  sp_rates[123] -= (fwd_rates[1094] - rev_rates[1064]);
  //sp 6
  sp_rates[5] += (fwd_rates[1094] - rev_rates[1064]);

  //rxn 1095
  //sp 4
  sp_rates[3] -= (fwd_rates[1095] - rev_rates[1065]);
  //sp 3
  sp_rates[2] += (fwd_rates[1095] - rev_rates[1065]);
  //sp 124
  sp_rates[123] -= (fwd_rates[1095] - rev_rates[1065]);
  //sp 121
  sp_rates[120] += (fwd_rates[1095] - rev_rates[1065]);

  //rxn 1096
  //sp 4
  sp_rates[3] -= (fwd_rates[1096] - rev_rates[1066]);
  //sp 124
  sp_rates[123] -= (fwd_rates[1096] - rev_rates[1066]);
  //sp 113
  sp_rates[112] += (fwd_rates[1096] - rev_rates[1066]);
  //sp 5
  sp_rates[4] += (fwd_rates[1096] - rev_rates[1066]);

  //rxn 1097
  //sp 113
  sp_rates[112] += (fwd_rates[1097] - rev_rates[1067]);
  //sp 10
  sp_rates[9] += (fwd_rates[1097] - rev_rates[1067]);
  //sp 124
  sp_rates[123] -= (fwd_rates[1097] - rev_rates[1067]);
  //sp 5
  sp_rates[4] -= (fwd_rates[1097] - rev_rates[1067]);

  //rxn 1098
  //sp 113
  sp_rates[112] += (fwd_rates[1098] - rev_rates[1068]);
  //sp 124
  sp_rates[123] -= (fwd_rates[1098] - rev_rates[1068]);
  //sp 21
  sp_rates[20] -= (fwd_rates[1098] - rev_rates[1068]);
  //sp 22
  sp_rates[21] += (fwd_rates[1098] - rev_rates[1068]);

  //sp 0
  sp_rates[0] = 0.0;
  //sp 1
  sp_rates[1] = 0.0;
} // end eval_spec_rates

