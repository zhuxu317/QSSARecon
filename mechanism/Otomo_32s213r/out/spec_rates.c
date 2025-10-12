#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 3
  sp_rates[3] = -(fwd_rates[0] - rev_rates[0]);
  //sp 4
  sp_rates[4] = -(fwd_rates[0] - rev_rates[0]);
  //sp 5
  sp_rates[5] = (fwd_rates[0] - rev_rates[0]);
  //sp 6
  sp_rates[6] = (fwd_rates[0] - rev_rates[0]);

  //rxn 1
  //sp 2
  sp_rates[2] = -(fwd_rates[1] - rev_rates[1]);
  //sp 4
  sp_rates[4] += (fwd_rates[1] - rev_rates[1]);
  //sp 5
  sp_rates[5] -= (fwd_rates[1] - rev_rates[1]);
  //sp 6
  sp_rates[6] += (fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 2
  sp_rates[2] -= (fwd_rates[2] - rev_rates[2]);
  //sp 4
  sp_rates[4] += (fwd_rates[2] - rev_rates[2]);
  //sp 5
  sp_rates[5] -= (fwd_rates[2] - rev_rates[2]);
  //sp 6
  sp_rates[6] += (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 9
  sp_rates[9] = (fwd_rates[3] - rev_rates[3]);
  //sp 2
  sp_rates[2] -= (fwd_rates[3] - rev_rates[3]);
  //sp 4
  sp_rates[4] += (fwd_rates[3] - rev_rates[3]);
  //sp 6
  sp_rates[6] -= (fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 9
  sp_rates[9] += (fwd_rates[4] - rev_rates[4]);
  //sp 5
  sp_rates[5] += (fwd_rates[4] - rev_rates[4]);
  //sp 6
  sp_rates[6] -= 2.0 * (fwd_rates[4] - rev_rates[4]);

  //rxn 5
  //sp 9
  sp_rates[9] += (fwd_rates[5] - rev_rates[5]);
  //sp 5
  sp_rates[5] += (fwd_rates[5] - rev_rates[5]);
  //sp 6
  sp_rates[6] -= 2.0 * (fwd_rates[5] - rev_rates[5]);

  //rxn 6
  //sp 2
  sp_rates[2] -= (fwd_rates[6] - rev_rates[6]) * pres_mod[0];
  //sp 4
  sp_rates[4] += 2.0 * (fwd_rates[6] - rev_rates[6]) * pres_mod[0];

  //rxn 7
  //sp 2
  sp_rates[2] -= (fwd_rates[7] - rev_rates[7]) * pres_mod[1];
  //sp 4
  sp_rates[4] += 2.0 * (fwd_rates[7] - rev_rates[7]) * pres_mod[1];

  //rxn 8
  //sp 4
  sp_rates[4] -= (fwd_rates[8] - rev_rates[8]) * pres_mod[2];
  //sp 5
  sp_rates[5] -= (fwd_rates[8] - rev_rates[8]) * pres_mod[2];
  //sp 6
  sp_rates[6] += (fwd_rates[8] - rev_rates[8]) * pres_mod[2];

  //rxn 9
  //sp 3
  sp_rates[3] += (fwd_rates[9] - rev_rates[9]) * pres_mod[3];
  //sp 5
  sp_rates[5] -= 2.0 * (fwd_rates[9] - rev_rates[9]) * pres_mod[3];

  //rxn 10
  //sp 9
  sp_rates[9] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[4];
  //sp 4
  sp_rates[4] += (fwd_rates[10] - rev_rates[10]) * pres_mod[4];
  //sp 6
  sp_rates[6] += (fwd_rates[10] - rev_rates[10]) * pres_mod[4];

  //rxn 11
  //sp 9
  sp_rates[9] -= (fwd_rates[11] - rev_rates[11]) * pres_mod[5];
  //sp 4
  sp_rates[4] += (fwd_rates[11] - rev_rates[11]) * pres_mod[5];
  //sp 6
  sp_rates[6] += (fwd_rates[11] - rev_rates[11]) * pres_mod[5];

  //rxn 12
  //sp 3
  sp_rates[3] -= (fwd_rates[12] - rev_rates[12]) * pres_mod[6];
  //sp 4
  sp_rates[4] -= (fwd_rates[12] - rev_rates[12]) * pres_mod[6];
  //sp 7
  sp_rates[7] = (fwd_rates[12] - rev_rates[12]) * pres_mod[6];

  //rxn 13
  //sp 2
  sp_rates[2] += (fwd_rates[13] - rev_rates[13]);
  //sp 3
  sp_rates[3] += (fwd_rates[13] - rev_rates[13]);
  //sp 4
  sp_rates[4] -= (fwd_rates[13] - rev_rates[13]);
  //sp 7
  sp_rates[7] -= (fwd_rates[13] - rev_rates[13]);

  //rxn 14
  //sp 4
  sp_rates[4] -= (fwd_rates[14] - rev_rates[14]);
  //sp 6
  sp_rates[6] += 2.0 * (fwd_rates[14] - rev_rates[14]);
  //sp 7
  sp_rates[7] -= (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 9
  sp_rates[9] += (fwd_rates[15] - rev_rates[15]);
  //sp 4
  sp_rates[4] -= (fwd_rates[15] - rev_rates[15]);
  //sp 5
  sp_rates[5] += (fwd_rates[15] - rev_rates[15]);
  //sp 7
  sp_rates[7] -= (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 3
  sp_rates[3] += (fwd_rates[16] - rev_rates[16]);
  //sp 5
  sp_rates[5] -= (fwd_rates[16] - rev_rates[16]);
  //sp 6
  sp_rates[6] += (fwd_rates[16] - rev_rates[16]);
  //sp 7
  sp_rates[7] -= (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 9
  sp_rates[9] += (fwd_rates[17] - rev_rates[17]);
  //sp 3
  sp_rates[3] += (fwd_rates[17] - rev_rates[17]);
  //sp 6
  sp_rates[6] -= (fwd_rates[17] - rev_rates[17]);
  //sp 7
  sp_rates[7] -= (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 9
  sp_rates[9] += (fwd_rates[18] - rev_rates[18]);
  //sp 3
  sp_rates[3] += (fwd_rates[18] - rev_rates[18]);
  //sp 6
  sp_rates[6] -= (fwd_rates[18] - rev_rates[18]);
  //sp 7
  sp_rates[7] -= (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 10
  sp_rates[10] = (fwd_rates[19] - rev_rates[19]);
  //sp 3
  sp_rates[3] += (fwd_rates[19] - rev_rates[19]);
  //sp 7
  sp_rates[7] -= 2.0 * (fwd_rates[19] - rev_rates[19]);

  //rxn 20
  //sp 10
  sp_rates[10] += (fwd_rates[20] - rev_rates[20]);
  //sp 3
  sp_rates[3] += (fwd_rates[20] - rev_rates[20]);
  //sp 7
  sp_rates[7] -= 2.0 * (fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 10
  sp_rates[10] -= (fwd_rates[21] - rev_rates[21]) * pres_mod[7];
  //sp 6
  sp_rates[6] += 2.0 * (fwd_rates[21] - rev_rates[21]) * pres_mod[7];

  //rxn 22
  //sp 9
  sp_rates[9] += (fwd_rates[22] - rev_rates[22]);
  //sp 10
  sp_rates[10] -= (fwd_rates[22] - rev_rates[22]);
  //sp 4
  sp_rates[4] -= (fwd_rates[22] - rev_rates[22]);
  //sp 6
  sp_rates[6] += (fwd_rates[22] - rev_rates[22]);

  //rxn 23
  //sp 10
  sp_rates[10] -= (fwd_rates[23] - rev_rates[23]);
  //sp 2
  sp_rates[2] += (fwd_rates[23] - rev_rates[23]);
  //sp 4
  sp_rates[4] -= (fwd_rates[23] - rev_rates[23]);
  //sp 7
  sp_rates[7] += (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 10
  sp_rates[10] -= (fwd_rates[24] - rev_rates[24]);
  //sp 5
  sp_rates[5] -= (fwd_rates[24] - rev_rates[24]);
  //sp 6
  sp_rates[6] += (fwd_rates[24] - rev_rates[24]);
  //sp 7
  sp_rates[7] += (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 9
  sp_rates[9] += (fwd_rates[25] - rev_rates[25]);
  //sp 10
  sp_rates[10] -= (fwd_rates[25] - rev_rates[25]);
  //sp 6
  sp_rates[6] -= (fwd_rates[25] - rev_rates[25]);
  //sp 7
  sp_rates[7] += (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 9
  sp_rates[9] += (fwd_rates[26] - rev_rates[26]);
  //sp 10
  sp_rates[10] -= (fwd_rates[26] - rev_rates[26]);
  //sp 6
  sp_rates[6] -= (fwd_rates[26] - rev_rates[26]);
  //sp 7
  sp_rates[7] += (fwd_rates[26] - rev_rates[26]);

  //rxn 27
  //sp 1
  sp_rates[1] = -(fwd_rates[27] - rev_rates[27]) * pres_mod[8];
  //sp 11
  sp_rates[11] = (fwd_rates[27] - rev_rates[27]) * pres_mod[8];
  //sp 4
  sp_rates[4] += (fwd_rates[27] - rev_rates[27]) * pres_mod[8];

  //rxn 28
  //sp 1
  sp_rates[1] -= (fwd_rates[28] - rev_rates[28]);
  //sp 2
  sp_rates[2] += (fwd_rates[28] - rev_rates[28]);
  //sp 11
  sp_rates[11] += (fwd_rates[28] - rev_rates[28]);
  //sp 4
  sp_rates[4] -= (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 1
  sp_rates[1] -= (fwd_rates[29] - rev_rates[29]);
  //sp 11
  sp_rates[11] += (fwd_rates[29] - rev_rates[29]);
  //sp 5
  sp_rates[5] -= (fwd_rates[29] - rev_rates[29]);
  //sp 6
  sp_rates[6] += (fwd_rates[29] - rev_rates[29]);

  //rxn 30
  //sp 1
  sp_rates[1] -= (fwd_rates[30] - rev_rates[30]);
  //sp 11
  sp_rates[11] += (fwd_rates[30] - rev_rates[30]);
  //sp 6
  sp_rates[6] -= (fwd_rates[30] - rev_rates[30]);
  //sp 9
  sp_rates[9] += (fwd_rates[30] - rev_rates[30]);

  //rxn 31
  //sp 1
  sp_rates[1] -= (fwd_rates[31] - rev_rates[31]);
  //sp 10
  sp_rates[10] += (fwd_rates[31] - rev_rates[31]);
  //sp 11
  sp_rates[11] += (fwd_rates[31] - rev_rates[31]);
  //sp 7
  sp_rates[7] -= (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 2
  sp_rates[2] += (fwd_rates[32] - rev_rates[32]);
  //sp 11
  sp_rates[11] -= (fwd_rates[32] - rev_rates[32]);
  //sp 4
  sp_rates[4] -= (fwd_rates[32] - rev_rates[32]);
  //sp 12
  sp_rates[12] = (fwd_rates[32] - rev_rates[32]);

  //rxn 33
  //sp 18
  sp_rates[18] = (fwd_rates[33] - rev_rates[33]);
  //sp 11
  sp_rates[11] -= (fwd_rates[33] - rev_rates[33]);
  //sp 4
  sp_rates[4] += (fwd_rates[33] - rev_rates[33]);
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
  //sp 11
  sp_rates[11] -= (fwd_rates[35] - rev_rates[35]);
  //sp 12
  sp_rates[12] += (fwd_rates[35] - rev_rates[35]);
  //sp 5
  sp_rates[5] -= (fwd_rates[35] - rev_rates[35]);
  //sp 6
  sp_rates[6] += (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 9
  sp_rates[9] += (fwd_rates[36] - rev_rates[36]);
  //sp 11
  sp_rates[11] -= (fwd_rates[36] - rev_rates[36]);
  //sp 12
  sp_rates[12] += (fwd_rates[36] - rev_rates[36]);
  //sp 6
  sp_rates[6] -= (fwd_rates[36] - rev_rates[36]);

  //rxn 37
  //sp 3
  sp_rates[3] += (fwd_rates[37] - rev_rates[37]);
  //sp 1
  sp_rates[1] += (fwd_rates[37] - rev_rates[37]);
  //sp 11
  sp_rates[11] -= (fwd_rates[37] - rev_rates[37]);
  //sp 7
  sp_rates[7] -= (fwd_rates[37] - rev_rates[37]);

  //rxn 38
  //sp 16
  sp_rates[16] = (fwd_rates[38] - rev_rates[38]);
  //sp 11
  sp_rates[11] -= (fwd_rates[38] - rev_rates[38]);
  //sp 6
  sp_rates[6] += (fwd_rates[38] - rev_rates[38]);
  //sp 7
  sp_rates[7] -= (fwd_rates[38] - rev_rates[38]);

  //rxn 39
  //sp 9
  sp_rates[9] += (fwd_rates[39] - rev_rates[39]);
  //sp 18
  sp_rates[18] += (fwd_rates[39] - rev_rates[39]);
  //sp 11
  sp_rates[11] -= (fwd_rates[39] - rev_rates[39]);
  //sp 7
  sp_rates[7] -= (fwd_rates[39] - rev_rates[39]);

  //rxn 40
  //sp 9
  sp_rates[9] += (fwd_rates[40] - rev_rates[40]);
  //sp 18
  sp_rates[18] += (fwd_rates[40] - rev_rates[40]);
  //sp 11
  sp_rates[11] -= (fwd_rates[40] - rev_rates[40]);
  //sp 7
  sp_rates[7] -= (fwd_rates[40] - rev_rates[40]);

  //rxn 41
  //sp 19
  sp_rates[19] = (fwd_rates[41] - rev_rates[41]);
  //sp 9
  sp_rates[9] += (fwd_rates[41] - rev_rates[41]);
  //sp 11
  sp_rates[11] -= (fwd_rates[41] - rev_rates[41]);
  //sp 7
  sp_rates[7] -= (fwd_rates[41] - rev_rates[41]);

  //rxn 42
  //sp 3
  sp_rates[3] -= (fwd_rates[42] - rev_rates[42]);
  //sp 16
  sp_rates[16] += (fwd_rates[42] - rev_rates[42]);
  //sp 11
  sp_rates[11] -= (fwd_rates[42] - rev_rates[42]);
  //sp 5
  sp_rates[5] += (fwd_rates[42] - rev_rates[42]);

  //rxn 43
  //sp 3
  sp_rates[3] -= (fwd_rates[43] - rev_rates[43]);
  //sp 18
  sp_rates[18] += (fwd_rates[43] - rev_rates[43]);
  //sp 11
  sp_rates[11] -= (fwd_rates[43] - rev_rates[43]);
  //sp 6
  sp_rates[6] += (fwd_rates[43] - rev_rates[43]);

  //rxn 44
  //sp 1
  sp_rates[1] += (fwd_rates[44] - rev_rates[44]);
  //sp 11
  sp_rates[11] -= 2.0 * (fwd_rates[44] - rev_rates[44]);
  //sp 12
  sp_rates[12] += (fwd_rates[44] - rev_rates[44]);

  //rxn 45
  //sp 1
  sp_rates[1] += (fwd_rates[45] - rev_rates[45]);
  //sp 11
  sp_rates[11] -= (fwd_rates[45] - rev_rates[45]);
  //sp 12
  sp_rates[12] -= (fwd_rates[45] - rev_rates[45]);
  //sp 13
  sp_rates[13] = (fwd_rates[45] - rev_rates[45]);

  //rxn 46
  //sp 8
  sp_rates[8] = (fwd_rates[46] - rev_rates[46]);
  //sp 11
  sp_rates[11] -= (fwd_rates[46] - rev_rates[46]);
  //sp 4
  sp_rates[4] += 2.0 * (fwd_rates[46] - rev_rates[46]);
  //sp 13
  sp_rates[13] -= (fwd_rates[46] - rev_rates[46]);

  //rxn 47
  //sp 0
  sp_rates[0] = (fwd_rates[47] - rev_rates[47]);
  //sp 1
  sp_rates[1] += (fwd_rates[47] - rev_rates[47]);
  //sp 18
  sp_rates[18] -= (fwd_rates[47] - rev_rates[47]);
  //sp 11
  sp_rates[11] -= (fwd_rates[47] - rev_rates[47]);

  //rxn 48
  //sp 0
  sp_rates[0] -= (fwd_rates[48] - rev_rates[48]);
  //sp 9
  sp_rates[9] += (fwd_rates[48] - rev_rates[48]);
  //sp 11
  sp_rates[11] -= (fwd_rates[48] - rev_rates[48]);
  //sp 8
  sp_rates[8] += (fwd_rates[48] - rev_rates[48]);

  //rxn 49
  //sp 0
  sp_rates[0] -= (fwd_rates[49] - rev_rates[49]);
  //sp 9
  sp_rates[9] += (fwd_rates[49] - rev_rates[49]);
  //sp 11
  sp_rates[11] -= (fwd_rates[49] - rev_rates[49]);
  //sp 8
  sp_rates[8] += (fwd_rates[49] - rev_rates[49]);

  //rxn 50
  //sp 0
  sp_rates[0] -= (fwd_rates[50] - rev_rates[50]);
  //sp 11
  sp_rates[11] -= (fwd_rates[50] - rev_rates[50]);
  //sp 14
  sp_rates[14] = (fwd_rates[50] - rev_rates[50]);
  //sp 6
  sp_rates[6] += (fwd_rates[50] - rev_rates[50]);

  //rxn 51
  //sp 1
  sp_rates[1] += (fwd_rates[51] - rev_rates[51]);
  //sp 11
  sp_rates[11] -= (fwd_rates[51] - rev_rates[51]);
  //sp 20
  sp_rates[20] = (fwd_rates[51] - rev_rates[51]);
  //sp 21
  sp_rates[21] = -(fwd_rates[51] - rev_rates[51]);

  //rxn 52
  //sp 9
  sp_rates[9] += (fwd_rates[52] - rev_rates[52]);
  //sp 11
  sp_rates[11] -= (fwd_rates[52] - rev_rates[52]);
  //sp 20
  sp_rates[20] -= (fwd_rates[52] - rev_rates[52]);
  //sp 25
  sp_rates[25] = (fwd_rates[52] - rev_rates[52]);

  //rxn 53
  //sp 16
  sp_rates[16] += (fwd_rates[53] - rev_rates[53]);
  //sp 0
  sp_rates[0] += (fwd_rates[53] - rev_rates[53]);
  //sp 11
  sp_rates[11] -= (fwd_rates[53] - rev_rates[53]);
  //sp 20
  sp_rates[20] -= (fwd_rates[53] - rev_rates[53]);

  //rxn 54
  //sp 2
  sp_rates[2] += (fwd_rates[54] - rev_rates[54]);
  //sp 13
  sp_rates[13] += (fwd_rates[54] - rev_rates[54]);
  //sp 4
  sp_rates[4] -= (fwd_rates[54] - rev_rates[54]);
  //sp 12
  sp_rates[12] -= (fwd_rates[54] - rev_rates[54]);

  //rxn 55
  //sp 0
  sp_rates[0] += (fwd_rates[55] - rev_rates[55]);
  //sp 4
  sp_rates[4] += (fwd_rates[55] - rev_rates[55]);
  //sp 12
  sp_rates[12] -= (fwd_rates[55] - rev_rates[55]);
  //sp 5
  sp_rates[5] -= (fwd_rates[55] - rev_rates[55]);

  //rxn 56
  //sp 18
  sp_rates[18] += (fwd_rates[56] - rev_rates[56]);
  //sp 12
  sp_rates[12] -= (fwd_rates[56] - rev_rates[56]);
  //sp 4
  sp_rates[4] += (fwd_rates[56] - rev_rates[56]);
  //sp 6
  sp_rates[6] -= (fwd_rates[56] - rev_rates[56]);

  //rxn 57
  //sp 9
  sp_rates[9] += (fwd_rates[57] - rev_rates[57]);
  //sp 12
  sp_rates[12] -= (fwd_rates[57] - rev_rates[57]);
  //sp 13
  sp_rates[13] += (fwd_rates[57] - rev_rates[57]);
  //sp 6
  sp_rates[6] -= (fwd_rates[57] - rev_rates[57]);

  //rxn 58
  //sp 18
  sp_rates[18] += (fwd_rates[58] - rev_rates[58]);
  //sp 3
  sp_rates[3] -= (fwd_rates[58] - rev_rates[58]);
  //sp 12
  sp_rates[12] -= (fwd_rates[58] - rev_rates[58]);
  //sp 5
  sp_rates[5] += (fwd_rates[58] - rev_rates[58]);

  //rxn 59
  //sp 0
  sp_rates[0] += (fwd_rates[59] - rev_rates[59]);
  //sp 3
  sp_rates[3] -= (fwd_rates[59] - rev_rates[59]);
  //sp 12
  sp_rates[12] -= (fwd_rates[59] - rev_rates[59]);
  //sp 6
  sp_rates[6] += (fwd_rates[59] - rev_rates[59]);

  //rxn 60
  //sp 11
  sp_rates[11] += (fwd_rates[60] - rev_rates[60]);
  //sp 12
  sp_rates[12] -= 2.0 * (fwd_rates[60] - rev_rates[60]);
  //sp 13
  sp_rates[13] += (fwd_rates[60] - rev_rates[60]);

  //rxn 61
  //sp 8
  sp_rates[8] += (fwd_rates[61] - rev_rates[61]);
  //sp 4
  sp_rates[4] += (fwd_rates[61] - rev_rates[61]);
  //sp 12
  sp_rates[12] -= (fwd_rates[61] - rev_rates[61]);
  //sp 13
  sp_rates[13] -= (fwd_rates[61] - rev_rates[61]);

  //rxn 62
  //sp 0
  sp_rates[0] -= (fwd_rates[62] - rev_rates[62]);
  //sp 25
  sp_rates[25] += (fwd_rates[62] - rev_rates[62]);
  //sp 12
  sp_rates[12] -= (fwd_rates[62] - rev_rates[62]);
  //sp 4
  sp_rates[4] += (fwd_rates[62] - rev_rates[62]);

  //rxn 63
  //sp 0
  sp_rates[0] -= (fwd_rates[63] - rev_rates[63]);
  //sp 8
  sp_rates[8] += (fwd_rates[63] - rev_rates[63]);
  //sp 12
  sp_rates[12] -= (fwd_rates[63] - rev_rates[63]);
  //sp 6
  sp_rates[6] += (fwd_rates[63] - rev_rates[63]);

  //rxn 64
  //sp 20
  sp_rates[20] += (fwd_rates[64] - rev_rates[64]);
  //sp 11
  sp_rates[11] += (fwd_rates[64] - rev_rates[64]);
  //sp 12
  sp_rates[12] -= (fwd_rates[64] - rev_rates[64]);
  //sp 21
  sp_rates[21] -= (fwd_rates[64] - rev_rates[64]);

  //rxn 65
  //sp 25
  sp_rates[25] += (fwd_rates[65] - rev_rates[65]);
  //sp 12
  sp_rates[12] -= (fwd_rates[65] - rev_rates[65]);
  //sp 20
  sp_rates[20] -= (fwd_rates[65] - rev_rates[65]);
  //sp 6
  sp_rates[6] += (fwd_rates[65] - rev_rates[65]);

  //rxn 66
  //sp 0
  sp_rates[0] += (fwd_rates[66] - rev_rates[66]);
  //sp 18
  sp_rates[18] += (fwd_rates[66] - rev_rates[66]);
  //sp 12
  sp_rates[12] -= (fwd_rates[66] - rev_rates[66]);
  //sp 20
  sp_rates[20] -= (fwd_rates[66] - rev_rates[66]);

  //rxn 67
  //sp 0
  sp_rates[0] += (fwd_rates[67] - rev_rates[67]);
  //sp 4
  sp_rates[4] += (fwd_rates[67] - rev_rates[67]);
  //sp 13
  sp_rates[13] -= (fwd_rates[67] - rev_rates[67]);
  //sp 6
  sp_rates[6] -= (fwd_rates[67] - rev_rates[67]);

  //rxn 68
  //sp 0
  sp_rates[0] += (fwd_rates[68] - rev_rates[68]);
  //sp 5
  sp_rates[5] += (fwd_rates[68] - rev_rates[68]);
  //sp 3
  sp_rates[3] -= (fwd_rates[68] - rev_rates[68]);
  //sp 13
  sp_rates[13] -= (fwd_rates[68] - rev_rates[68]);

  //rxn 69
  //sp 0
  sp_rates[0] -= (fwd_rates[69] - rev_rates[69]);
  //sp 8
  sp_rates[8] += (fwd_rates[69] - rev_rates[69]);
  //sp 5
  sp_rates[5] += (fwd_rates[69] - rev_rates[69]);
  //sp 13
  sp_rates[13] -= (fwd_rates[69] - rev_rates[69]);

  //rxn 70
  //sp 8
  sp_rates[8] += (fwd_rates[70] - rev_rates[70]);
  //sp 4
  sp_rates[4] += (fwd_rates[70] - rev_rates[70]);
  //sp 14
  sp_rates[14] -= (fwd_rates[70] - rev_rates[70]);

  //rxn 71
  //sp 8
  sp_rates[8] += (fwd_rates[71] - rev_rates[71]);
  //sp 2
  sp_rates[2] += (fwd_rates[71] - rev_rates[71]);
  //sp 4
  sp_rates[4] -= (fwd_rates[71] - rev_rates[71]);
  //sp 14
  sp_rates[14] -= (fwd_rates[71] - rev_rates[71]);

  //rxn 72
  //sp 25
  sp_rates[25] += (fwd_rates[72] - rev_rates[72]);
  //sp 4
  sp_rates[4] += (fwd_rates[72] - rev_rates[72]);
  //sp 5
  sp_rates[5] -= (fwd_rates[72] - rev_rates[72]);
  //sp 14
  sp_rates[14] -= (fwd_rates[72] - rev_rates[72]);

  //rxn 73
  //sp 8
  sp_rates[8] += (fwd_rates[73] - rev_rates[73]);
  //sp 5
  sp_rates[5] -= (fwd_rates[73] - rev_rates[73]);
  //sp 14
  sp_rates[14] -= (fwd_rates[73] - rev_rates[73]);
  //sp 6
  sp_rates[6] += (fwd_rates[73] - rev_rates[73]);

  //rxn 74
  //sp 0
  sp_rates[0] += (fwd_rates[74] - rev_rates[74]);
  //sp 12
  sp_rates[12] += (fwd_rates[74] - rev_rates[74]);
  //sp 5
  sp_rates[5] -= (fwd_rates[74] - rev_rates[74]);
  //sp 14
  sp_rates[14] -= (fwd_rates[74] - rev_rates[74]);

  //rxn 75
  //sp 8
  sp_rates[8] += (fwd_rates[75] - rev_rates[75]);
  //sp 9
  sp_rates[9] += (fwd_rates[75] - rev_rates[75]);
  //sp 14
  sp_rates[14] -= (fwd_rates[75] - rev_rates[75]);
  //sp 6
  sp_rates[6] -= (fwd_rates[75] - rev_rates[75]);

  //rxn 76
  //sp 8
  sp_rates[8] += (fwd_rates[76] - rev_rates[76]);
  //sp 3
  sp_rates[3] -= (fwd_rates[76] - rev_rates[76]);
  //sp 14
  sp_rates[14] -= (fwd_rates[76] - rev_rates[76]);
  //sp 7
  sp_rates[7] += (fwd_rates[76] - rev_rates[76]);

  //rxn 77
  //sp 8
  sp_rates[8] += (fwd_rates[77] - rev_rates[77]);
  //sp 11
  sp_rates[11] += (fwd_rates[77] - rev_rates[77]);
  //sp 12
  sp_rates[12] -= (fwd_rates[77] - rev_rates[77]);
  //sp 14
  sp_rates[14] -= (fwd_rates[77] - rev_rates[77]);

  //rxn 78
  //sp 8
  sp_rates[8] += (fwd_rates[78] - rev_rates[78]);
  //sp 1
  sp_rates[1] += (fwd_rates[78] - rev_rates[78]);
  //sp 11
  sp_rates[11] -= (fwd_rates[78] - rev_rates[78]);
  //sp 14
  sp_rates[14] -= (fwd_rates[78] - rev_rates[78]);

  //rxn 79
  //sp 0
  sp_rates[0] -= (fwd_rates[79] - rev_rates[79]);
  //sp 8
  sp_rates[8] += (fwd_rates[79] - rev_rates[79]);
  //sp 18
  sp_rates[18] += (fwd_rates[79] - rev_rates[79]);
  //sp 14
  sp_rates[14] -= (fwd_rates[79] - rev_rates[79]);

  //rxn 80
  //sp 11
  sp_rates[11] += (fwd_rates[80] - rev_rates[80]) * pres_mod[9];
  //sp 6
  sp_rates[6] += (fwd_rates[80] - rev_rates[80]) * pres_mod[9];
  //sp 15
  sp_rates[15] = -(fwd_rates[80] - rev_rates[80]) * pres_mod[9];

  //rxn 81
  //sp 17
  sp_rates[17] = (fwd_rates[81] - rev_rates[81]);
  //sp 2
  sp_rates[2] += (fwd_rates[81] - rev_rates[81]);
  //sp 4
  sp_rates[4] -= (fwd_rates[81] - rev_rates[81]);
  //sp 15
  sp_rates[15] -= (fwd_rates[81] - rev_rates[81]);

  //rxn 82
  //sp 16
  sp_rates[16] += (fwd_rates[82] - rev_rates[82]);
  //sp 2
  sp_rates[2] += (fwd_rates[82] - rev_rates[82]);
  //sp 4
  sp_rates[4] -= (fwd_rates[82] - rev_rates[82]);
  //sp 15
  sp_rates[15] -= (fwd_rates[82] - rev_rates[82]);

  //rxn 83
  //sp 17
  sp_rates[17] += (fwd_rates[83] - rev_rates[83]);
  //sp 5
  sp_rates[5] -= (fwd_rates[83] - rev_rates[83]);
  //sp 6
  sp_rates[6] += (fwd_rates[83] - rev_rates[83]);
  //sp 15
  sp_rates[15] -= (fwd_rates[83] - rev_rates[83]);

  //rxn 84
  //sp 16
  sp_rates[16] += (fwd_rates[84] - rev_rates[84]);
  //sp 5
  sp_rates[5] -= (fwd_rates[84] - rev_rates[84]);
  //sp 6
  sp_rates[6] += (fwd_rates[84] - rev_rates[84]);
  //sp 15
  sp_rates[15] -= (fwd_rates[84] - rev_rates[84]);

  //rxn 85
  //sp 9
  sp_rates[9] += (fwd_rates[85] - rev_rates[85]);
  //sp 17
  sp_rates[17] += (fwd_rates[85] - rev_rates[85]);
  //sp 6
  sp_rates[6] -= (fwd_rates[85] - rev_rates[85]);
  //sp 15
  sp_rates[15] -= (fwd_rates[85] - rev_rates[85]);

  //rxn 86
  //sp 16
  sp_rates[16] += (fwd_rates[86] - rev_rates[86]);
  //sp 9
  sp_rates[9] += (fwd_rates[86] - rev_rates[86]);
  //sp 6
  sp_rates[6] -= (fwd_rates[86] - rev_rates[86]);
  //sp 15
  sp_rates[15] -= (fwd_rates[86] - rev_rates[86]);

  //rxn 87
  //sp 17
  sp_rates[17] += (fwd_rates[87] - rev_rates[87]);
  //sp 11
  sp_rates[11] -= (fwd_rates[87] - rev_rates[87]);
  //sp 1
  sp_rates[1] += (fwd_rates[87] - rev_rates[87]);
  //sp 15
  sp_rates[15] -= (fwd_rates[87] - rev_rates[87]);

  //rxn 88
  //sp 16
  sp_rates[16] += (fwd_rates[88] - rev_rates[88]);
  //sp 1
  sp_rates[1] += (fwd_rates[88] - rev_rates[88]);
  //sp 11
  sp_rates[11] -= (fwd_rates[88] - rev_rates[88]);
  //sp 15
  sp_rates[15] -= (fwd_rates[88] - rev_rates[88]);

  //rxn 89
  //sp 17
  sp_rates[17] += (fwd_rates[89] - rev_rates[89]);
  //sp 11
  sp_rates[11] += (fwd_rates[89] - rev_rates[89]);
  //sp 12
  sp_rates[12] -= (fwd_rates[89] - rev_rates[89]);
  //sp 15
  sp_rates[15] -= (fwd_rates[89] - rev_rates[89]);

  //rxn 90
  //sp 16
  sp_rates[16] += (fwd_rates[90] - rev_rates[90]);
  //sp 11
  sp_rates[11] += (fwd_rates[90] - rev_rates[90]);
  //sp 12
  sp_rates[12] -= (fwd_rates[90] - rev_rates[90]);
  //sp 15
  sp_rates[15] -= (fwd_rates[90] - rev_rates[90]);

  //rxn 91
  //sp 17
  sp_rates[17] += (fwd_rates[91] - rev_rates[91]);
  //sp 10
  sp_rates[10] += (fwd_rates[91] - rev_rates[91]);
  //sp 15
  sp_rates[15] -= (fwd_rates[91] - rev_rates[91]);
  //sp 7
  sp_rates[7] -= (fwd_rates[91] - rev_rates[91]);

  //rxn 92
  //sp 16
  sp_rates[16] += (fwd_rates[92] - rev_rates[92]);
  //sp 10
  sp_rates[10] += (fwd_rates[92] - rev_rates[92]);
  //sp 15
  sp_rates[15] -= (fwd_rates[92] - rev_rates[92]);
  //sp 7
  sp_rates[7] -= (fwd_rates[92] - rev_rates[92]);

  //rxn 93
  //sp 16
  sp_rates[16] -= (fwd_rates[93] - rev_rates[93]) * pres_mod[10];
  //sp 18
  sp_rates[18] += (fwd_rates[93] - rev_rates[93]) * pres_mod[10];
  //sp 4
  sp_rates[4] += (fwd_rates[93] - rev_rates[93]) * pres_mod[10];

  //rxn 94
  //sp 16
  sp_rates[16] -= (fwd_rates[94] - rev_rates[94]) * pres_mod[11];
  //sp 17
  sp_rates[17] += (fwd_rates[94] - rev_rates[94]) * pres_mod[11];

  //rxn 95
  //sp 16
  sp_rates[16] -= (fwd_rates[95] - rev_rates[95]) * pres_mod[12];
  //sp 0
  sp_rates[0] += (fwd_rates[95] - rev_rates[95]) * pres_mod[12];
  //sp 2
  sp_rates[2] += (fwd_rates[95] - rev_rates[95]) * pres_mod[12];

  //rxn 96
  //sp 16
  sp_rates[16] -= (fwd_rates[96] - rev_rates[96]);
  //sp 2
  sp_rates[2] += (fwd_rates[96] - rev_rates[96]);
  //sp 18
  sp_rates[18] += (fwd_rates[96] - rev_rates[96]);
  //sp 4
  sp_rates[4] -= (fwd_rates[96] - rev_rates[96]);

  //rxn 97
  //sp 16
  sp_rates[16] -= (fwd_rates[97] - rev_rates[97]);
  //sp 11
  sp_rates[11] += (fwd_rates[97] - rev_rates[97]);
  //sp 4
  sp_rates[4] -= (fwd_rates[97] - rev_rates[97]);
  //sp 6
  sp_rates[6] += (fwd_rates[97] - rev_rates[97]);

  //rxn 98
  //sp 16
  sp_rates[16] -= (fwd_rates[98] - rev_rates[98]);
  //sp 18
  sp_rates[18] += (fwd_rates[98] - rev_rates[98]);
  //sp 5
  sp_rates[5] -= (fwd_rates[98] - rev_rates[98]);
  //sp 6
  sp_rates[6] += (fwd_rates[98] - rev_rates[98]);

  //rxn 99
  //sp 16
  sp_rates[16] -= (fwd_rates[99] - rev_rates[99]);
  //sp 9
  sp_rates[9] += (fwd_rates[99] - rev_rates[99]);
  //sp 18
  sp_rates[18] += (fwd_rates[99] - rev_rates[99]);
  //sp 6
  sp_rates[6] -= (fwd_rates[99] - rev_rates[99]);

  //rxn 100
  //sp 16
  sp_rates[16] -= (fwd_rates[100] - rev_rates[100]);
  //sp 10
  sp_rates[10] += (fwd_rates[100] - rev_rates[100]);
  //sp 18
  sp_rates[18] += (fwd_rates[100] - rev_rates[100]);
  //sp 7
  sp_rates[7] -= (fwd_rates[100] - rev_rates[100]);

  //rxn 101
  //sp 16
  sp_rates[16] -= (fwd_rates[101] - rev_rates[101]);
  //sp 18
  sp_rates[18] += (fwd_rates[101] - rev_rates[101]);
  //sp 3
  sp_rates[3] -= (fwd_rates[101] - rev_rates[101]);
  //sp 7
  sp_rates[7] += (fwd_rates[101] - rev_rates[101]);

  //rxn 102
  //sp 16
  sp_rates[16] -= (fwd_rates[102] - rev_rates[102]);
  //sp 1
  sp_rates[1] += (fwd_rates[102] - rev_rates[102]);
  //sp 18
  sp_rates[18] += (fwd_rates[102] - rev_rates[102]);
  //sp 11
  sp_rates[11] -= (fwd_rates[102] - rev_rates[102]);

  //rxn 103
  //sp 16
  sp_rates[16] -= (fwd_rates[103] - rev_rates[103]);
  //sp 0
  sp_rates[0] -= (fwd_rates[103] - rev_rates[103]);
  //sp 18
  sp_rates[18] += 2.0 * (fwd_rates[103] - rev_rates[103]);

  //rxn 104
  //sp 16
  sp_rates[16] -= (fwd_rates[104] - rev_rates[104]);
  //sp 18
  sp_rates[18] += (fwd_rates[104] - rev_rates[104]);
  //sp 20
  sp_rates[20] -= (fwd_rates[104] - rev_rates[104]);
  //sp 21
  sp_rates[21] += (fwd_rates[104] - rev_rates[104]);

  //rxn 105
  //sp 17
  sp_rates[17] -= (fwd_rates[105] - rev_rates[105]) * pres_mod[13];
  //sp 18
  sp_rates[18] += (fwd_rates[105] - rev_rates[105]) * pres_mod[13];
  //sp 4
  sp_rates[4] += (fwd_rates[105] - rev_rates[105]) * pres_mod[13];

  //rxn 106
  //sp 17
  sp_rates[17] -= (fwd_rates[106] - rev_rates[106]);
  //sp 11
  sp_rates[11] += (fwd_rates[106] - rev_rates[106]);
  //sp 4
  sp_rates[4] -= (fwd_rates[106] - rev_rates[106]);
  //sp 6
  sp_rates[6] += (fwd_rates[106] - rev_rates[106]);

  //rxn 107
  //sp 17
  sp_rates[17] -= (fwd_rates[107] - rev_rates[107]);
  //sp 2
  sp_rates[2] += (fwd_rates[107] - rev_rates[107]);
  //sp 18
  sp_rates[18] += (fwd_rates[107] - rev_rates[107]);
  //sp 4
  sp_rates[4] -= (fwd_rates[107] - rev_rates[107]);

  //rxn 108
  //sp 17
  sp_rates[17] -= (fwd_rates[108] - rev_rates[108]);
  //sp 18
  sp_rates[18] += (fwd_rates[108] - rev_rates[108]);
  //sp 5
  sp_rates[5] -= (fwd_rates[108] - rev_rates[108]);
  //sp 6
  sp_rates[6] += (fwd_rates[108] - rev_rates[108]);

  //rxn 109
  //sp 17
  sp_rates[17] -= (fwd_rates[109] - rev_rates[109]);
  //sp 18
  sp_rates[18] += (fwd_rates[109] - rev_rates[109]);
  //sp 5
  sp_rates[5] -= (fwd_rates[109] - rev_rates[109]);
  //sp 6
  sp_rates[6] += (fwd_rates[109] - rev_rates[109]);

  //rxn 110
  //sp 17
  sp_rates[17] -= (fwd_rates[110] - rev_rates[110]);
  //sp 18
  sp_rates[18] += (fwd_rates[110] - rev_rates[110]);
  //sp 6
  sp_rates[6] -= (fwd_rates[110] - rev_rates[110]);
  //sp 9
  sp_rates[9] += (fwd_rates[110] - rev_rates[110]);

  //rxn 111
  //sp 17
  sp_rates[17] -= (fwd_rates[111] - rev_rates[111]);
  //sp 10
  sp_rates[10] += (fwd_rates[111] - rev_rates[111]);
  //sp 18
  sp_rates[18] += (fwd_rates[111] - rev_rates[111]);
  //sp 7
  sp_rates[7] -= (fwd_rates[111] - rev_rates[111]);

  //rxn 112
  //sp 17
  sp_rates[17] -= (fwd_rates[112] - rev_rates[112]);
  //sp 3
  sp_rates[3] += (fwd_rates[112] - rev_rates[112]);
  //sp 15
  sp_rates[15] += (fwd_rates[112] - rev_rates[112]);
  //sp 7
  sp_rates[7] -= (fwd_rates[112] - rev_rates[112]);

  //rxn 113
  //sp 17
  sp_rates[17] -= (fwd_rates[113] - rev_rates[113]);
  //sp 18
  sp_rates[18] += (fwd_rates[113] - rev_rates[113]);
  //sp 3
  sp_rates[3] -= (fwd_rates[113] - rev_rates[113]);
  //sp 7
  sp_rates[7] += (fwd_rates[113] - rev_rates[113]);

  //rxn 114
  //sp 17
  sp_rates[17] -= (fwd_rates[114] - rev_rates[114]);
  //sp 18
  sp_rates[18] += (fwd_rates[114] - rev_rates[114]);
  //sp 11
  sp_rates[11] -= (fwd_rates[114] - rev_rates[114]);
  //sp 1
  sp_rates[1] += (fwd_rates[114] - rev_rates[114]);

  //rxn 115
  //sp 17
  sp_rates[17] -= (fwd_rates[115] - rev_rates[115]);
  //sp 18
  sp_rates[18] += (fwd_rates[115] - rev_rates[115]);
  //sp 20
  sp_rates[20] -= (fwd_rates[115] - rev_rates[115]);
  //sp 21
  sp_rates[21] += (fwd_rates[115] - rev_rates[115]);

  //rxn 116
  //sp 0
  sp_rates[0] -= (fwd_rates[116] - rev_rates[116]) * pres_mod[14];
  //sp 18
  sp_rates[18] += (fwd_rates[116] - rev_rates[116]) * pres_mod[14];
  //sp 4
  sp_rates[4] -= (fwd_rates[116] - rev_rates[116]) * pres_mod[14];

  //rxn 117
  //sp 0
  sp_rates[0] += (fwd_rates[117] - rev_rates[117]);
  //sp 18
  sp_rates[18] -= (fwd_rates[117] - rev_rates[117]);
  //sp 2
  sp_rates[2] += (fwd_rates[117] - rev_rates[117]);
  //sp 4
  sp_rates[4] -= (fwd_rates[117] - rev_rates[117]);

  //rxn 118
  //sp 0
  sp_rates[0] += (fwd_rates[118] - rev_rates[118]);
  //sp 18
  sp_rates[18] -= (fwd_rates[118] - rev_rates[118]);
  //sp 5
  sp_rates[5] -= (fwd_rates[118] - rev_rates[118]);
  //sp 6
  sp_rates[6] += (fwd_rates[118] - rev_rates[118]);

  //rxn 119
  //sp 0
  sp_rates[0] += (fwd_rates[119] - rev_rates[119]);
  //sp 9
  sp_rates[9] += (fwd_rates[119] - rev_rates[119]);
  //sp 18
  sp_rates[18] -= (fwd_rates[119] - rev_rates[119]);
  //sp 6
  sp_rates[6] -= (fwd_rates[119] - rev_rates[119]);

  //rxn 120
  //sp 0
  sp_rates[0] += (fwd_rates[120] - rev_rates[120]);
  //sp 18
  sp_rates[18] -= (fwd_rates[120] - rev_rates[120]);
  //sp 3
  sp_rates[3] -= (fwd_rates[120] - rev_rates[120]);
  //sp 7
  sp_rates[7] += (fwd_rates[120] - rev_rates[120]);

  //rxn 121
  //sp 9
  sp_rates[9] += (fwd_rates[121] - rev_rates[121]);
  //sp 18
  sp_rates[18] -= 2.0 * (fwd_rates[121] - rev_rates[121]);
  //sp 25
  sp_rates[25] += (fwd_rates[121] - rev_rates[121]);

  //rxn 122
  //sp 0
  sp_rates[0] += (fwd_rates[122] - rev_rates[122]);
  //sp 18
  sp_rates[18] -= (fwd_rates[122] - rev_rates[122]);
  //sp 20
  sp_rates[20] -= (fwd_rates[122] - rev_rates[122]);
  //sp 21
  sp_rates[21] += (fwd_rates[122] - rev_rates[122]);

  //rxn 123
  //sp 0
  sp_rates[0] -= (fwd_rates[123] - rev_rates[123]);
  //sp 20
  sp_rates[20] += (fwd_rates[123] - rev_rates[123]);
  //sp 6
  sp_rates[6] += (fwd_rates[123] - rev_rates[123]);
  //sp 7
  sp_rates[7] -= (fwd_rates[123] - rev_rates[123]);

  //rxn 124
  //sp 0
  sp_rates[0] -= (fwd_rates[124] - rev_rates[124]) * pres_mod[15];
  //sp 20
  sp_rates[20] += (fwd_rates[124] - rev_rates[124]) * pres_mod[15];
  //sp 5
  sp_rates[5] -= (fwd_rates[124] - rev_rates[124]) * pres_mod[15];

  //rxn 125
  //sp 0
  sp_rates[0] -= (fwd_rates[125] - rev_rates[125]) * pres_mod[16];
  //sp 24
  sp_rates[24] = (fwd_rates[125] - rev_rates[125]) * pres_mod[16];
  //sp 7
  sp_rates[7] -= (fwd_rates[125] - rev_rates[125]) * pres_mod[16];

  //rxn 126
  //sp 0
  sp_rates[0] += (fwd_rates[126] - rev_rates[126]);
  //sp 4
  sp_rates[4] -= (fwd_rates[126] - rev_rates[126]);
  //sp 20
  sp_rates[20] -= (fwd_rates[126] - rev_rates[126]);
  //sp 6
  sp_rates[6] += (fwd_rates[126] - rev_rates[126]);

  //rxn 127
  //sp 0
  sp_rates[0] += (fwd_rates[127] - rev_rates[127]);
  //sp 3
  sp_rates[3] += (fwd_rates[127] - rev_rates[127]);
  //sp 20
  sp_rates[20] -= (fwd_rates[127] - rev_rates[127]);
  //sp 5
  sp_rates[5] -= (fwd_rates[127] - rev_rates[127]);

  //rxn 128
  //sp 3
  sp_rates[3] += (fwd_rates[128] - rev_rates[128]);
  //sp 20
  sp_rates[20] -= (fwd_rates[128] - rev_rates[128]);
  //sp 21
  sp_rates[21] += (fwd_rates[128] - rev_rates[128]);
  //sp 7
  sp_rates[7] -= (fwd_rates[128] - rev_rates[128]);

  //rxn 129
  //sp 3
  sp_rates[3] += (fwd_rates[129] - rev_rates[129]);
  //sp 20
  sp_rates[20] -= (fwd_rates[129] - rev_rates[129]);
  //sp 22
  sp_rates[22] = (fwd_rates[129] - rev_rates[129]);
  //sp 7
  sp_rates[7] -= (fwd_rates[129] - rev_rates[129]);

  //rxn 130
  //sp 0
  sp_rates[0] += 2.0 * (fwd_rates[130] - rev_rates[130]);
  //sp 3
  sp_rates[3] += (fwd_rates[130] - rev_rates[130]);
  //sp 20
  sp_rates[20] -= 2.0 * (fwd_rates[130] - rev_rates[130]);

  //rxn 131
  //sp 0
  sp_rates[0] += (fwd_rates[131] - rev_rates[131]);
  //sp 20
  sp_rates[20] -= 2.0 * (fwd_rates[131] - rev_rates[131]);
  //sp 23
  sp_rates[23] = (fwd_rates[131] - rev_rates[131]);

  //rxn 132
  //sp 0
  sp_rates[0] -= (fwd_rates[132] - rev_rates[132]);
  //sp 25
  sp_rates[25] += (fwd_rates[132] - rev_rates[132]);
  //sp 3
  sp_rates[3] += (fwd_rates[132] - rev_rates[132]);
  //sp 20
  sp_rates[20] -= (fwd_rates[132] - rev_rates[132]);

  //rxn 133
  //sp 0
  sp_rates[0] -= (fwd_rates[133] - rev_rates[133]) * pres_mod[17];
  //sp 21
  sp_rates[21] += (fwd_rates[133] - rev_rates[133]) * pres_mod[17];
  //sp 6
  sp_rates[6] -= (fwd_rates[133] - rev_rates[133]) * pres_mod[17];

  //rxn 134
  //sp 2
  sp_rates[2] -= (fwd_rates[134] - rev_rates[134]);
  //sp 21
  sp_rates[21] += (fwd_rates[134] - rev_rates[134]);
  //sp 20
  sp_rates[20] -= (fwd_rates[134] - rev_rates[134]);
  //sp 4
  sp_rates[4] += (fwd_rates[134] - rev_rates[134]);

  //rxn 135
  //sp 25
  sp_rates[25] += (fwd_rates[135] - rev_rates[135]);
  //sp 5
  sp_rates[5] += (fwd_rates[135] - rev_rates[135]);
  //sp 20
  sp_rates[20] -= (fwd_rates[135] - rev_rates[135]);
  //sp 13
  sp_rates[13] -= (fwd_rates[135] - rev_rates[135]);

  //rxn 136
  //sp 18
  sp_rates[18] += (fwd_rates[136] - rev_rates[136]);
  //sp 4
  sp_rates[4] -= (fwd_rates[136] - rev_rates[136]);
  //sp 21
  sp_rates[21] -= (fwd_rates[136] - rev_rates[136]);
  //sp 6
  sp_rates[6] += (fwd_rates[136] - rev_rates[136]);

  //rxn 137
  //sp 0
  sp_rates[0] += (fwd_rates[137] - rev_rates[137]);
  //sp 9
  sp_rates[9] += (fwd_rates[137] - rev_rates[137]);
  //sp 4
  sp_rates[4] -= (fwd_rates[137] - rev_rates[137]);
  //sp 21
  sp_rates[21] -= (fwd_rates[137] - rev_rates[137]);

  //rxn 138
  //sp 5
  sp_rates[5] -= (fwd_rates[138] - rev_rates[138]);
  //sp 20
  sp_rates[20] += (fwd_rates[138] - rev_rates[138]);
  //sp 21
  sp_rates[21] -= (fwd_rates[138] - rev_rates[138]);
  //sp 6
  sp_rates[6] += (fwd_rates[138] - rev_rates[138]);

  //rxn 139
  //sp 9
  sp_rates[9] += (fwd_rates[139] - rev_rates[139]);
  //sp 20
  sp_rates[20] += (fwd_rates[139] - rev_rates[139]);
  //sp 21
  sp_rates[21] -= (fwd_rates[139] - rev_rates[139]);
  //sp 6
  sp_rates[6] -= (fwd_rates[139] - rev_rates[139]);

  //rxn 140
  //sp 24
  sp_rates[24] += (fwd_rates[140] - rev_rates[140]);
  //sp 0
  sp_rates[0] += (fwd_rates[140] - rev_rates[140]);
  //sp 20
  sp_rates[20] -= (fwd_rates[140] - rev_rates[140]);
  //sp 21
  sp_rates[21] -= (fwd_rates[140] - rev_rates[140]);

  //rxn 141
  //sp 0
  sp_rates[0] += (fwd_rates[141] - rev_rates[141]);
  //sp 9
  sp_rates[9] += (fwd_rates[141] - rev_rates[141]);
  //sp 20
  sp_rates[20] += (fwd_rates[141] - rev_rates[141]);
  //sp 21
  sp_rates[21] -= 2.0 * (fwd_rates[141] - rev_rates[141]);

  //rxn 142
  //sp 21
  sp_rates[21] += (fwd_rates[142] - rev_rates[142]) * pres_mod[18];
  //sp 22
  sp_rates[22] -= (fwd_rates[142] - rev_rates[142]) * pres_mod[18];

  //rxn 143
  //sp 2
  sp_rates[2] += (fwd_rates[143] - rev_rates[143]);
  //sp 4
  sp_rates[4] -= (fwd_rates[143] - rev_rates[143]);
  //sp 20
  sp_rates[20] += (fwd_rates[143] - rev_rates[143]);
  //sp 22
  sp_rates[22] -= (fwd_rates[143] - rev_rates[143]);

  //rxn 144
  //sp 20
  sp_rates[20] += (fwd_rates[144] - rev_rates[144]);
  //sp 5
  sp_rates[5] -= (fwd_rates[144] - rev_rates[144]);
  //sp 22
  sp_rates[22] -= (fwd_rates[144] - rev_rates[144]);
  //sp 6
  sp_rates[6] += (fwd_rates[144] - rev_rates[144]);

  //rxn 145
  //sp 9
  sp_rates[9] += (fwd_rates[145] - rev_rates[145]);
  //sp 20
  sp_rates[20] += (fwd_rates[145] - rev_rates[145]);
  //sp 22
  sp_rates[22] -= (fwd_rates[145] - rev_rates[145]);
  //sp 6
  sp_rates[6] -= (fwd_rates[145] - rev_rates[145]);

  //rxn 146
  //sp 20
  sp_rates[20] -= (fwd_rates[146] - rev_rates[146]) * pres_mod[19];
  //sp 5
  sp_rates[5] -= (fwd_rates[146] - rev_rates[146]) * pres_mod[19];
  //sp 23
  sp_rates[23] += (fwd_rates[146] - rev_rates[146]) * pres_mod[19];

  //rxn 147
  //sp 4
  sp_rates[4] -= (fwd_rates[147] - rev_rates[147]);
  //sp 20
  sp_rates[20] += (fwd_rates[147] - rev_rates[147]);
  //sp 6
  sp_rates[6] += (fwd_rates[147] - rev_rates[147]);
  //sp 23
  sp_rates[23] -= (fwd_rates[147] - rev_rates[147]);

  //rxn 148
  //sp 3
  sp_rates[3] += (fwd_rates[148] - rev_rates[148]);
  //sp 20
  sp_rates[20] += (fwd_rates[148] - rev_rates[148]);
  //sp 5
  sp_rates[5] -= (fwd_rates[148] - rev_rates[148]);
  //sp 23
  sp_rates[23] -= (fwd_rates[148] - rev_rates[148]);

  //rxn 149
  //sp 7
  sp_rates[7] += (fwd_rates[149] - rev_rates[149]);
  //sp 20
  sp_rates[20] += (fwd_rates[149] - rev_rates[149]);
  //sp 6
  sp_rates[6] -= (fwd_rates[149] - rev_rates[149]);
  //sp 23
  sp_rates[23] -= (fwd_rates[149] - rev_rates[149]);

  //rxn 150
  //sp 3
  sp_rates[3] += (fwd_rates[150] - rev_rates[150]);
  //sp 6
  sp_rates[6] += (fwd_rates[150] - rev_rates[150]);
  //sp 7
  sp_rates[7] -= (fwd_rates[150] - rev_rates[150]);
  //sp 20
  sp_rates[20] += (fwd_rates[150] - rev_rates[150]);
  //sp 23
  sp_rates[23] -= (fwd_rates[150] - rev_rates[150]);

  //rxn 151
  //sp 0
  sp_rates[0] += (fwd_rates[151] - rev_rates[151]) * pres_mod[20];
  //sp 3
  sp_rates[3] += (fwd_rates[151] - rev_rates[151]) * pres_mod[20];
  //sp 23
  sp_rates[23] -= (fwd_rates[151] - rev_rates[151]) * pres_mod[20];

  //rxn 152
  //sp 24
  sp_rates[24] += (fwd_rates[152] - rev_rates[152]) * pres_mod[21];
  //sp 20
  sp_rates[20] -= (fwd_rates[152] - rev_rates[152]) * pres_mod[21];
  //sp 6
  sp_rates[6] -= (fwd_rates[152] - rev_rates[152]) * pres_mod[21];

  //rxn 153
  //sp 24
  sp_rates[24] -= (fwd_rates[153] - rev_rates[153]);
  //sp 2
  sp_rates[2] += (fwd_rates[153] - rev_rates[153]);
  //sp 4
  sp_rates[4] -= (fwd_rates[153] - rev_rates[153]);
  //sp 23
  sp_rates[23] += (fwd_rates[153] - rev_rates[153]);

  //rxn 154
  //sp 24
  sp_rates[24] -= (fwd_rates[154] - rev_rates[154]);
  //sp 9
  sp_rates[9] += (fwd_rates[154] - rev_rates[154]);
  //sp 4
  sp_rates[4] -= (fwd_rates[154] - rev_rates[154]);
  //sp 20
  sp_rates[20] += (fwd_rates[154] - rev_rates[154]);

  //rxn 155
  //sp 24
  sp_rates[24] -= (fwd_rates[155] - rev_rates[155]);
  //sp 4
  sp_rates[4] -= (fwd_rates[155] - rev_rates[155]);
  //sp 21
  sp_rates[21] += (fwd_rates[155] - rev_rates[155]);
  //sp 6
  sp_rates[6] += (fwd_rates[155] - rev_rates[155]);

  //rxn 156
  //sp 24
  sp_rates[24] -= (fwd_rates[156] - rev_rates[156]);
  //sp 9
  sp_rates[9] += (fwd_rates[156] - rev_rates[156]);
  //sp 6
  sp_rates[6] -= (fwd_rates[156] - rev_rates[156]);
  //sp 23
  sp_rates[23] += (fwd_rates[156] - rev_rates[156]);

  //rxn 157
  //sp 8
  sp_rates[8] += (fwd_rates[157] - rev_rates[157]) * pres_mod[22];
  //sp 25
  sp_rates[25] -= (fwd_rates[157] - rev_rates[157]) * pres_mod[22];
  //sp 5
  sp_rates[5] += (fwd_rates[157] - rev_rates[157]) * pres_mod[22];

  //rxn 158
  //sp 8
  sp_rates[8] += (fwd_rates[158] - rev_rates[158]);
  //sp 25
  sp_rates[25] -= (fwd_rates[158] - rev_rates[158]);
  //sp 4
  sp_rates[4] -= (fwd_rates[158] - rev_rates[158]);
  //sp 6
  sp_rates[6] += (fwd_rates[158] - rev_rates[158]);

  //rxn 159
  //sp 8
  sp_rates[8] += (fwd_rates[159] - rev_rates[159]);
  //sp 25
  sp_rates[25] -= (fwd_rates[159] - rev_rates[159]);
  //sp 4
  sp_rates[4] -= (fwd_rates[159] - rev_rates[159]);
  //sp 6
  sp_rates[6] += (fwd_rates[159] - rev_rates[159]);

  //rxn 160
  //sp 0
  sp_rates[0] += 2.0 * (fwd_rates[160] - rev_rates[160]);
  //sp 25
  sp_rates[25] -= (fwd_rates[160] - rev_rates[160]);
  //sp 5
  sp_rates[5] -= (fwd_rates[160] - rev_rates[160]);

  //rxn 161
  //sp 8
  sp_rates[8] += (fwd_rates[161] - rev_rates[161]);
  //sp 25
  sp_rates[25] -= (fwd_rates[161] - rev_rates[161]);
  //sp 3
  sp_rates[3] += (fwd_rates[161] - rev_rates[161]);
  //sp 5
  sp_rates[5] -= (fwd_rates[161] - rev_rates[161]);

  //rxn 162
  //sp 8
  sp_rates[8] += (fwd_rates[162] - rev_rates[162]);
  //sp 25
  sp_rates[25] -= (fwd_rates[162] - rev_rates[162]);
  //sp 6
  sp_rates[6] -= (fwd_rates[162] - rev_rates[162]);
  //sp 7
  sp_rates[7] += (fwd_rates[162] - rev_rates[162]);

  //rxn 163
  //sp 0
  sp_rates[0] -= (fwd_rates[163] - rev_rates[163]);
  //sp 25
  sp_rates[25] += (fwd_rates[163] - rev_rates[163]);
  //sp 18
  sp_rates[18] -= (fwd_rates[163] - rev_rates[163]);
  //sp 6
  sp_rates[6] += (fwd_rates[163] - rev_rates[163]);

  //rxn 164
  //sp 0
  sp_rates[0] += (fwd_rates[164] - rev_rates[164]);
  //sp 18
  sp_rates[18] -= (fwd_rates[164] - rev_rates[164]);
  //sp 12
  sp_rates[12] += (fwd_rates[164] - rev_rates[164]);
  //sp 13
  sp_rates[13] -= (fwd_rates[164] - rev_rates[164]);

  //rxn 165
  //sp 0
  sp_rates[0] -= (fwd_rates[165] - rev_rates[165]);
  //sp 25
  sp_rates[25] -= (fwd_rates[165] - rev_rates[165]);
  //sp 20
  sp_rates[20] += (fwd_rates[165] - rev_rates[165]);
  //sp 8
  sp_rates[8] += (fwd_rates[165] - rev_rates[165]);

  //rxn 166
  //sp 8
  sp_rates[8] += (fwd_rates[166] - rev_rates[166]);
  //sp 25
  sp_rates[25] -= (fwd_rates[166] - rev_rates[166]);
  //sp 13
  sp_rates[13] -= (fwd_rates[166] - rev_rates[166]);
  //sp 0
  sp_rates[0] += (fwd_rates[166] - rev_rates[166]);

  //rxn 167
  //sp 2
  sp_rates[2] += (fwd_rates[167] - rev_rates[167]);
  //sp 11
  sp_rates[11] -= 2.0 * (fwd_rates[167] - rev_rates[167]);
  //sp 28
  sp_rates[28] = (fwd_rates[167] - rev_rates[167]);

  //rxn 168
  //sp 2
  sp_rates[2] += (fwd_rates[168] - rev_rates[168]);
  //sp 11
  sp_rates[11] -= 2.0 * (fwd_rates[168] - rev_rates[168]);
  //sp 29
  sp_rates[29] = (fwd_rates[168] - rev_rates[168]);

  //rxn 169
  //sp 28
  sp_rates[28] += (fwd_rates[169] - rev_rates[169]);
  //sp 11
  sp_rates[11] -= (fwd_rates[169] - rev_rates[169]);
  //sp 12
  sp_rates[12] -= (fwd_rates[169] - rev_rates[169]);
  //sp 4
  sp_rates[4] += (fwd_rates[169] - rev_rates[169]);

  //rxn 170
  //sp 27
  sp_rates[27] = (fwd_rates[170] - rev_rates[170]);
  //sp 17
  sp_rates[17] -= (fwd_rates[170] - rev_rates[170]);
  //sp 11
  sp_rates[11] -= (fwd_rates[170] - rev_rates[170]);
  //sp 6
  sp_rates[6] += (fwd_rates[170] - rev_rates[170]);

  //rxn 171
  //sp 17
  sp_rates[17] -= (fwd_rates[171] - rev_rates[171]);
  //sp 11
  sp_rates[11] -= (fwd_rates[171] - rev_rates[171]);
  //sp 29
  sp_rates[29] += (fwd_rates[171] - rev_rates[171]);
  //sp 9
  sp_rates[9] += (fwd_rates[171] - rev_rates[171]);

  //rxn 172
  //sp 26
  sp_rates[26] = (fwd_rates[172] - rev_rates[172]) * pres_mod[23];
  //sp 11
  sp_rates[11] -= 2.0 * (fwd_rates[172] - rev_rates[172]) * pres_mod[23];

  //rxn 173
  //sp 27
  sp_rates[27] += (fwd_rates[173] - rev_rates[173]);
  //sp 26
  sp_rates[26] -= (fwd_rates[173] - rev_rates[173]);
  //sp 2
  sp_rates[2] += (fwd_rates[173] - rev_rates[173]);
  //sp 4
  sp_rates[4] -= (fwd_rates[173] - rev_rates[173]);

  //rxn 174
  //sp 26
  sp_rates[26] -= (fwd_rates[174] - rev_rates[174]);
  //sp 12
  sp_rates[12] += (fwd_rates[174] - rev_rates[174]);
  //sp 5
  sp_rates[5] -= (fwd_rates[174] - rev_rates[174]);
  //sp 15
  sp_rates[15] += (fwd_rates[174] - rev_rates[174]);

  //rxn 175
  //sp 26
  sp_rates[26] -= (fwd_rates[175] - rev_rates[175]);
  //sp 27
  sp_rates[27] += (fwd_rates[175] - rev_rates[175]);
  //sp 5
  sp_rates[5] -= (fwd_rates[175] - rev_rates[175]);
  //sp 6
  sp_rates[6] += (fwd_rates[175] - rev_rates[175]);

  //rxn 176
  //sp 9
  sp_rates[9] += (fwd_rates[176] - rev_rates[176]);
  //sp 26
  sp_rates[26] -= (fwd_rates[176] - rev_rates[176]);
  //sp 27
  sp_rates[27] += (fwd_rates[176] - rev_rates[176]);
  //sp 6
  sp_rates[6] -= (fwd_rates[176] - rev_rates[176]);

  //rxn 177
  //sp 27
  sp_rates[27] += (fwd_rates[177] - rev_rates[177]);
  //sp 1
  sp_rates[1] += (fwd_rates[177] - rev_rates[177]);
  //sp 26
  sp_rates[26] -= (fwd_rates[177] - rev_rates[177]);
  //sp 11
  sp_rates[11] -= (fwd_rates[177] - rev_rates[177]);

  //rxn 178
  //sp 27
  sp_rates[27] -= (fwd_rates[178] - rev_rates[178]);
  //sp 4
  sp_rates[4] += (fwd_rates[178] - rev_rates[178]);
  //sp 28
  sp_rates[28] += (fwd_rates[178] - rev_rates[178]);

  //rxn 179
  //sp 2
  sp_rates[2] += (fwd_rates[179] - rev_rates[179]);
  //sp 27
  sp_rates[27] -= (fwd_rates[179] - rev_rates[179]);
  //sp 4
  sp_rates[4] -= (fwd_rates[179] - rev_rates[179]);
  //sp 28
  sp_rates[28] += (fwd_rates[179] - rev_rates[179]);

  //rxn 180
  //sp 27
  sp_rates[27] -= (fwd_rates[180] - rev_rates[180]);
  //sp 28
  sp_rates[28] += (fwd_rates[180] - rev_rates[180]);
  //sp 5
  sp_rates[5] -= (fwd_rates[180] - rev_rates[180]);
  //sp 6
  sp_rates[6] += (fwd_rates[180] - rev_rates[180]);

  //rxn 181
  //sp 11
  sp_rates[11] += (fwd_rates[181] - rev_rates[181]);
  //sp 18
  sp_rates[18] += (fwd_rates[181] - rev_rates[181]);
  //sp 27
  sp_rates[27] -= (fwd_rates[181] - rev_rates[181]);
  //sp 5
  sp_rates[5] -= (fwd_rates[181] - rev_rates[181]);

  //rxn 182
  //sp 0
  sp_rates[0] += fwd_rates[182];
  //sp 4
  sp_rates[4] += fwd_rates[182];
  //sp 5
  sp_rates[5] -= fwd_rates[182];
  //sp 11
  sp_rates[11] += fwd_rates[182];
  //sp 27
  sp_rates[27] -= fwd_rates[182];

  //rxn 183
  //sp 9
  sp_rates[9] += (fwd_rates[183] - rev_rates[182]);
  //sp 27
  sp_rates[27] -= (fwd_rates[183] - rev_rates[182]);
  //sp 28
  sp_rates[28] += (fwd_rates[183] - rev_rates[182]);
  //sp 6
  sp_rates[6] -= (fwd_rates[183] - rev_rates[182]);

  //rxn 184
  //sp 9
  sp_rates[9] += (fwd_rates[184] - rev_rates[183]);
  //sp 27
  sp_rates[27] -= (fwd_rates[184] - rev_rates[183]);
  //sp 29
  sp_rates[29] += (fwd_rates[184] - rev_rates[183]);
  //sp 6
  sp_rates[6] -= (fwd_rates[184] - rev_rates[183]);

  //rxn 185
  //sp 1
  sp_rates[1] += (fwd_rates[185] - rev_rates[184]);
  //sp 18
  sp_rates[18] += (fwd_rates[185] - rev_rates[184]);
  //sp 27
  sp_rates[27] -= (fwd_rates[185] - rev_rates[184]);
  //sp 6
  sp_rates[6] -= (fwd_rates[185] - rev_rates[184]);

  //rxn 186
  //sp 10
  sp_rates[10] += (fwd_rates[186] - rev_rates[185]);
  //sp 27
  sp_rates[27] -= (fwd_rates[186] - rev_rates[185]);
  //sp 28
  sp_rates[28] += (fwd_rates[186] - rev_rates[185]);
  //sp 7
  sp_rates[7] -= (fwd_rates[186] - rev_rates[185]);

  //rxn 187
  //sp 3
  sp_rates[3] += (fwd_rates[187] - rev_rates[186]);
  //sp 26
  sp_rates[26] += (fwd_rates[187] - rev_rates[186]);
  //sp 27
  sp_rates[27] -= (fwd_rates[187] - rev_rates[186]);
  //sp 7
  sp_rates[7] -= (fwd_rates[187] - rev_rates[186]);

  //rxn 188
  //sp 11
  sp_rates[11] -= (fwd_rates[188] - rev_rates[187]);
  //sp 1
  sp_rates[1] += (fwd_rates[188] - rev_rates[187]);
  //sp 27
  sp_rates[27] -= (fwd_rates[188] - rev_rates[187]);
  //sp 28
  sp_rates[28] += (fwd_rates[188] - rev_rates[187]);

  //rxn 189
  //sp 11
  sp_rates[11] -= (fwd_rates[189] - rev_rates[188]);
  //sp 1
  sp_rates[1] += (fwd_rates[189] - rev_rates[188]);
  //sp 27
  sp_rates[27] -= (fwd_rates[189] - rev_rates[188]);
  //sp 29
  sp_rates[29] += (fwd_rates[189] - rev_rates[188]);

  //rxn 190
  //sp 11
  sp_rates[11] += (fwd_rates[190] - rev_rates[189]);
  //sp 27
  sp_rates[27] -= (fwd_rates[190] - rev_rates[189]);
  //sp 12
  sp_rates[12] -= (fwd_rates[190] - rev_rates[189]);
  //sp 28
  sp_rates[28] += (fwd_rates[190] - rev_rates[189]);

  //rxn 191
  //sp 28
  sp_rates[28] -= (fwd_rates[191] - rev_rates[190]);
  //sp 4
  sp_rates[4] += (fwd_rates[191] - rev_rates[190]);
  //sp 14
  sp_rates[14] += (fwd_rates[191] - rev_rates[190]);

  //rxn 192
  //sp 2
  sp_rates[2] += (fwd_rates[192] - rev_rates[191]);
  //sp 4
  sp_rates[4] -= (fwd_rates[192] - rev_rates[191]);
  //sp 28
  sp_rates[28] -= (fwd_rates[192] - rev_rates[191]);
  //sp 14
  sp_rates[14] += (fwd_rates[192] - rev_rates[191]);

  //rxn 193
  //sp 28
  sp_rates[28] -= (fwd_rates[193] - rev_rates[192]);
  //sp 5
  sp_rates[5] -= (fwd_rates[193] - rev_rates[192]);
  //sp 14
  sp_rates[14] += (fwd_rates[193] - rev_rates[192]);
  //sp 6
  sp_rates[6] += (fwd_rates[193] - rev_rates[192]);

  //rxn 194
  //sp 0
  sp_rates[0] += (fwd_rates[194] - rev_rates[193]);
  //sp 11
  sp_rates[11] += (fwd_rates[194] - rev_rates[193]);
  //sp 28
  sp_rates[28] -= (fwd_rates[194] - rev_rates[193]);
  //sp 5
  sp_rates[5] -= (fwd_rates[194] - rev_rates[193]);

  //rxn 195
  //sp 9
  sp_rates[9] += (fwd_rates[195] - rev_rates[194]);
  //sp 28
  sp_rates[28] -= (fwd_rates[195] - rev_rates[194]);
  //sp 6
  sp_rates[6] -= (fwd_rates[195] - rev_rates[194]);
  //sp 14
  sp_rates[14] += (fwd_rates[195] - rev_rates[194]);

  //rxn 196
  //sp 1
  sp_rates[1] += (fwd_rates[196] - rev_rates[195]);
  //sp 11
  sp_rates[11] -= (fwd_rates[196] - rev_rates[195]);
  //sp 28
  sp_rates[28] -= (fwd_rates[196] - rev_rates[195]);
  //sp 14
  sp_rates[14] += (fwd_rates[196] - rev_rates[195]);

  //rxn 197
  //sp 11
  sp_rates[11] += (fwd_rates[197] - rev_rates[196]);
  //sp 28
  sp_rates[28] -= (fwd_rates[197] - rev_rates[196]);
  //sp 12
  sp_rates[12] -= (fwd_rates[197] - rev_rates[196]);
  //sp 14
  sp_rates[14] += (fwd_rates[197] - rev_rates[196]);

  //rxn 198
  //sp 0
  sp_rates[0] -= (fwd_rates[198] - rev_rates[197]);
  //sp 25
  sp_rates[25] += (fwd_rates[198] - rev_rates[197]);
  //sp 11
  sp_rates[11] += (fwd_rates[198] - rev_rates[197]);
  //sp 28
  sp_rates[28] -= (fwd_rates[198] - rev_rates[197]);

  //rxn 199
  //sp 4
  sp_rates[4] += (fwd_rates[199] - rev_rates[198]);
  //sp 29
  sp_rates[29] -= (fwd_rates[199] - rev_rates[198]);
  //sp 14
  sp_rates[14] += (fwd_rates[199] - rev_rates[198]);

  //rxn 200
  //sp 4
  sp_rates[4] += (fwd_rates[200] - rev_rates[199]);
  //sp 29
  sp_rates[29] -= (fwd_rates[200] - rev_rates[199]);
  //sp 14
  sp_rates[14] += (fwd_rates[200] - rev_rates[199]);

  //rxn 201
  //sp 2
  sp_rates[2] += (fwd_rates[201] - rev_rates[200]);
  //sp 4
  sp_rates[4] -= (fwd_rates[201] - rev_rates[200]);
  //sp 29
  sp_rates[29] -= (fwd_rates[201] - rev_rates[200]);
  //sp 14
  sp_rates[14] += (fwd_rates[201] - rev_rates[200]);

  //rxn 202
  //sp 28
  sp_rates[28] += (fwd_rates[202] - rev_rates[201]);
  //sp 29
  sp_rates[29] -= (fwd_rates[202] - rev_rates[201]);

  //rxn 203
  //sp 5
  sp_rates[5] -= (fwd_rates[203] - rev_rates[202]);
  //sp 29
  sp_rates[29] -= (fwd_rates[203] - rev_rates[202]);
  //sp 14
  sp_rates[14] += (fwd_rates[203] - rev_rates[202]);
  //sp 6
  sp_rates[6] += (fwd_rates[203] - rev_rates[202]);

  //rxn 204
  //sp 0
  sp_rates[0] += (fwd_rates[204] - rev_rates[203]);
  //sp 5
  sp_rates[5] -= (fwd_rates[204] - rev_rates[203]);
  //sp 11
  sp_rates[11] += (fwd_rates[204] - rev_rates[203]);
  //sp 29
  sp_rates[29] -= (fwd_rates[204] - rev_rates[203]);

  //rxn 205
  //sp 9
  sp_rates[9] += (fwd_rates[205] - rev_rates[204]);
  //sp 29
  sp_rates[29] -= (fwd_rates[205] - rev_rates[204]);
  //sp 6
  sp_rates[6] -= (fwd_rates[205] - rev_rates[204]);
  //sp 14
  sp_rates[14] += (fwd_rates[205] - rev_rates[204]);

  //rxn 206
  //sp 0
  sp_rates[0] += fwd_rates[206];
  //sp 4
  sp_rates[4] += fwd_rates[206];
  //sp 6
  sp_rates[6] -= fwd_rates[206];
  //sp 11
  sp_rates[11] += fwd_rates[206];
  //sp 29
  sp_rates[29] -= fwd_rates[206];

  //rxn 207
  //sp 0
  sp_rates[0] += fwd_rates[207];
  //sp 6
  sp_rates[6] += fwd_rates[207];
  //sp 7
  sp_rates[7] -= fwd_rates[207];
  //sp 11
  sp_rates[11] += fwd_rates[207];
  //sp 29
  sp_rates[29] -= fwd_rates[207];

  //rxn 208
  //sp 10
  sp_rates[10] += (fwd_rates[208] - rev_rates[205]);
  //sp 29
  sp_rates[29] -= (fwd_rates[208] - rev_rates[205]);
  //sp 14
  sp_rates[14] += (fwd_rates[208] - rev_rates[205]);
  //sp 7
  sp_rates[7] -= (fwd_rates[208] - rev_rates[205]);

  //rxn 209
  //sp 11
  sp_rates[11] += (fwd_rates[209] - rev_rates[206]);
  //sp 3
  sp_rates[3] -= (fwd_rates[209] - rev_rates[206]);
  //sp 20
  sp_rates[20] += (fwd_rates[209] - rev_rates[206]);
  //sp 29
  sp_rates[29] -= (fwd_rates[209] - rev_rates[206]);

  //rxn 210
  //sp 1
  sp_rates[1] += (fwd_rates[210] - rev_rates[207]);
  //sp 11
  sp_rates[11] -= (fwd_rates[210] - rev_rates[207]);
  //sp 29
  sp_rates[29] -= (fwd_rates[210] - rev_rates[207]);
  //sp 14
  sp_rates[14] += (fwd_rates[210] - rev_rates[207]);

  //rxn 211
  //sp 8
  sp_rates[8] -= (fwd_rates[211] - rev_rates[208]) * pres_mod[24];
  //sp 13
  sp_rates[13] += 2.0 * (fwd_rates[211] - rev_rates[208]) * pres_mod[24];

  //rxn 212
  //sp 0
  sp_rates[0] += (fwd_rates[212] - rev_rates[209]) * pres_mod[25];
  //sp 5
  sp_rates[5] -= (fwd_rates[212] - rev_rates[209]) * pres_mod[25];
  //sp 13
  sp_rates[13] -= (fwd_rates[212] - rev_rates[209]) * pres_mod[25];

  //sp 30
  (*dy_N) = 0.0;
} // end eval_spec_rates

