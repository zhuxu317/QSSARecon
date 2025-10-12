#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 10
  sp_rates[9] = 2.0 * (fwd_rates[0] - rev_rates[0]) * pres_mod[0];
  //sp 6
  sp_rates[5] = -(fwd_rates[0] - rev_rates[0]) * pres_mod[0];

  //rxn 1
  //sp 10
  sp_rates[9] += (fwd_rates[1] - rev_rates[1]);
  //sp 11
  sp_rates[10] = -(fwd_rates[1] - rev_rates[1]);
  //sp 12
  sp_rates[11] = (fwd_rates[1] - rev_rates[1]);
  //sp 6
  sp_rates[5] -= (fwd_rates[1] - rev_rates[1]);

  //rxn 2
  //sp 10
  sp_rates[9] += (fwd_rates[2] - rev_rates[2]);
  //sp 11
  sp_rates[10] -= (fwd_rates[2] - rev_rates[2]);
  //sp 12
  sp_rates[11] += (fwd_rates[2] - rev_rates[2]);
  //sp 6
  sp_rates[5] -= (fwd_rates[2] - rev_rates[2]);

  //rxn 3
  //sp 0
  sp_rates[0] = (fwd_rates[3] - rev_rates[3]);
  //sp 10
  sp_rates[9] += (fwd_rates[3] - rev_rates[3]);
  //sp 12
  sp_rates[11] -= (fwd_rates[3] - rev_rates[3]);
  //sp 6
  sp_rates[5] -= (fwd_rates[3] - rev_rates[3]);

  //rxn 4
  //sp 10
  sp_rates[9] -= 2.0 * (fwd_rates[4] - rev_rates[4]) * pres_mod[1];
  //sp 6
  sp_rates[5] += (fwd_rates[4] - rev_rates[4]) * pres_mod[1];

  //rxn 5
  //sp 10
  sp_rates[9] -= 2.0 * (fwd_rates[5] - rev_rates[5]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[5] - rev_rates[5]);
  //sp 5
  sp_rates[4] = -(fwd_rates[5] - rev_rates[5]);

  //rxn 6
  //sp 10
  sp_rates[9] -= (fwd_rates[6] - rev_rates[6]) * pres_mod[2];
  //sp 11
  sp_rates[10] -= (fwd_rates[6] - rev_rates[6]) * pres_mod[2];
  //sp 12
  sp_rates[11] += (fwd_rates[6] - rev_rates[6]) * pres_mod[2];

  //rxn 7
  //sp 0
  sp_rates[0] += (fwd_rates[7] - rev_rates[7]) * pres_mod[3];
  //sp 10
  sp_rates[9] -= (fwd_rates[7] - rev_rates[7]) * pres_mod[3];
  //sp 12
  sp_rates[11] -= (fwd_rates[7] - rev_rates[7]) * pres_mod[3];

  //rxn 8
  //sp 11
  sp_rates[10] -= 2.0 * (fwd_rates[8] - rev_rates[8]) * pres_mod[4];
  //sp 5
  sp_rates[4] += (fwd_rates[8] - rev_rates[8]) * pres_mod[4];

  //rxn 9
  //sp 10
  sp_rates[9] -= (fwd_rates[9] - rev_rates[9]);
  //sp 11
  sp_rates[10] += (fwd_rates[9] - rev_rates[9]);
  //sp 12
  sp_rates[11] += (fwd_rates[9] - rev_rates[9]);
  //sp 5
  sp_rates[4] -= (fwd_rates[9] - rev_rates[9]);

  //rxn 10
  //sp 0
  sp_rates[0] += (fwd_rates[10] - rev_rates[10]) * pres_mod[5];
  //sp 10
  sp_rates[9] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[5];
  //sp 12
  sp_rates[11] -= (fwd_rates[10] - rev_rates[10]) * pres_mod[5];

  //rxn 11
  //sp 0
  sp_rates[0] += (fwd_rates[11] - rev_rates[11]);
  //sp 11
  sp_rates[10] += (fwd_rates[11] - rev_rates[11]);
  //sp 12
  sp_rates[11] -= 2.0 * (fwd_rates[11] - rev_rates[11]);

  //rxn 12
  //sp 10
  sp_rates[9] -= (fwd_rates[12] - rev_rates[12]) * pres_mod[6];
  //sp 11
  sp_rates[10] -= (fwd_rates[12] - rev_rates[12]) * pres_mod[6];
  //sp 12
  sp_rates[11] += (fwd_rates[12] - rev_rates[12]) * pres_mod[6];

  //rxn 13
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[13] - rev_rates[13]) * pres_mod[7];
  //sp 14
  sp_rates[13] = -(fwd_rates[13] - rev_rates[13]) * pres_mod[7];

  //rxn 14
  //sp 0
  sp_rates[0] += (fwd_rates[14] - rev_rates[14]);
  //sp 10
  sp_rates[9] -= (fwd_rates[14] - rev_rates[14]);
  //sp 12
  sp_rates[11] += (fwd_rates[14] - rev_rates[14]);
  //sp 14
  sp_rates[13] -= (fwd_rates[14] - rev_rates[14]);

  //rxn 15
  //sp 10
  sp_rates[9] -= (fwd_rates[15] - rev_rates[15]);
  //sp 6
  sp_rates[5] += (fwd_rates[15] - rev_rates[15]);
  //sp 13
  sp_rates[12] = (fwd_rates[15] - rev_rates[15]);
  //sp 14
  sp_rates[13] -= (fwd_rates[15] - rev_rates[15]);

  //rxn 16
  //sp 11
  sp_rates[10] -= (fwd_rates[16] - rev_rates[16]);
  //sp 12
  sp_rates[11] += (fwd_rates[16] - rev_rates[16]);
  //sp 13
  sp_rates[12] += (fwd_rates[16] - rev_rates[16]);
  //sp 14
  sp_rates[13] -= (fwd_rates[16] - rev_rates[16]);

  //rxn 17
  //sp 0
  sp_rates[0] += (fwd_rates[17] - rev_rates[17]);
  //sp 12
  sp_rates[11] -= (fwd_rates[17] - rev_rates[17]);
  //sp 13
  sp_rates[12] += (fwd_rates[17] - rev_rates[17]);
  //sp 14
  sp_rates[13] -= (fwd_rates[17] - rev_rates[17]);

  //rxn 18
  //sp 0
  sp_rates[0] += (fwd_rates[18] - rev_rates[18]);
  //sp 12
  sp_rates[11] -= (fwd_rates[18] - rev_rates[18]);
  //sp 13
  sp_rates[12] += (fwd_rates[18] - rev_rates[18]);
  //sp 14
  sp_rates[13] -= (fwd_rates[18] - rev_rates[18]);

  //rxn 19
  //sp 10
  sp_rates[9] -= (fwd_rates[19] - rev_rates[19]);
  //sp 5
  sp_rates[4] += (fwd_rates[19] - rev_rates[19]);
  //sp 13
  sp_rates[12] -= (fwd_rates[19] - rev_rates[19]);
  //sp 6
  sp_rates[5] += (fwd_rates[19] - rev_rates[19]);

  //rxn 20
  //sp 10
  sp_rates[9] -= (fwd_rates[20] - rev_rates[20]);
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[20] - rev_rates[20]);
  //sp 13
  sp_rates[12] -= (fwd_rates[20] - rev_rates[20]);

  //rxn 21
  //sp 11
  sp_rates[10] -= (fwd_rates[21] - rev_rates[21]);
  //sp 12
  sp_rates[11] += (fwd_rates[21] - rev_rates[21]);
  //sp 13
  sp_rates[12] -= (fwd_rates[21] - rev_rates[21]);
  //sp 5
  sp_rates[4] += (fwd_rates[21] - rev_rates[21]);

  //rxn 22
  //sp 0
  sp_rates[0] += (fwd_rates[22] - rev_rates[22]);
  //sp 12
  sp_rates[11] -= (fwd_rates[22] - rev_rates[22]);
  //sp 13
  sp_rates[12] -= (fwd_rates[22] - rev_rates[22]);
  //sp 5
  sp_rates[4] += (fwd_rates[22] - rev_rates[22]);

  //rxn 23
  //sp 0
  sp_rates[0] += (fwd_rates[23] - rev_rates[23]);
  //sp 12
  sp_rates[11] -= (fwd_rates[23] - rev_rates[23]);
  //sp 13
  sp_rates[12] -= (fwd_rates[23] - rev_rates[23]);
  //sp 5
  sp_rates[4] += (fwd_rates[23] - rev_rates[23]);

  //rxn 24
  //sp 5
  sp_rates[4] += (fwd_rates[24] - rev_rates[24]);
  //sp 13
  sp_rates[12] -= 2.0 * (fwd_rates[24] - rev_rates[24]);
  //sp 14
  sp_rates[13] += (fwd_rates[24] - rev_rates[24]);

  //rxn 25
  //sp 12
  sp_rates[11] += 2.0 * (fwd_rates[25] - rev_rates[25]);
  //sp 13
  sp_rates[12] -= 2.0 * (fwd_rates[25] - rev_rates[25]);
  //sp 5
  sp_rates[4] += (fwd_rates[25] - rev_rates[25]);

  //rxn 26
  //sp 10
  sp_rates[9] -= (fwd_rates[26] - rev_rates[26]) * pres_mod[8];
  //sp 5
  sp_rates[4] -= (fwd_rates[26] - rev_rates[26]) * pres_mod[8];
  //sp 13
  sp_rates[12] += (fwd_rates[26] - rev_rates[26]) * pres_mod[8];

  //rxn 27
  //sp 10
  sp_rates[9] -= (fwd_rates[27] - rev_rates[27]) * pres_mod[9];
  //sp 5
  sp_rates[4] -= (fwd_rates[27] - rev_rates[27]) * pres_mod[9];
  //sp 13
  sp_rates[12] += (fwd_rates[27] - rev_rates[27]) * pres_mod[9];

  //rxn 28
  //sp 10
  sp_rates[9] -= (fwd_rates[28] - rev_rates[28]);
  //sp 5
  sp_rates[4] -= (fwd_rates[28] - rev_rates[28]);
  //sp 13
  sp_rates[12] += (fwd_rates[28] - rev_rates[28]);

  //rxn 29
  //sp 10
  sp_rates[9] -= (fwd_rates[29] - rev_rates[29]) * pres_mod[10];
  //sp 11
  sp_rates[10] -= (fwd_rates[29] - rev_rates[29]) * pres_mod[10];
  //sp 15
  sp_rates[14] = (fwd_rates[29] - rev_rates[29]) * pres_mod[10];

  //rxn 30
  //sp 12
  sp_rates[11] += (fwd_rates[30] - rev_rates[30]);
  //sp 15
  sp_rates[14] -= (fwd_rates[30] - rev_rates[30]);

  //rxn 31
  //sp 12
  sp_rates[11] += (fwd_rates[31] - rev_rates[31]);
  //sp 15
  sp_rates[14] -= (fwd_rates[31] - rev_rates[31]);

  //rxn 32
  //sp 12
  sp_rates[11] += (fwd_rates[32] - rev_rates[32]);
  //sp 15
  sp_rates[14] -= (fwd_rates[32] - rev_rates[32]);

  //rxn 33
  //sp 12
  sp_rates[11] += (fwd_rates[33] - rev_rates[33]);
  //sp 15
  sp_rates[14] -= (fwd_rates[33] - rev_rates[33]);

  //rxn 34
  //sp 12
  sp_rates[11] += (fwd_rates[34] - rev_rates[34]);
  //sp 15
  sp_rates[14] -= (fwd_rates[34] - rev_rates[34]);

  //rxn 35
  //sp 12
  sp_rates[11] += (fwd_rates[35] - rev_rates[35]);
  //sp 15
  sp_rates[14] -= (fwd_rates[35] - rev_rates[35]);

  //rxn 36
  //sp 12
  sp_rates[11] += (fwd_rates[36] - rev_rates[36]);
  //sp 15
  sp_rates[14] -= (fwd_rates[36] - rev_rates[36]);

  //rxn 37
  //sp 12
  sp_rates[11] += (fwd_rates[37] - rev_rates[37]);
  //sp 15
  sp_rates[14] -= (fwd_rates[37] - rev_rates[37]);

  //rxn 38
  //sp 12
  sp_rates[11] += (fwd_rates[38] - rev_rates[38]);
  //sp 15
  sp_rates[14] -= (fwd_rates[38] - rev_rates[38]);

  //rxn 39
  //sp 12
  sp_rates[11] += (fwd_rates[39] - rev_rates[39]);
  //sp 15
  sp_rates[14] -= (fwd_rates[39] - rev_rates[39]);

  //rxn 40
  //sp 10
  sp_rates[9] -= (fwd_rates[40] - rev_rates[40]);
  //sp 2
  sp_rates[2] = (fwd_rates[40] - rev_rates[40]);
  //sp 7
  sp_rates[6] = -(fwd_rates[40] - rev_rates[40]);
  //sp 15
  sp_rates[14] += (fwd_rates[40] - rev_rates[40]);

  //rxn 41
  //sp 10
  sp_rates[9] -= (fwd_rates[41] - rev_rates[41]);
  //sp 4
  (*dy_N) = (fwd_rates[41] - rev_rates[41]);
  //sp 3
  sp_rates[3] = -(fwd_rates[41] - rev_rates[41]);
  //sp 15
  sp_rates[14] += (fwd_rates[41] - rev_rates[41]);

  //rxn 42
  //sp 12
  sp_rates[11] += (fwd_rates[42] - rev_rates[42]);
  //sp 15
  sp_rates[14] -= (fwd_rates[42] - rev_rates[42]);

  //rxn 43
  //sp 12
  sp_rates[11] += (fwd_rates[43] - rev_rates[43]);
  //sp 15
  sp_rates[14] -= (fwd_rates[43] - rev_rates[43]);

  //rxn 44
  //sp 1
  sp_rates[1] = -(fwd_rates[44] - rev_rates[44]) * pres_mod[11];
  //sp 10
  sp_rates[9] += (fwd_rates[44] - rev_rates[44]) * pres_mod[11];
  //sp 16
  sp_rates[15] = (fwd_rates[44] - rev_rates[44]) * pres_mod[11];

  //rxn 45
  //sp 10
  sp_rates[9] -= (fwd_rates[45] - rev_rates[45]);
  //sp 6
  sp_rates[5] += (fwd_rates[45] - rev_rates[45]);
  //sp 1
  sp_rates[1] -= (fwd_rates[45] - rev_rates[45]);
  //sp 16
  sp_rates[15] += (fwd_rates[45] - rev_rates[45]);

  //rxn 46
  //sp 1
  sp_rates[1] -= (fwd_rates[46] - rev_rates[46]);
  //sp 11
  sp_rates[10] -= (fwd_rates[46] - rev_rates[46]);
  //sp 12
  sp_rates[11] += (fwd_rates[46] - rev_rates[46]);
  //sp 16
  sp_rates[15] += (fwd_rates[46] - rev_rates[46]);

  //rxn 47
  //sp 0
  sp_rates[0] += (fwd_rates[47] - rev_rates[47]);
  //sp 1
  sp_rates[1] -= (fwd_rates[47] - rev_rates[47]);
  //sp 12
  sp_rates[11] -= (fwd_rates[47] - rev_rates[47]);
  //sp 16
  sp_rates[15] += (fwd_rates[47] - rev_rates[47]);

  //rxn 48
  //sp 1
  sp_rates[1] -= (fwd_rates[48] - rev_rates[48]);
  //sp 13
  sp_rates[12] -= (fwd_rates[48] - rev_rates[48]);
  //sp 14
  sp_rates[13] += (fwd_rates[48] - rev_rates[48]);
  //sp 16
  sp_rates[15] += (fwd_rates[48] - rev_rates[48]);

  //rxn 49
  //sp 1
  sp_rates[1] -= (fwd_rates[49] - rev_rates[49]);
  //sp 30
  sp_rates[29] = (fwd_rates[49] - rev_rates[49]);
  //sp 7
  sp_rates[6] -= (fwd_rates[49] - rev_rates[49]);
  //sp 16
  sp_rates[15] += (fwd_rates[49] - rev_rates[49]);

  //rxn 50
  //sp 1
  sp_rates[1] -= (fwd_rates[50] - rev_rates[50]);
  //sp 30
  sp_rates[29] += (fwd_rates[50] - rev_rates[50]);
  //sp 7
  sp_rates[6] -= (fwd_rates[50] - rev_rates[50]);
  //sp 16
  sp_rates[15] += (fwd_rates[50] - rev_rates[50]);

  //rxn 51
  //sp 1
  sp_rates[1] -= (fwd_rates[51] - rev_rates[51]);
  //sp 16
  sp_rates[15] += (fwd_rates[51] - rev_rates[51]);
  //sp 7
  sp_rates[6] -= (fwd_rates[51] - rev_rates[51]);
  //sp 31
  sp_rates[30] = (fwd_rates[51] - rev_rates[51]);

  //rxn 52
  //sp 17
  sp_rates[16] = (fwd_rates[52] - rev_rates[52]) * pres_mod[12];
  //sp 10
  sp_rates[9] += (fwd_rates[52] - rev_rates[52]) * pres_mod[12];
  //sp 16
  sp_rates[15] -= (fwd_rates[52] - rev_rates[52]) * pres_mod[12];

  //rxn 53
  //sp 17
  sp_rates[16] += (fwd_rates[53] - rev_rates[53]);
  //sp 10
  sp_rates[9] -= (fwd_rates[53] - rev_rates[53]);
  //sp 6
  sp_rates[5] += (fwd_rates[53] - rev_rates[53]);
  //sp 16
  sp_rates[15] -= (fwd_rates[53] - rev_rates[53]);

  //rxn 54
  //sp 17
  sp_rates[16] += (fwd_rates[54] - rev_rates[54]);
  //sp 11
  sp_rates[10] -= (fwd_rates[54] - rev_rates[54]);
  //sp 12
  sp_rates[11] += (fwd_rates[54] - rev_rates[54]);
  //sp 16
  sp_rates[15] -= (fwd_rates[54] - rev_rates[54]);

  //rxn 55
  //sp 0
  sp_rates[0] += (fwd_rates[55] - rev_rates[55]);
  //sp 17
  sp_rates[16] += (fwd_rates[55] - rev_rates[55]);
  //sp 12
  sp_rates[11] -= (fwd_rates[55] - rev_rates[55]);
  //sp 16
  sp_rates[15] -= (fwd_rates[55] - rev_rates[55]);

  //rxn 56
  //sp 17
  sp_rates[16] += (fwd_rates[56] - rev_rates[56]);
  //sp 1
  sp_rates[1] += (fwd_rates[56] - rev_rates[56]);
  //sp 16
  sp_rates[15] -= 2.0 * (fwd_rates[56] - rev_rates[56]);

  //rxn 57
  //sp 17
  sp_rates[16] -= 2.0 * (fwd_rates[57] - rev_rates[57]);
  //sp 18
  sp_rates[17] = (fwd_rates[57] - rev_rates[57]);
  //sp 16
  sp_rates[15] += (fwd_rates[57] - rev_rates[57]);

  //rxn 58
  //sp 25
  sp_rates[24] = (fwd_rates[58] - rev_rates[58]);
  //sp 10
  sp_rates[9] += (fwd_rates[58] - rev_rates[58]);
  //sp 11
  sp_rates[10] -= (fwd_rates[58] - rev_rates[58]);
  //sp 16
  sp_rates[15] -= (fwd_rates[58] - rev_rates[58]);

  //rxn 59
  //sp 11
  sp_rates[10] -= (fwd_rates[59] - rev_rates[59]);
  //sp 2
  sp_rates[2] += (fwd_rates[59] - rev_rates[59]);
  //sp 6
  sp_rates[5] += (fwd_rates[59] - rev_rates[59]);
  //sp 16
  sp_rates[15] -= (fwd_rates[59] - rev_rates[59]);

  //rxn 60
  //sp 10
  sp_rates[9] += (fwd_rates[60] - rev_rates[60]);
  //sp 27
  sp_rates[26] = (fwd_rates[60] - rev_rates[60]);
  //sp 12
  sp_rates[11] -= (fwd_rates[60] - rev_rates[60]);
  //sp 16
  sp_rates[15] -= (fwd_rates[60] - rev_rates[60]);

  //rxn 61
  //sp 28
  sp_rates[27] = (fwd_rates[61] - rev_rates[61]);
  //sp 10
  sp_rates[9] += (fwd_rates[61] - rev_rates[61]);
  //sp 12
  sp_rates[11] -= (fwd_rates[61] - rev_rates[61]);
  //sp 16
  sp_rates[15] -= (fwd_rates[61] - rev_rates[61]);

  //rxn 62
  //sp 28
  sp_rates[27] += (fwd_rates[62] - rev_rates[62]);
  //sp 10
  sp_rates[9] += (fwd_rates[62] - rev_rates[62]);
  //sp 12
  sp_rates[11] -= (fwd_rates[62] - rev_rates[62]);
  //sp 16
  sp_rates[15] -= (fwd_rates[62] - rev_rates[62]);

  //rxn 63
  //sp 25
  sp_rates[24] += (fwd_rates[63] - rev_rates[63]);
  //sp 12
  sp_rates[11] += (fwd_rates[63] - rev_rates[63]);
  //sp 5
  sp_rates[4] -= (fwd_rates[63] - rev_rates[63]);
  //sp 16
  sp_rates[15] -= (fwd_rates[63] - rev_rates[63]);

  //rxn 64
  //sp 27
  sp_rates[26] += (fwd_rates[64] - rev_rates[64]);
  //sp 11
  sp_rates[10] += (fwd_rates[64] - rev_rates[64]);
  //sp 5
  sp_rates[4] -= (fwd_rates[64] - rev_rates[64]);
  //sp 16
  sp_rates[15] -= (fwd_rates[64] - rev_rates[64]);

  //rxn 65
  //sp 1
  sp_rates[1] += (fwd_rates[65] - rev_rates[65]);
  //sp 13
  sp_rates[12] -= (fwd_rates[65] - rev_rates[65]);
  //sp 5
  sp_rates[4] += (fwd_rates[65] - rev_rates[65]);
  //sp 16
  sp_rates[15] -= (fwd_rates[65] - rev_rates[65]);

  //rxn 66
  //sp 1
  sp_rates[1] += (fwd_rates[66] - rev_rates[66]);
  //sp 13
  sp_rates[12] -= (fwd_rates[66] - rev_rates[66]);
  //sp 5
  sp_rates[4] += (fwd_rates[66] - rev_rates[66]);
  //sp 16
  sp_rates[15] -= (fwd_rates[66] - rev_rates[66]);

  //rxn 67
  //sp 27
  sp_rates[26] += (fwd_rates[67] - rev_rates[67]);
  //sp 12
  sp_rates[11] += (fwd_rates[67] - rev_rates[67]);
  //sp 13
  sp_rates[12] -= (fwd_rates[67] - rev_rates[67]);
  //sp 16
  sp_rates[15] -= (fwd_rates[67] - rev_rates[67]);

  //rxn 68
  //sp 0
  sp_rates[0] += (fwd_rates[68] - rev_rates[68]);
  //sp 25
  sp_rates[24] += (fwd_rates[68] - rev_rates[68]);
  //sp 13
  sp_rates[12] -= (fwd_rates[68] - rev_rates[68]);
  //sp 16
  sp_rates[15] -= (fwd_rates[68] - rev_rates[68]);

  //rxn 69
  //sp 2
  sp_rates[2] -= (fwd_rates[69] - rev_rates[69]);
  //sp 12
  sp_rates[11] += (fwd_rates[69] - rev_rates[69]);
  //sp 23
  sp_rates[22] = (fwd_rates[69] - rev_rates[69]);
  //sp 16
  sp_rates[15] -= (fwd_rates[69] - rev_rates[69]);

  //rxn 70
  //sp 0
  sp_rates[0] += (fwd_rates[70] - rev_rates[70]);
  //sp 2
  sp_rates[2] -= (fwd_rates[70] - rev_rates[70]);
  //sp 4
  (*dy_N) += (fwd_rates[70] - rev_rates[70]);
  //sp 16
  sp_rates[15] -= (fwd_rates[70] - rev_rates[70]);

  //rxn 71
  //sp 27
  sp_rates[26] += (fwd_rates[71] - rev_rates[71]);
  //sp 2
  sp_rates[2] += (fwd_rates[71] - rev_rates[71]);
  //sp 7
  sp_rates[6] -= (fwd_rates[71] - rev_rates[71]);
  //sp 16
  sp_rates[15] -= (fwd_rates[71] - rev_rates[71]);

  //rxn 72
  //sp 0
  sp_rates[0] += (fwd_rates[72] - rev_rates[72]);
  //sp 3
  sp_rates[3] += (fwd_rates[72] - rev_rates[72]);
  //sp 7
  sp_rates[6] -= (fwd_rates[72] - rev_rates[72]);
  //sp 16
  sp_rates[15] -= (fwd_rates[72] - rev_rates[72]);

  //rxn 73
  //sp 17
  sp_rates[16] -= (fwd_rates[73] - rev_rates[73]);
  //sp 10
  sp_rates[9] += (fwd_rates[73] - rev_rates[73]);
  //sp 21
  sp_rates[20] = (fwd_rates[73] - rev_rates[73]);
  //sp 16
  sp_rates[15] -= (fwd_rates[73] - rev_rates[73]);

  //rxn 74
  //sp 18
  sp_rates[17] -= fwd_rates[74];
  //sp 4
  (*dy_N) += fwd_rates[74];
  //sp 10
  sp_rates[9] += 2.0 * fwd_rates[74];
  //sp 16
  sp_rates[15] -= fwd_rates[74];

  //rxn 75
  //sp 18
  sp_rates[17] -= (fwd_rates[75] - rev_rates[74]);
  //sp 4
  (*dy_N) += (fwd_rates[75] - rev_rates[74]);
  //sp 6
  sp_rates[5] += (fwd_rates[75] - rev_rates[74]);
  //sp 16
  sp_rates[15] -= (fwd_rates[75] - rev_rates[74]);

  //rxn 76
  //sp 21
  sp_rates[20] += (fwd_rates[76] - rev_rates[75]);
  //sp 6
  sp_rates[5] += (fwd_rates[76] - rev_rates[75]);
  //sp 16
  sp_rates[15] -= 2.0 * (fwd_rates[76] - rev_rates[75]);

  //rxn 77
  //sp 22
  sp_rates[21] = (fwd_rates[77] - rev_rates[76]);
  //sp 6
  sp_rates[5] += (fwd_rates[77] - rev_rates[76]);
  //sp 16
  sp_rates[15] -= 2.0 * (fwd_rates[77] - rev_rates[76]);

  //rxn 78
  //sp 1
  sp_rates[1] -= (fwd_rates[78] - rev_rates[77]);
  //sp 20
  sp_rates[19] = (fwd_rates[78] - rev_rates[77]);
  //sp 6
  sp_rates[5] += (fwd_rates[78] - rev_rates[77]);
  //sp 16
  sp_rates[15] -= (fwd_rates[78] - rev_rates[77]);

  //rxn 79
  //sp 17
  sp_rates[16] -= (fwd_rates[79] - rev_rates[78]) * pres_mod[13];
  //sp 10
  sp_rates[9] += (fwd_rates[79] - rev_rates[78]) * pres_mod[13];
  //sp 18
  sp_rates[17] += (fwd_rates[79] - rev_rates[78]) * pres_mod[13];

  //rxn 80
  //sp 17
  sp_rates[16] -= (fwd_rates[80] - rev_rates[79]);
  //sp 10
  sp_rates[9] -= (fwd_rates[80] - rev_rates[79]);
  //sp 6
  sp_rates[5] += (fwd_rates[80] - rev_rates[79]);
  //sp 18
  sp_rates[17] += (fwd_rates[80] - rev_rates[79]);

  //rxn 81
  //sp 17
  sp_rates[16] -= (fwd_rates[81] - rev_rates[80]);
  //sp 18
  sp_rates[17] += (fwd_rates[81] - rev_rates[80]);
  //sp 11
  sp_rates[10] -= (fwd_rates[81] - rev_rates[80]);
  //sp 12
  sp_rates[11] += (fwd_rates[81] - rev_rates[80]);

  //rxn 82
  //sp 17
  sp_rates[16] -= (fwd_rates[82] - rev_rates[81]);
  //sp 0
  sp_rates[0] += (fwd_rates[82] - rev_rates[81]);
  //sp 12
  sp_rates[11] -= (fwd_rates[82] - rev_rates[81]);
  //sp 18
  sp_rates[17] += (fwd_rates[82] - rev_rates[81]);

  //rxn 83
  //sp 17
  sp_rates[16] -= (fwd_rates[83] - rev_rates[82]);
  //sp 18
  sp_rates[17] += (fwd_rates[83] - rev_rates[82]);
  //sp 1
  sp_rates[1] += (fwd_rates[83] - rev_rates[82]);
  //sp 16
  sp_rates[15] -= (fwd_rates[83] - rev_rates[82]);

  //rxn 84
  //sp 17
  sp_rates[16] -= (fwd_rates[84] - rev_rates[83]);
  //sp 10
  sp_rates[9] += (fwd_rates[84] - rev_rates[83]);
  //sp 11
  sp_rates[10] -= (fwd_rates[84] - rev_rates[83]);
  //sp 2
  sp_rates[2] += (fwd_rates[84] - rev_rates[83]);

  //rxn 85
  //sp 17
  sp_rates[16] -= (fwd_rates[85] - rev_rates[84]);
  //sp 10
  sp_rates[9] += (fwd_rates[85] - rev_rates[84]);
  //sp 12
  sp_rates[11] -= (fwd_rates[85] - rev_rates[84]);
  //sp 25
  sp_rates[24] += (fwd_rates[85] - rev_rates[84]);

  //rxn 86
  //sp 17
  sp_rates[16] -= (fwd_rates[86] - rev_rates[85]);
  //sp 2
  sp_rates[2] += (fwd_rates[86] - rev_rates[85]);
  //sp 12
  sp_rates[11] -= (fwd_rates[86] - rev_rates[85]);
  //sp 6
  sp_rates[5] += (fwd_rates[86] - rev_rates[85]);

  //rxn 87
  //sp 17
  sp_rates[16] -= (fwd_rates[87] - rev_rates[86]);
  //sp 25
  sp_rates[24] += (fwd_rates[87] - rev_rates[86]);
  //sp 11
  sp_rates[10] += (fwd_rates[87] - rev_rates[86]);
  //sp 5
  sp_rates[4] -= (fwd_rates[87] - rev_rates[86]);

  //rxn 88
  //sp 17
  sp_rates[16] -= (fwd_rates[88] - rev_rates[87]);
  //sp 2
  sp_rates[2] += (fwd_rates[88] - rev_rates[87]);
  //sp 12
  sp_rates[11] += (fwd_rates[88] - rev_rates[87]);
  //sp 5
  sp_rates[4] -= (fwd_rates[88] - rev_rates[87]);

  //rxn 89
  //sp 17
  sp_rates[16] -= 2.0 * fwd_rates[89];
  //sp 4
  (*dy_N) += fwd_rates[89];
  //sp 6
  sp_rates[5] += fwd_rates[89];

  //rxn 90
  //sp 17
  sp_rates[16] -= 2.0 * fwd_rates[90];
  //sp 10
  sp_rates[9] += 2.0 * fwd_rates[90];
  //sp 4
  (*dy_N) += fwd_rates[90];

  //rxn 91
  //sp 17
  sp_rates[16] -= (fwd_rates[91] - rev_rates[88]);
  //sp 18
  sp_rates[17] -= (fwd_rates[91] - rev_rates[88]);
  //sp 4
  (*dy_N) += (fwd_rates[91] - rev_rates[88]);
  //sp 10
  sp_rates[9] += (fwd_rates[91] - rev_rates[88]);

  //rxn 92
  //sp 17
  sp_rates[16] -= (fwd_rates[92] - rev_rates[89]);
  //sp 10
  sp_rates[9] += (fwd_rates[92] - rev_rates[89]);
  //sp 2
  sp_rates[2] -= (fwd_rates[92] - rev_rates[89]);
  //sp 3
  sp_rates[3] += (fwd_rates[92] - rev_rates[89]);

  //rxn 93
  //sp 17
  sp_rates[16] -= (fwd_rates[93] - rev_rates[90]);
  //sp 2
  sp_rates[2] -= (fwd_rates[93] - rev_rates[90]);
  //sp 12
  sp_rates[11] += (fwd_rates[93] - rev_rates[90]);
  //sp 4
  (*dy_N) += (fwd_rates[93] - rev_rates[90]);

  //rxn 94
  //sp 17
  sp_rates[16] -= (fwd_rates[94] - rev_rates[91]);
  //sp 25
  sp_rates[24] += (fwd_rates[94] - rev_rates[91]);
  //sp 2
  sp_rates[2] += (fwd_rates[94] - rev_rates[91]);
  //sp 7
  sp_rates[6] -= (fwd_rates[94] - rev_rates[91]);

  //rxn 95
  //sp 17
  sp_rates[16] -= (fwd_rates[95] - rev_rates[92]);
  //sp 12
  sp_rates[11] += (fwd_rates[95] - rev_rates[92]);
  //sp 3
  sp_rates[3] += (fwd_rates[95] - rev_rates[92]);
  //sp 7
  sp_rates[6] -= (fwd_rates[95] - rev_rates[92]);

  //rxn 96
  //sp 18
  sp_rates[17] += 2.0 * (fwd_rates[96] - rev_rates[93]) * pres_mod[14];
  //sp 4
  (*dy_N) -= (fwd_rates[96] - rev_rates[93]) * pres_mod[14];

  //rxn 97
  //sp 18
  sp_rates[17] -= (fwd_rates[97] - rev_rates[94]);
  //sp 2
  sp_rates[2] += (fwd_rates[97] - rev_rates[94]);
  //sp 12
  sp_rates[11] -= (fwd_rates[97] - rev_rates[94]);
  //sp 10
  sp_rates[9] += (fwd_rates[97] - rev_rates[94]);

  //rxn 98
  //sp 18
  sp_rates[17] -= (fwd_rates[98] - rev_rates[95]);
  //sp 2
  sp_rates[2] += (fwd_rates[98] - rev_rates[95]);
  //sp 11
  sp_rates[10] += (fwd_rates[98] - rev_rates[95]);
  //sp 5
  sp_rates[4] -= (fwd_rates[98] - rev_rates[95]);

  //rxn 99
  //sp 18
  sp_rates[17] -= (fwd_rates[99] - rev_rates[96]);
  //sp 2
  sp_rates[2] -= (fwd_rates[99] - rev_rates[96]);
  //sp 11
  sp_rates[10] += (fwd_rates[99] - rev_rates[96]);
  //sp 4
  (*dy_N) += (fwd_rates[99] - rev_rates[96]);

  //rxn 100
  //sp 18
  sp_rates[17] -= (fwd_rates[100] - rev_rates[97]);
  //sp 2
  sp_rates[2] -= (fwd_rates[100] - rev_rates[97]);
  //sp 11
  sp_rates[10] += (fwd_rates[100] - rev_rates[97]);
  //sp 4
  (*dy_N) += (fwd_rates[100] - rev_rates[97]);

  //rxn 101
  //sp 19
  sp_rates[18] = (fwd_rates[101] - rev_rates[98]) * pres_mod[15];
  //sp 16
  sp_rates[15] -= 2.0 * (fwd_rates[101] - rev_rates[98]) * pres_mod[15];

  //rxn 102
  //sp 19
  sp_rates[18] -= (fwd_rates[102] - rev_rates[99]);
  //sp 22
  sp_rates[21] += (fwd_rates[102] - rev_rates[99]);
  //sp 6
  sp_rates[5] += (fwd_rates[102] - rev_rates[99]);

  //rxn 103
  //sp 10
  sp_rates[9] -= (fwd_rates[103] - rev_rates[100]);
  //sp 19
  sp_rates[18] -= (fwd_rates[103] - rev_rates[100]);
  //sp 20
  sp_rates[19] += (fwd_rates[103] - rev_rates[100]);
  //sp 6
  sp_rates[5] += (fwd_rates[103] - rev_rates[100]);

  //rxn 104
  //sp 20
  sp_rates[19] += (fwd_rates[104] - rev_rates[101]);
  //sp 12
  sp_rates[11] += (fwd_rates[104] - rev_rates[101]);
  //sp 19
  sp_rates[18] -= (fwd_rates[104] - rev_rates[101]);
  //sp 11
  sp_rates[10] -= (fwd_rates[104] - rev_rates[101]);

  //rxn 105
  //sp 0
  sp_rates[0] += (fwd_rates[105] - rev_rates[102]);
  //sp 20
  sp_rates[19] += (fwd_rates[105] - rev_rates[102]);
  //sp 19
  sp_rates[18] -= (fwd_rates[105] - rev_rates[102]);
  //sp 12
  sp_rates[11] -= (fwd_rates[105] - rev_rates[102]);

  //rxn 106
  //sp 19
  sp_rates[18] -= (fwd_rates[106] - rev_rates[103]);
  //sp 20
  sp_rates[19] += (fwd_rates[106] - rev_rates[103]);
  //sp 5
  sp_rates[4] -= (fwd_rates[106] - rev_rates[103]);
  //sp 13
  sp_rates[12] += (fwd_rates[106] - rev_rates[103]);

  //rxn 107
  //sp 1
  sp_rates[1] += (fwd_rates[107] - rev_rates[104]);
  //sp 19
  sp_rates[18] -= (fwd_rates[107] - rev_rates[104]);
  //sp 20
  sp_rates[19] += (fwd_rates[107] - rev_rates[104]);
  //sp 16
  sp_rates[15] -= (fwd_rates[107] - rev_rates[104]);

  //rxn 108
  //sp 17
  sp_rates[16] -= (fwd_rates[108] - rev_rates[105]);
  //sp 19
  sp_rates[18] -= (fwd_rates[108] - rev_rates[105]);
  //sp 20
  sp_rates[19] += (fwd_rates[108] - rev_rates[105]);
  //sp 16
  sp_rates[15] += (fwd_rates[108] - rev_rates[105]);

  //rxn 109
  //sp 17
  sp_rates[16] += (fwd_rates[109] - rev_rates[106]);
  //sp 18
  sp_rates[17] -= (fwd_rates[109] - rev_rates[106]);
  //sp 19
  sp_rates[18] -= (fwd_rates[109] - rev_rates[106]);
  //sp 20
  sp_rates[19] += (fwd_rates[109] - rev_rates[106]);

  //rxn 110
  //sp 25
  sp_rates[24] += (fwd_rates[110] - rev_rates[107]);
  //sp 20
  sp_rates[19] += (fwd_rates[110] - rev_rates[107]);
  //sp 19
  sp_rates[18] -= (fwd_rates[110] - rev_rates[107]);
  //sp 2
  sp_rates[2] -= (fwd_rates[110] - rev_rates[107]);

  //rxn 111
  //sp 19
  sp_rates[18] -= (fwd_rates[111] - rev_rates[108]);
  //sp 20
  sp_rates[19] += (fwd_rates[111] - rev_rates[108]);
  //sp 7
  sp_rates[6] -= (fwd_rates[111] - rev_rates[108]);
  //sp 31
  sp_rates[30] += (fwd_rates[111] - rev_rates[108]);

  //rxn 112
  //sp 19
  sp_rates[18] -= (fwd_rates[112] - rev_rates[109]);
  //sp 20
  sp_rates[19] += (fwd_rates[112] - rev_rates[109]);
  //sp 30
  sp_rates[29] += (fwd_rates[112] - rev_rates[109]);
  //sp 7
  sp_rates[6] -= (fwd_rates[112] - rev_rates[109]);

  //rxn 113
  //sp 0
  sp_rates[0] += (fwd_rates[113] - rev_rates[110]);
  //sp 19
  sp_rates[18] -= (fwd_rates[113] - rev_rates[110]);
  //sp 11
  sp_rates[10] -= (fwd_rates[113] - rev_rates[110]);
  //sp 21
  sp_rates[20] += (fwd_rates[113] - rev_rates[110]);

  //rxn 114
  //sp 10
  sp_rates[9] -= (fwd_rates[114] - rev_rates[111]);
  //sp 19
  sp_rates[18] -= (fwd_rates[114] - rev_rates[111]);
  //sp 1
  sp_rates[1] += (fwd_rates[114] - rev_rates[111]);
  //sp 16
  sp_rates[15] += (fwd_rates[114] - rev_rates[111]);

  //rxn 115
  //sp 17
  sp_rates[16] += (fwd_rates[115] - rev_rates[112]) * pres_mod[16];
  //sp 20
  sp_rates[19] -= (fwd_rates[115] - rev_rates[112]) * pres_mod[16];
  //sp 16
  sp_rates[15] += (fwd_rates[115] - rev_rates[112]) * pres_mod[16];

  //rxn 116
  //sp 10
  sp_rates[9] += (fwd_rates[116] - rev_rates[113]) * pres_mod[17];
  //sp 20
  sp_rates[19] -= (fwd_rates[116] - rev_rates[113]) * pres_mod[17];
  //sp 21
  sp_rates[20] += (fwd_rates[116] - rev_rates[113]) * pres_mod[17];

  //rxn 117
  //sp 10
  sp_rates[9] -= (fwd_rates[117] - rev_rates[114]);
  //sp 20
  sp_rates[19] -= (fwd_rates[117] - rev_rates[114]);
  //sp 21
  sp_rates[20] += (fwd_rates[117] - rev_rates[114]);
  //sp 6
  sp_rates[5] += (fwd_rates[117] - rev_rates[114]);

  //rxn 118
  //sp 10
  sp_rates[9] -= (fwd_rates[118] - rev_rates[115]);
  //sp 22
  sp_rates[21] += (fwd_rates[118] - rev_rates[115]);
  //sp 20
  sp_rates[19] -= (fwd_rates[118] - rev_rates[115]);
  //sp 6
  sp_rates[5] += (fwd_rates[118] - rev_rates[115]);

  //rxn 119
  //sp 12
  sp_rates[11] += (fwd_rates[119] - rev_rates[116]);
  //sp 11
  sp_rates[10] -= (fwd_rates[119] - rev_rates[116]);
  //sp 20
  sp_rates[19] -= (fwd_rates[119] - rev_rates[116]);
  //sp 21
  sp_rates[20] += (fwd_rates[119] - rev_rates[116]);

  //rxn 120
  //sp 12
  sp_rates[11] -= (fwd_rates[120] - rev_rates[117]);
  //sp 0
  sp_rates[0] += (fwd_rates[120] - rev_rates[117]);
  //sp 20
  sp_rates[19] -= (fwd_rates[120] - rev_rates[117]);
  //sp 21
  sp_rates[20] += (fwd_rates[120] - rev_rates[117]);

  //rxn 121
  //sp 12
  sp_rates[11] -= (fwd_rates[121] - rev_rates[118]);
  //sp 0
  sp_rates[0] += (fwd_rates[121] - rev_rates[118]);
  //sp 20
  sp_rates[19] -= (fwd_rates[121] - rev_rates[118]);
  //sp 22
  sp_rates[21] += (fwd_rates[121] - rev_rates[118]);

  //rxn 122
  //sp 21
  sp_rates[20] += (fwd_rates[122] - rev_rates[119]);
  //sp 20
  sp_rates[19] -= (fwd_rates[122] - rev_rates[119]);
  //sp 13
  sp_rates[12] -= (fwd_rates[122] - rev_rates[119]);
  //sp 14
  sp_rates[13] += (fwd_rates[122] - rev_rates[119]);

  //rxn 123
  //sp 17
  sp_rates[16] -= (fwd_rates[123] - rev_rates[120]);
  //sp 20
  sp_rates[19] -= (fwd_rates[123] - rev_rates[120]);
  //sp 21
  sp_rates[20] += (fwd_rates[123] - rev_rates[120]);
  //sp 16
  sp_rates[15] += (fwd_rates[123] - rev_rates[120]);

  //rxn 124
  //sp 1
  sp_rates[1] += (fwd_rates[124] - rev_rates[121]);
  //sp 20
  sp_rates[19] -= (fwd_rates[124] - rev_rates[121]);
  //sp 21
  sp_rates[20] += (fwd_rates[124] - rev_rates[121]);
  //sp 16
  sp_rates[15] -= (fwd_rates[124] - rev_rates[121]);

  //rxn 125
  //sp 1
  sp_rates[1] += (fwd_rates[125] - rev_rates[122]);
  //sp 20
  sp_rates[19] -= (fwd_rates[125] - rev_rates[122]);
  //sp 22
  sp_rates[21] += (fwd_rates[125] - rev_rates[122]);
  //sp 16
  sp_rates[15] -= (fwd_rates[125] - rev_rates[122]);

  //rxn 126
  //sp 17
  sp_rates[16] += (fwd_rates[126] - rev_rates[123]);
  //sp 10
  sp_rates[9] -= (fwd_rates[126] - rev_rates[123]);
  //sp 20
  sp_rates[19] -= (fwd_rates[126] - rev_rates[123]);
  //sp 1
  sp_rates[1] += (fwd_rates[126] - rev_rates[123]);

  //rxn 127
  //sp 25
  sp_rates[24] += (fwd_rates[127] - rev_rates[124]);
  //sp 11
  sp_rates[10] -= (fwd_rates[127] - rev_rates[124]);
  //sp 20
  sp_rates[19] -= (fwd_rates[127] - rev_rates[124]);
  //sp 16
  sp_rates[15] += (fwd_rates[127] - rev_rates[124]);

  //rxn 128
  //sp 2
  sp_rates[2] += fwd_rates[128];
  //sp 10
  sp_rates[9] += fwd_rates[128];
  //sp 11
  sp_rates[10] -= fwd_rates[128];
  //sp 16
  sp_rates[15] += fwd_rates[128];
  //sp 20
  sp_rates[19] -= fwd_rates[128];

  //rxn 129
  //sp 12
  sp_rates[11] -= (fwd_rates[129] - rev_rates[125]);
  //sp 25
  sp_rates[24] += (fwd_rates[129] - rev_rates[125]);
  //sp 20
  sp_rates[19] -= (fwd_rates[129] - rev_rates[125]);
  //sp 1
  sp_rates[1] += (fwd_rates[129] - rev_rates[125]);

  //rxn 130
  //sp 20
  sp_rates[19] += (fwd_rates[130] - rev_rates[126]);
  //sp 12
  sp_rates[11] += (fwd_rates[130] - rev_rates[126]);
  //sp 28
  sp_rates[27] -= (fwd_rates[130] - rev_rates[126]);
  //sp 16
  sp_rates[15] -= (fwd_rates[130] - rev_rates[126]);

  //rxn 131
  //sp 10
  sp_rates[9] += (fwd_rates[131] - rev_rates[127]);
  //sp 22
  sp_rates[21] -= (fwd_rates[131] - rev_rates[127]);
  //sp 23
  sp_rates[22] += (fwd_rates[131] - rev_rates[127]);

  //rxn 132
  //sp 10
  sp_rates[9] += (fwd_rates[132] - rev_rates[128]);
  //sp 22
  sp_rates[21] -= (fwd_rates[132] - rev_rates[128]);
  //sp 23
  sp_rates[22] += (fwd_rates[132] - rev_rates[128]);

  //rxn 133
  //sp 6
  sp_rates[5] += (fwd_rates[133] - rev_rates[129]);
  //sp 4
  (*dy_N) += (fwd_rates[133] - rev_rates[129]);
  //sp 22
  sp_rates[21] -= (fwd_rates[133] - rev_rates[129]);

  //rxn 134
  //sp 21
  sp_rates[20] -= (fwd_rates[134] - rev_rates[130]);
  //sp 22
  sp_rates[21] += (fwd_rates[134] - rev_rates[130]);

  //rxn 135
  //sp 10
  sp_rates[9] -= (fwd_rates[135] - rev_rates[131]);
  //sp 6
  sp_rates[5] += (fwd_rates[135] - rev_rates[131]);
  //sp 22
  sp_rates[21] -= (fwd_rates[135] - rev_rates[131]);
  //sp 23
  sp_rates[22] += (fwd_rates[135] - rev_rates[131]);

  //rxn 136
  //sp 11
  sp_rates[10] -= (fwd_rates[136] - rev_rates[132]);
  //sp 12
  sp_rates[11] += (fwd_rates[136] - rev_rates[132]);
  //sp 22
  sp_rates[21] -= (fwd_rates[136] - rev_rates[132]);
  //sp 23
  sp_rates[22] += (fwd_rates[136] - rev_rates[132]);

  //rxn 137
  //sp 0
  sp_rates[0] += (fwd_rates[137] - rev_rates[133]);
  //sp 12
  sp_rates[11] -= (fwd_rates[137] - rev_rates[133]);
  //sp 22
  sp_rates[21] -= (fwd_rates[137] - rev_rates[133]);
  //sp 23
  sp_rates[22] += (fwd_rates[137] - rev_rates[133]);

  //rxn 138
  //sp 14
  sp_rates[13] += (fwd_rates[138] - rev_rates[134]);
  //sp 13
  sp_rates[12] -= (fwd_rates[138] - rev_rates[134]);
  //sp 22
  sp_rates[21] -= (fwd_rates[138] - rev_rates[134]);
  //sp 23
  sp_rates[22] += (fwd_rates[138] - rev_rates[134]);

  //rxn 139
  //sp 1
  sp_rates[1] += (fwd_rates[139] - rev_rates[135]);
  //sp 22
  sp_rates[21] -= (fwd_rates[139] - rev_rates[135]);
  //sp 23
  sp_rates[22] += (fwd_rates[139] - rev_rates[135]);
  //sp 16
  sp_rates[15] -= (fwd_rates[139] - rev_rates[135]);

  //rxn 140
  //sp 21
  sp_rates[20] += (fwd_rates[140] - rev_rates[136]);
  //sp 22
  sp_rates[21] -= (fwd_rates[140] - rev_rates[136]);

  //rxn 141
  //sp 11
  sp_rates[10] -= (fwd_rates[141] - rev_rates[137]);
  //sp 2
  sp_rates[2] += (fwd_rates[141] - rev_rates[137]);
  //sp 22
  sp_rates[21] -= (fwd_rates[141] - rev_rates[137]);
  //sp 16
  sp_rates[15] += (fwd_rates[141] - rev_rates[137]);

  //rxn 142
  //sp 2
  sp_rates[2] += fwd_rates[142];
  //sp 10
  sp_rates[9] += fwd_rates[142];
  //sp 12
  sp_rates[11] -= fwd_rates[142];
  //sp 16
  sp_rates[15] += fwd_rates[142];
  //sp 22
  sp_rates[21] -= fwd_rates[142];

  //rxn 143
  //sp 5
  sp_rates[4] -= (fwd_rates[143] - rev_rates[138]);
  //sp 22
  sp_rates[21] -= (fwd_rates[143] - rev_rates[138]);
  //sp 7
  sp_rates[6] += (fwd_rates[143] - rev_rates[138]);
  //sp 16
  sp_rates[15] += (fwd_rates[143] - rev_rates[138]);

  //rxn 144
  //sp 2
  sp_rates[2] += fwd_rates[144];
  //sp 12
  sp_rates[11] += fwd_rates[144];
  //sp 13
  sp_rates[12] -= fwd_rates[144];
  //sp 16
  sp_rates[15] += fwd_rates[144];
  //sp 22
  sp_rates[21] -= fwd_rates[144];

  //rxn 145
  //sp 10
  sp_rates[9] += (fwd_rates[145] - rev_rates[139]) * pres_mod[18];
  //sp 21
  sp_rates[20] -= (fwd_rates[145] - rev_rates[139]) * pres_mod[18];
  //sp 23
  sp_rates[22] += (fwd_rates[145] - rev_rates[139]) * pres_mod[18];

  //rxn 146
  //sp 10
  sp_rates[9] -= (fwd_rates[146] - rev_rates[140]);
  //sp 21
  sp_rates[20] -= (fwd_rates[146] - rev_rates[140]);
  //sp 6
  sp_rates[5] += (fwd_rates[146] - rev_rates[140]);
  //sp 23
  sp_rates[22] += (fwd_rates[146] - rev_rates[140]);

  //rxn 147
  //sp 11
  sp_rates[10] -= (fwd_rates[147] - rev_rates[141]);
  //sp 12
  sp_rates[11] += (fwd_rates[147] - rev_rates[141]);
  //sp 21
  sp_rates[20] -= (fwd_rates[147] - rev_rates[141]);
  //sp 23
  sp_rates[22] += (fwd_rates[147] - rev_rates[141]);

  //rxn 148
  //sp 0
  sp_rates[0] += (fwd_rates[148] - rev_rates[142]);
  //sp 12
  sp_rates[11] -= (fwd_rates[148] - rev_rates[142]);
  //sp 21
  sp_rates[20] -= (fwd_rates[148] - rev_rates[142]);
  //sp 23
  sp_rates[22] += (fwd_rates[148] - rev_rates[142]);

  //rxn 149
  //sp 17
  sp_rates[16] -= (fwd_rates[149] - rev_rates[143]);
  //sp 21
  sp_rates[20] -= (fwd_rates[149] - rev_rates[143]);
  //sp 23
  sp_rates[22] += (fwd_rates[149] - rev_rates[143]);
  //sp 16
  sp_rates[15] += (fwd_rates[149] - rev_rates[143]);

  //rxn 150
  //sp 1
  sp_rates[1] += (fwd_rates[150] - rev_rates[144]);
  //sp 21
  sp_rates[20] -= (fwd_rates[150] - rev_rates[144]);
  //sp 23
  sp_rates[22] += (fwd_rates[150] - rev_rates[144]);
  //sp 16
  sp_rates[15] -= (fwd_rates[150] - rev_rates[144]);

  //rxn 151
  //sp 2
  sp_rates[2] -= (fwd_rates[151] - rev_rates[145]);
  //sp 3
  sp_rates[3] += (fwd_rates[151] - rev_rates[145]);
  //sp 21
  sp_rates[20] -= (fwd_rates[151] - rev_rates[145]);
  //sp 16
  sp_rates[15] += (fwd_rates[151] - rev_rates[145]);

  //rxn 152
  //sp 10
  sp_rates[9] += (fwd_rates[152] - rev_rates[146]);
  //sp 4
  (*dy_N) += (fwd_rates[152] - rev_rates[146]);
  //sp 23
  sp_rates[22] -= (fwd_rates[152] - rev_rates[146]);

  //rxn 153
  //sp 11
  sp_rates[10] -= (fwd_rates[153] - rev_rates[147]);
  //sp 12
  sp_rates[11] += (fwd_rates[153] - rev_rates[147]);
  //sp 4
  (*dy_N) += (fwd_rates[153] - rev_rates[147]);
  //sp 23
  sp_rates[22] -= (fwd_rates[153] - rev_rates[147]);

  //rxn 154
  //sp 4
  (*dy_N) += (fwd_rates[154] - rev_rates[148]);
  //sp 5
  sp_rates[4] -= (fwd_rates[154] - rev_rates[148]);
  //sp 13
  sp_rates[12] += (fwd_rates[154] - rev_rates[148]);
  //sp 23
  sp_rates[22] -= (fwd_rates[154] - rev_rates[148]);

  //rxn 155
  //sp 4
  (*dy_N) += (fwd_rates[155] - rev_rates[149]);
  //sp 13
  sp_rates[12] -= (fwd_rates[155] - rev_rates[149]);
  //sp 14
  sp_rates[13] += (fwd_rates[155] - rev_rates[149]);
  //sp 23
  sp_rates[22] -= (fwd_rates[155] - rev_rates[149]);

  //rxn 156
  //sp 17
  sp_rates[16] -= (fwd_rates[156] - rev_rates[150]);
  //sp 4
  (*dy_N) += (fwd_rates[156] - rev_rates[150]);
  //sp 23
  sp_rates[22] -= (fwd_rates[156] - rev_rates[150]);
  //sp 16
  sp_rates[15] += (fwd_rates[156] - rev_rates[150]);

  //rxn 157
  //sp 10
  sp_rates[9] -= (fwd_rates[157] - rev_rates[151]);
  //sp 4
  (*dy_N) += (fwd_rates[157] - rev_rates[151]);
  //sp 6
  sp_rates[5] += (fwd_rates[157] - rev_rates[151]);
  //sp 23
  sp_rates[22] -= (fwd_rates[157] - rev_rates[151]);

  //rxn 158
  //sp 0
  sp_rates[0] += (fwd_rates[158] - rev_rates[152]);
  //sp 4
  (*dy_N) += (fwd_rates[158] - rev_rates[152]);
  //sp 12
  sp_rates[11] -= (fwd_rates[158] - rev_rates[152]);
  //sp 23
  sp_rates[22] -= (fwd_rates[158] - rev_rates[152]);

  //rxn 159
  //sp 25
  sp_rates[24] += (fwd_rates[159] - rev_rates[153]);
  //sp 2
  sp_rates[2] -= (fwd_rates[159] - rev_rates[153]);
  //sp 4
  (*dy_N) += (fwd_rates[159] - rev_rates[153]);
  //sp 23
  sp_rates[22] -= (fwd_rates[159] - rev_rates[153]);

  //rxn 160
  //sp 1
  sp_rates[1] += (fwd_rates[160] - rev_rates[154]);
  //sp 4
  (*dy_N) += (fwd_rates[160] - rev_rates[154]);
  //sp 23
  sp_rates[22] -= (fwd_rates[160] - rev_rates[154]);
  //sp 16
  sp_rates[15] -= (fwd_rates[160] - rev_rates[154]);

  //rxn 161
  //sp 17
  sp_rates[16] += (fwd_rates[161] - rev_rates[155]);
  //sp 11
  sp_rates[10] -= (fwd_rates[161] - rev_rates[155]);
  //sp 2
  sp_rates[2] += (fwd_rates[161] - rev_rates[155]);
  //sp 23
  sp_rates[22] -= (fwd_rates[161] - rev_rates[155]);

  //rxn 162
  //sp 10
  sp_rates[9] += (fwd_rates[162] - rev_rates[156]);
  //sp 11
  sp_rates[10] -= (fwd_rates[162] - rev_rates[156]);
  //sp 3
  sp_rates[3] += (fwd_rates[162] - rev_rates[156]);
  //sp 23
  sp_rates[22] -= (fwd_rates[162] - rev_rates[156]);

  //rxn 163
  //sp 18
  sp_rates[17] -= (fwd_rates[163] - rev_rates[157]) * pres_mod[19];
  //sp 11
  sp_rates[10] -= (fwd_rates[163] - rev_rates[157]) * pres_mod[19];
  //sp 2
  sp_rates[2] += (fwd_rates[163] - rev_rates[157]) * pres_mod[19];

  //rxn 164
  //sp 2
  sp_rates[2] -= (fwd_rates[164] - rev_rates[158]);
  //sp 12
  sp_rates[11] += (fwd_rates[164] - rev_rates[158]);
  //sp 13
  sp_rates[12] -= (fwd_rates[164] - rev_rates[158]);
  //sp 7
  sp_rates[6] += (fwd_rates[164] - rev_rates[158]);

  //rxn 165
  //sp 2
  sp_rates[2] -= (fwd_rates[165] - rev_rates[159]) * pres_mod[20];
  //sp 11
  sp_rates[10] -= (fwd_rates[165] - rev_rates[159]) * pres_mod[20];
  //sp 7
  sp_rates[6] += (fwd_rates[165] - rev_rates[159]) * pres_mod[20];

  //rxn 166
  //sp 10
  sp_rates[9] -= (fwd_rates[166] - rev_rates[160]);
  //sp 2
  sp_rates[2] += (fwd_rates[166] - rev_rates[160]);
  //sp 12
  sp_rates[11] += (fwd_rates[166] - rev_rates[160]);
  //sp 7
  sp_rates[6] -= (fwd_rates[166] - rev_rates[160]);

  //rxn 167
  //sp 11
  sp_rates[10] -= (fwd_rates[167] - rev_rates[161]);
  //sp 2
  sp_rates[2] += (fwd_rates[167] - rev_rates[161]);
  //sp 5
  sp_rates[4] += (fwd_rates[167] - rev_rates[161]);
  //sp 7
  sp_rates[6] -= (fwd_rates[167] - rev_rates[161]);

  //rxn 168
  //sp 11
  sp_rates[10] -= (fwd_rates[168] - rev_rates[162]);
  //sp 2
  sp_rates[2] += (fwd_rates[168] - rev_rates[162]);
  //sp 5
  sp_rates[4] += (fwd_rates[168] - rev_rates[162]);
  //sp 7
  sp_rates[6] -= (fwd_rates[168] - rev_rates[162]);

  //rxn 169
  //sp 10
  sp_rates[9] -= (fwd_rates[169] - rev_rates[163]);
  //sp 12
  sp_rates[11] += (fwd_rates[169] - rev_rates[163]);
  //sp 7
  sp_rates[6] += (fwd_rates[169] - rev_rates[163]);
  //sp 24
  sp_rates[23] = -(fwd_rates[169] - rev_rates[163]);

  //rxn 170
  //sp 11
  sp_rates[10] -= (fwd_rates[170] - rev_rates[164]);
  //sp 12
  sp_rates[11] += (fwd_rates[170] - rev_rates[164]);
  //sp 31
  sp_rates[30] -= (fwd_rates[170] - rev_rates[164]);
  //sp 7
  sp_rates[6] += (fwd_rates[170] - rev_rates[164]);

  //rxn 171
  //sp 2
  sp_rates[2] += 2.0 * (fwd_rates[171] - rev_rates[165]);
  //sp 5
  sp_rates[4] += (fwd_rates[171] - rev_rates[165]);
  //sp 7
  sp_rates[6] -= 2.0 * (fwd_rates[171] - rev_rates[165]);

  //rxn 172
  //sp 2
  sp_rates[2] += (fwd_rates[172] - rev_rates[166]);
  //sp 7
  sp_rates[6] -= 2.0 * (fwd_rates[172] - rev_rates[166]);
  //sp 24
  sp_rates[23] += (fwd_rates[172] - rev_rates[166]);

  //rxn 173
  //sp 12
  sp_rates[11] -= (fwd_rates[173] - rev_rates[167]);
  //sp 13
  sp_rates[12] += (fwd_rates[173] - rev_rates[167]);
  //sp 7
  sp_rates[6] += (fwd_rates[173] - rev_rates[167]);
  //sp 24
  sp_rates[23] -= (fwd_rates[173] - rev_rates[167]);

  //rxn 174
  //sp 10
  sp_rates[9] += (fwd_rates[174] - rev_rates[168]);
  //sp 30
  sp_rates[29] += (fwd_rates[174] - rev_rates[168]);
  //sp 6
  sp_rates[5] -= (fwd_rates[174] - rev_rates[168]);
  //sp 7
  sp_rates[6] -= (fwd_rates[174] - rev_rates[168]);

  //rxn 175
  //sp 10
  sp_rates[9] -= (fwd_rates[175] - rev_rates[169]);
  //sp 6
  sp_rates[5] += (fwd_rates[175] - rev_rates[169]);
  //sp 31
  sp_rates[30] -= (fwd_rates[175] - rev_rates[169]);
  //sp 7
  sp_rates[6] += (fwd_rates[175] - rev_rates[169]);

  //rxn 176
  //sp 5
  sp_rates[4] += (fwd_rates[176] - rev_rates[170]);
  //sp 13
  sp_rates[12] -= (fwd_rates[176] - rev_rates[170]);
  //sp 30
  sp_rates[29] += (fwd_rates[176] - rev_rates[170]);
  //sp 7
  sp_rates[6] -= (fwd_rates[176] - rev_rates[170]);

  //rxn 177
  //sp 13
  sp_rates[12] -= (fwd_rates[177] - rev_rates[171]);
  //sp 5
  sp_rates[4] += (fwd_rates[177] - rev_rates[171]);
  //sp 7
  sp_rates[6] -= (fwd_rates[177] - rev_rates[171]);
  //sp 31
  sp_rates[30] += (fwd_rates[177] - rev_rates[171]);

  //rxn 178
  //sp 11
  sp_rates[10] -= (fwd_rates[178] - rev_rates[172]) * pres_mod[21];
  //sp 7
  sp_rates[6] -= (fwd_rates[178] - rev_rates[172]) * pres_mod[21];
  //sp 24
  sp_rates[23] += (fwd_rates[178] - rev_rates[172]) * pres_mod[21];

  //rxn 179
  //sp 11
  sp_rates[10] -= (fwd_rates[179] - rev_rates[173]);
  //sp 5
  sp_rates[4] += (fwd_rates[179] - rev_rates[173]);
  //sp 7
  sp_rates[6] += (fwd_rates[179] - rev_rates[173]);
  //sp 24
  sp_rates[23] -= (fwd_rates[179] - rev_rates[173]);

  //rxn 180
  //sp 5
  sp_rates[4] += (fwd_rates[180] - rev_rates[174]);
  //sp 7
  sp_rates[6] += (fwd_rates[180] - rev_rates[174]);
  //sp 12
  sp_rates[11] += (fwd_rates[180] - rev_rates[174]);
  //sp 13
  sp_rates[12] -= (fwd_rates[180] - rev_rates[174]);
  //sp 24
  sp_rates[23] -= (fwd_rates[180] - rev_rates[174]);

  //rxn 181
  //sp 2
  sp_rates[2] += (fwd_rates[181] - rev_rates[175]) * pres_mod[22];
  //sp 5
  sp_rates[4] += (fwd_rates[181] - rev_rates[175]) * pres_mod[22];
  //sp 24
  sp_rates[23] -= (fwd_rates[181] - rev_rates[175]) * pres_mod[22];

  //rxn 182
  //sp 11
  sp_rates[10] += (fwd_rates[182] - rev_rates[176]) * pres_mod[23];
  //sp 4
  (*dy_N) += (fwd_rates[182] - rev_rates[176]) * pres_mod[23];
  //sp 3
  sp_rates[3] -= (fwd_rates[182] - rev_rates[176]) * pres_mod[23];

  //rxn 183
  //sp 12
  sp_rates[11] += (fwd_rates[183] - rev_rates[177]);
  //sp 10
  sp_rates[9] -= (fwd_rates[183] - rev_rates[177]);
  //sp 4
  (*dy_N) += (fwd_rates[183] - rev_rates[177]);
  //sp 3
  sp_rates[3] -= (fwd_rates[183] - rev_rates[177]);

  //rxn 184
  //sp 12
  sp_rates[11] += (fwd_rates[184] - rev_rates[178]);
  //sp 10
  sp_rates[9] -= (fwd_rates[184] - rev_rates[178]);
  //sp 4
  (*dy_N) += (fwd_rates[184] - rev_rates[178]);
  //sp 3
  sp_rates[3] -= (fwd_rates[184] - rev_rates[178]);

  //rxn 185
  //sp 2
  sp_rates[2] += 2.0 * (fwd_rates[185] - rev_rates[179]);
  //sp 11
  sp_rates[10] -= (fwd_rates[185] - rev_rates[179]);
  //sp 3
  sp_rates[3] -= (fwd_rates[185] - rev_rates[179]);

  //rxn 186
  //sp 11
  sp_rates[10] -= (fwd_rates[186] - rev_rates[180]);
  //sp 3
  sp_rates[3] -= (fwd_rates[186] - rev_rates[180]);
  //sp 4
  (*dy_N) += (fwd_rates[186] - rev_rates[180]);
  //sp 5
  sp_rates[4] += (fwd_rates[186] - rev_rates[180]);

  //rxn 187
  //sp 12
  sp_rates[11] -= (fwd_rates[187] - rev_rates[181]);
  //sp 4
  (*dy_N) += (fwd_rates[187] - rev_rates[181]);
  //sp 3
  sp_rates[3] -= (fwd_rates[187] - rev_rates[181]);
  //sp 13
  sp_rates[12] += (fwd_rates[187] - rev_rates[181]);

  //rxn 188
  //sp 12
  sp_rates[11] -= (fwd_rates[188] - rev_rates[182]);
  //sp 25
  sp_rates[24] += (fwd_rates[188] - rev_rates[182]);
  //sp 2
  sp_rates[2] += (fwd_rates[188] - rev_rates[182]);
  //sp 3
  sp_rates[3] -= (fwd_rates[188] - rev_rates[182]);

  //rxn 189
  //sp 2
  sp_rates[2] -= (fwd_rates[189] - rev_rates[183]);
  //sp 3
  sp_rates[3] -= (fwd_rates[189] - rev_rates[183]);
  //sp 4
  (*dy_N) += (fwd_rates[189] - rev_rates[183]);
  //sp 7
  sp_rates[6] += (fwd_rates[189] - rev_rates[183]);

  //rxn 190
  //sp 25
  sp_rates[24] -= (fwd_rates[190] - rev_rates[184]);
  //sp 10
  sp_rates[9] += (fwd_rates[190] - rev_rates[184]);
  //sp 2
  sp_rates[2] += (fwd_rates[190] - rev_rates[184]);

  //rxn 191
  //sp 25
  sp_rates[24] -= (fwd_rates[191] - rev_rates[185]);
  //sp 10
  sp_rates[9] -= (fwd_rates[191] - rev_rates[185]);
  //sp 2
  sp_rates[2] += (fwd_rates[191] - rev_rates[185]);
  //sp 6
  sp_rates[5] += (fwd_rates[191] - rev_rates[185]);

  //rxn 192
  //sp 25
  sp_rates[24] += (fwd_rates[192] - rev_rates[186]);
  //sp 1
  sp_rates[1] -= (fwd_rates[192] - rev_rates[186]);
  //sp 2
  sp_rates[2] -= (fwd_rates[192] - rev_rates[186]);
  //sp 16
  sp_rates[15] += (fwd_rates[192] - rev_rates[186]);

  //rxn 193
  //sp 25
  sp_rates[24] -= (fwd_rates[193] - rev_rates[187]);
  //sp 12
  sp_rates[11] += (fwd_rates[193] - rev_rates[187]);
  //sp 11
  sp_rates[10] -= (fwd_rates[193] - rev_rates[187]);
  //sp 2
  sp_rates[2] += (fwd_rates[193] - rev_rates[187]);

  //rxn 194
  //sp 25
  sp_rates[24] -= (fwd_rates[194] - rev_rates[188]);
  //sp 0
  sp_rates[0] += (fwd_rates[194] - rev_rates[188]);
  //sp 2
  sp_rates[2] += (fwd_rates[194] - rev_rates[188]);
  //sp 12
  sp_rates[11] -= (fwd_rates[194] - rev_rates[188]);

  //rxn 195
  //sp 25
  sp_rates[24] -= (fwd_rates[195] - rev_rates[189]);
  //sp 0
  sp_rates[0] += (fwd_rates[195] - rev_rates[189]);
  //sp 2
  sp_rates[2] += (fwd_rates[195] - rev_rates[189]);
  //sp 12
  sp_rates[11] -= (fwd_rates[195] - rev_rates[189]);

  //rxn 196
  //sp 25
  sp_rates[24] -= (fwd_rates[196] - rev_rates[190]);
  //sp 2
  sp_rates[2] += (fwd_rates[196] - rev_rates[190]);
  //sp 5
  sp_rates[4] -= (fwd_rates[196] - rev_rates[190]);
  //sp 13
  sp_rates[12] += (fwd_rates[196] - rev_rates[190]);

  //rxn 197
  //sp 25
  sp_rates[24] -= (fwd_rates[197] - rev_rates[191]);
  //sp 2
  sp_rates[2] += (fwd_rates[197] - rev_rates[191]);
  //sp 30
  sp_rates[29] += (fwd_rates[197] - rev_rates[191]);
  //sp 7
  sp_rates[6] -= (fwd_rates[197] - rev_rates[191]);

  //rxn 198
  //sp 25
  sp_rates[24] -= (fwd_rates[198] - rev_rates[192]);
  //sp 10
  sp_rates[9] += (fwd_rates[198] - rev_rates[192]);
  //sp 12
  sp_rates[11] -= (fwd_rates[198] - rev_rates[192]);
  //sp 30
  sp_rates[29] += (fwd_rates[198] - rev_rates[192]);

  //rxn 199
  //sp 25
  sp_rates[24] -= 2.0 * (fwd_rates[199] - rev_rates[193]);
  //sp 0
  sp_rates[0] += (fwd_rates[199] - rev_rates[193]);
  //sp 3
  sp_rates[3] += (fwd_rates[199] - rev_rates[193]);

  //rxn 200
  //sp 26
  sp_rates[25] = -(fwd_rates[200] - rev_rates[194]) * pres_mod[24];
  //sp 2
  sp_rates[2] += (fwd_rates[200] - rev_rates[194]) * pres_mod[24];
  //sp 10
  sp_rates[9] += (fwd_rates[200] - rev_rates[194]) * pres_mod[24];

  //rxn 201
  //sp 12
  sp_rates[11] += (fwd_rates[201] - rev_rates[195]);
  //sp 26
  sp_rates[25] -= (fwd_rates[201] - rev_rates[195]);
  //sp 11
  sp_rates[10] -= (fwd_rates[201] - rev_rates[195]);
  //sp 2
  sp_rates[2] += (fwd_rates[201] - rev_rates[195]);

  //rxn 202
  //sp 25
  sp_rates[24] += (fwd_rates[202] - rev_rates[196]);
  //sp 26
  sp_rates[25] -= (fwd_rates[202] - rev_rates[196]);

  //rxn 203
  //sp 26
  sp_rates[25] -= (fwd_rates[203] - rev_rates[197]);
  //sp 12
  sp_rates[11] -= (fwd_rates[203] - rev_rates[197]);
  //sp 30
  sp_rates[29] += (fwd_rates[203] - rev_rates[197]);
  //sp 10
  sp_rates[9] += (fwd_rates[203] - rev_rates[197]);

  //rxn 204
  //sp 26
  sp_rates[25] -= (fwd_rates[204] - rev_rates[198]);
  //sp 12
  sp_rates[11] += (fwd_rates[204] - rev_rates[198]);
  //sp 5
  sp_rates[4] -= (fwd_rates[204] - rev_rates[198]);
  //sp 7
  sp_rates[6] += (fwd_rates[204] - rev_rates[198]);

  //rxn 205
  //sp 25
  sp_rates[24] += (fwd_rates[205] - rev_rates[199]) * pres_mod[25];
  //sp 10
  sp_rates[9] += (fwd_rates[205] - rev_rates[199]) * pres_mod[25];
  //sp 27
  sp_rates[26] -= (fwd_rates[205] - rev_rates[199]) * pres_mod[25];

  //rxn 206
  //sp 27
  sp_rates[26] -= (fwd_rates[206] - rev_rates[200]) * pres_mod[26];
  //sp 28
  sp_rates[27] += (fwd_rates[206] - rev_rates[200]) * pres_mod[26];

  //rxn 207
  //sp 25
  sp_rates[24] += (fwd_rates[207] - rev_rates[201]);
  //sp 10
  sp_rates[9] -= (fwd_rates[207] - rev_rates[201]);
  //sp 27
  sp_rates[26] -= (fwd_rates[207] - rev_rates[201]);
  //sp 6
  sp_rates[5] += (fwd_rates[207] - rev_rates[201]);

  //rxn 208
  //sp 25
  sp_rates[24] += (fwd_rates[208] - rev_rates[202]);
  //sp 12
  sp_rates[11] += (fwd_rates[208] - rev_rates[202]);
  //sp 27
  sp_rates[26] -= (fwd_rates[208] - rev_rates[202]);
  //sp 11
  sp_rates[10] -= (fwd_rates[208] - rev_rates[202]);

  //rxn 209
  //sp 0
  sp_rates[0] += (fwd_rates[209] - rev_rates[203]);
  //sp 25
  sp_rates[24] += (fwd_rates[209] - rev_rates[203]);
  //sp 27
  sp_rates[26] -= (fwd_rates[209] - rev_rates[203]);
  //sp 12
  sp_rates[11] -= (fwd_rates[209] - rev_rates[203]);

  //rxn 210
  //sp 25
  sp_rates[24] += (fwd_rates[210] - rev_rates[204]);
  //sp 27
  sp_rates[26] -= (fwd_rates[210] - rev_rates[204]);
  //sp 5
  sp_rates[4] -= (fwd_rates[210] - rev_rates[204]);
  //sp 13
  sp_rates[12] += (fwd_rates[210] - rev_rates[204]);

  //rxn 211
  //sp 25
  sp_rates[24] += (fwd_rates[211] - rev_rates[205]);
  //sp 27
  sp_rates[26] -= (fwd_rates[211] - rev_rates[205]);
  //sp 13
  sp_rates[12] -= (fwd_rates[211] - rev_rates[205]);
  //sp 14
  sp_rates[13] += (fwd_rates[211] - rev_rates[205]);

  //rxn 212
  //sp 25
  sp_rates[24] += (fwd_rates[212] - rev_rates[206]);
  //sp 27
  sp_rates[26] -= (fwd_rates[212] - rev_rates[206]);
  //sp 30
  sp_rates[29] += (fwd_rates[212] - rev_rates[206]);
  //sp 7
  sp_rates[6] -= (fwd_rates[212] - rev_rates[206]);

  //rxn 213
  //sp 25
  sp_rates[24] += (fwd_rates[213] - rev_rates[207]);
  //sp 1
  sp_rates[1] += (fwd_rates[213] - rev_rates[207]);
  //sp 27
  sp_rates[26] -= (fwd_rates[213] - rev_rates[207]);
  //sp 16
  sp_rates[15] -= (fwd_rates[213] - rev_rates[207]);

  //rxn 214
  //sp 25
  sp_rates[24] += 2.0 * (fwd_rates[214] - rev_rates[208]);
  //sp 27
  sp_rates[26] -= (fwd_rates[214] - rev_rates[208]);
  //sp 2
  sp_rates[2] -= (fwd_rates[214] - rev_rates[208]);

  //rxn 215
  //sp 25
  sp_rates[24] += (fwd_rates[215] - rev_rates[209]) * pres_mod[27];
  //sp 10
  sp_rates[9] += (fwd_rates[215] - rev_rates[209]) * pres_mod[27];
  //sp 28
  sp_rates[27] -= (fwd_rates[215] - rev_rates[209]) * pres_mod[27];

  //rxn 216
  //sp 25
  sp_rates[24] += (fwd_rates[216] - rev_rates[210]);
  //sp 10
  sp_rates[9] -= (fwd_rates[216] - rev_rates[210]);
  //sp 28
  sp_rates[27] -= (fwd_rates[216] - rev_rates[210]);
  //sp 6
  sp_rates[5] += (fwd_rates[216] - rev_rates[210]);

  //rxn 217
  //sp 25
  sp_rates[24] += (fwd_rates[217] - rev_rates[211]);
  //sp 12
  sp_rates[11] += (fwd_rates[217] - rev_rates[211]);
  //sp 11
  sp_rates[10] -= (fwd_rates[217] - rev_rates[211]);
  //sp 28
  sp_rates[27] -= (fwd_rates[217] - rev_rates[211]);

  //rxn 218
  //sp 25
  sp_rates[24] += (fwd_rates[218] - rev_rates[212]);
  //sp 12
  sp_rates[11] += (fwd_rates[218] - rev_rates[212]);
  //sp 11
  sp_rates[10] -= (fwd_rates[218] - rev_rates[212]);
  //sp 28
  sp_rates[27] -= (fwd_rates[218] - rev_rates[212]);

  //rxn 219
  //sp 12
  sp_rates[11] -= (fwd_rates[219] - rev_rates[213]);
  //sp 0
  sp_rates[0] += (fwd_rates[219] - rev_rates[213]);
  //sp 28
  sp_rates[27] -= (fwd_rates[219] - rev_rates[213]);
  //sp 25
  sp_rates[24] += (fwd_rates[219] - rev_rates[213]);

  //rxn 220
  //sp 25
  sp_rates[24] += (fwd_rates[220] - rev_rates[214]);
  //sp 28
  sp_rates[27] -= (fwd_rates[220] - rev_rates[214]);
  //sp 5
  sp_rates[4] -= (fwd_rates[220] - rev_rates[214]);
  //sp 13
  sp_rates[12] += (fwd_rates[220] - rev_rates[214]);

  //rxn 221
  //sp 25
  sp_rates[24] += (fwd_rates[221] - rev_rates[215]);
  //sp 28
  sp_rates[27] -= (fwd_rates[221] - rev_rates[215]);
  //sp 13
  sp_rates[12] -= (fwd_rates[221] - rev_rates[215]);
  //sp 14
  sp_rates[13] += (fwd_rates[221] - rev_rates[215]);

  //rxn 222
  //sp 25
  sp_rates[24] += (fwd_rates[222] - rev_rates[216]);
  //sp 1
  sp_rates[1] += (fwd_rates[222] - rev_rates[216]);
  //sp 28
  sp_rates[27] -= (fwd_rates[222] - rev_rates[216]);
  //sp 16
  sp_rates[15] -= (fwd_rates[222] - rev_rates[216]);

  //rxn 223
  //sp 25
  sp_rates[24] += (fwd_rates[223] - rev_rates[217]);
  //sp 28
  sp_rates[27] -= (fwd_rates[223] - rev_rates[217]);
  //sp 30
  sp_rates[29] += (fwd_rates[223] - rev_rates[217]);
  //sp 7
  sp_rates[6] -= (fwd_rates[223] - rev_rates[217]);

  //rxn 224
  //sp 0
  sp_rates[0] += (fwd_rates[224] - rev_rates[218]);
  //sp 28
  sp_rates[27] -= (fwd_rates[224] - rev_rates[218]);
  //sp 22
  sp_rates[21] += (fwd_rates[224] - rev_rates[218]);
  //sp 16
  sp_rates[15] -= (fwd_rates[224] - rev_rates[218]);

  //rxn 225
  //sp 2
  sp_rates[2] -= (fwd_rates[225] - rev_rates[219]);
  //sp 12
  sp_rates[11] -= (fwd_rates[225] - rev_rates[219]);
  //sp 30
  sp_rates[29] += (fwd_rates[225] - rev_rates[219]);

  //rxn 226
  //sp 11
  sp_rates[10] -= (fwd_rates[226] - rev_rates[220]);
  //sp 12
  sp_rates[11] += (fwd_rates[226] - rev_rates[220]);
  //sp 30
  sp_rates[29] -= (fwd_rates[226] - rev_rates[220]);
  //sp 7
  sp_rates[6] += (fwd_rates[226] - rev_rates[220]);

  //rxn 227
  //sp 0
  sp_rates[0] += (fwd_rates[227] - rev_rates[221]);
  //sp 12
  sp_rates[11] -= (fwd_rates[227] - rev_rates[221]);
  //sp 30
  sp_rates[29] -= (fwd_rates[227] - rev_rates[221]);
  //sp 7
  sp_rates[6] += (fwd_rates[227] - rev_rates[221]);

  //rxn 228
  //sp 2
  sp_rates[2] += (fwd_rates[228] - rev_rates[222]);
  //sp 30
  sp_rates[29] -= (fwd_rates[228] - rev_rates[222]);
  //sp 7
  sp_rates[6] -= (fwd_rates[228] - rev_rates[222]);
  //sp 32
  sp_rates[31] = (fwd_rates[228] - rev_rates[222]);

  //rxn 229
  //sp 0
  sp_rates[0] += (fwd_rates[229] - rev_rates[223]);
  //sp 2
  sp_rates[2] += (fwd_rates[229] - rev_rates[223]);
  //sp 30
  sp_rates[29] -= 2.0 * (fwd_rates[229] - rev_rates[223]);
  //sp 7
  sp_rates[6] += (fwd_rates[229] - rev_rates[223]);

  //rxn 230
  //sp 0
  sp_rates[0] += (fwd_rates[230] - rev_rates[224]);
  //sp 10
  sp_rates[9] -= (fwd_rates[230] - rev_rates[224]);
  //sp 2
  sp_rates[2] += (fwd_rates[230] - rev_rates[224]);
  //sp 30
  sp_rates[29] -= (fwd_rates[230] - rev_rates[224]);

  //rxn 231
  //sp 2
  sp_rates[2] -= (fwd_rates[231] - rev_rates[225]);
  //sp 12
  sp_rates[11] -= (fwd_rates[231] - rev_rates[225]);
  //sp 31
  sp_rates[30] += (fwd_rates[231] - rev_rates[225]);

  //rxn 232
  //sp 30
  sp_rates[29] += (fwd_rates[232] - rev_rates[226]);
  //sp 31
  sp_rates[30] -= (fwd_rates[232] - rev_rates[226]);

  //rxn 233
  //sp 0
  sp_rates[0] += (fwd_rates[233] - rev_rates[227]);
  //sp 12
  sp_rates[11] -= (fwd_rates[233] - rev_rates[227]);
  //sp 31
  sp_rates[30] -= (fwd_rates[233] - rev_rates[227]);
  //sp 7
  sp_rates[6] += (fwd_rates[233] - rev_rates[227]);

  //rxn 234
  //sp 0
  sp_rates[0] += (fwd_rates[234] - rev_rates[228]);
  //sp 10
  sp_rates[9] -= (fwd_rates[234] - rev_rates[228]);
  //sp 2
  sp_rates[2] += (fwd_rates[234] - rev_rates[228]);
  //sp 31
  sp_rates[30] -= (fwd_rates[234] - rev_rates[228]);

  //rxn 235
  //sp 25
  sp_rates[24] += (fwd_rates[235] - rev_rates[229]);
  //sp 10
  sp_rates[9] -= (fwd_rates[235] - rev_rates[229]);
  //sp 12
  sp_rates[11] += (fwd_rates[235] - rev_rates[229]);
  //sp 31
  sp_rates[30] -= (fwd_rates[235] - rev_rates[229]);

  //rxn 236
  //sp 12
  sp_rates[11] -= (fwd_rates[236] - rev_rates[230]) * pres_mod[28];
  //sp 7
  sp_rates[6] -= (fwd_rates[236] - rev_rates[230]) * pres_mod[28];
  //sp 32
  sp_rates[31] += (fwd_rates[236] - rev_rates[230]) * pres_mod[28];

  //rxn 237
  //sp 10
  sp_rates[9] -= (fwd_rates[237] - rev_rates[231]);
  //sp 24
  sp_rates[23] += (fwd_rates[237] - rev_rates[231]);
  //sp 6
  sp_rates[5] += (fwd_rates[237] - rev_rates[231]);
  //sp 32
  sp_rates[31] -= (fwd_rates[237] - rev_rates[231]);

  //rxn 238
  //sp 0
  sp_rates[0] += (fwd_rates[238] - rev_rates[232]);
  //sp 10
  sp_rates[9] -= (fwd_rates[238] - rev_rates[232]);
  //sp 7
  sp_rates[6] += (fwd_rates[238] - rev_rates[232]);
  //sp 32
  sp_rates[31] -= (fwd_rates[238] - rev_rates[232]);

  //rxn 239
  //sp 10
  sp_rates[9] -= (fwd_rates[239] - rev_rates[233]);
  //sp 12
  sp_rates[11] += (fwd_rates[239] - rev_rates[233]);
  //sp 30
  sp_rates[29] += (fwd_rates[239] - rev_rates[233]);
  //sp 32
  sp_rates[31] -= (fwd_rates[239] - rev_rates[233]);

  //rxn 240
  //sp 0
  sp_rates[0] += (fwd_rates[240] - rev_rates[234]);
  //sp 12
  sp_rates[11] -= (fwd_rates[240] - rev_rates[234]);
  //sp 24
  sp_rates[23] += (fwd_rates[240] - rev_rates[234]);
  //sp 32
  sp_rates[31] -= (fwd_rates[240] - rev_rates[234]);

  //rxn 241
  //sp 12
  sp_rates[11] += (fwd_rates[241] - rev_rates[235]) * pres_mod[29];
  //sp 29
  sp_rates[28] = -(fwd_rates[241] - rev_rates[235]) * pres_mod[29];
  //sp 16
  sp_rates[15] += (fwd_rates[241] - rev_rates[235]) * pres_mod[29];

  //rxn 242
  //sp 10
  sp_rates[9] -= (fwd_rates[242] - rev_rates[236]);
  //sp 28
  sp_rates[27] += (fwd_rates[242] - rev_rates[236]);
  //sp 29
  sp_rates[28] -= (fwd_rates[242] - rev_rates[236]);
  //sp 6
  sp_rates[5] += (fwd_rates[242] - rev_rates[236]);

  //rxn 243
  //sp 10
  sp_rates[9] -= (fwd_rates[243] - rev_rates[237]);
  //sp 27
  sp_rates[26] += (fwd_rates[243] - rev_rates[237]);
  //sp 29
  sp_rates[28] -= (fwd_rates[243] - rev_rates[237]);
  //sp 6
  sp_rates[5] += (fwd_rates[243] - rev_rates[237]);

  //rxn 244
  //sp 12
  sp_rates[11] += (fwd_rates[244] - rev_rates[238]);
  //sp 11
  sp_rates[10] -= (fwd_rates[244] - rev_rates[238]);
  //sp 28
  sp_rates[27] += (fwd_rates[244] - rev_rates[238]);
  //sp 29
  sp_rates[28] -= (fwd_rates[244] - rev_rates[238]);

  //rxn 245
  //sp 12
  sp_rates[11] += (fwd_rates[245] - rev_rates[239]);
  //sp 11
  sp_rates[10] -= (fwd_rates[245] - rev_rates[239]);
  //sp 27
  sp_rates[26] += (fwd_rates[245] - rev_rates[239]);
  //sp 29
  sp_rates[28] -= (fwd_rates[245] - rev_rates[239]);

  //rxn 246
  //sp 0
  sp_rates[0] += (fwd_rates[246] - rev_rates[240]);
  //sp 28
  sp_rates[27] += (fwd_rates[246] - rev_rates[240]);
  //sp 12
  sp_rates[11] -= (fwd_rates[246] - rev_rates[240]);
  //sp 29
  sp_rates[28] -= (fwd_rates[246] - rev_rates[240]);

  //rxn 247
  //sp 0
  sp_rates[0] += (fwd_rates[247] - rev_rates[241]);
  //sp 27
  sp_rates[26] += (fwd_rates[247] - rev_rates[241]);
  //sp 12
  sp_rates[11] -= (fwd_rates[247] - rev_rates[241]);
  //sp 29
  sp_rates[28] -= (fwd_rates[247] - rev_rates[241]);

  //rxn 248
  //sp 14
  sp_rates[13] += (fwd_rates[248] - rev_rates[242]);
  //sp 28
  sp_rates[27] += (fwd_rates[248] - rev_rates[242]);
  //sp 13
  sp_rates[12] -= (fwd_rates[248] - rev_rates[242]);
  //sp 29
  sp_rates[28] -= (fwd_rates[248] - rev_rates[242]);

  //rxn 249
  //sp 27
  sp_rates[26] += (fwd_rates[249] - rev_rates[243]);
  //sp 14
  sp_rates[13] += (fwd_rates[249] - rev_rates[243]);
  //sp 13
  sp_rates[12] -= (fwd_rates[249] - rev_rates[243]);
  //sp 29
  sp_rates[28] -= (fwd_rates[249] - rev_rates[243]);

  //rxn 250
  //sp 1
  sp_rates[1] += (fwd_rates[250] - rev_rates[244]);
  //sp 28
  sp_rates[27] += (fwd_rates[250] - rev_rates[244]);
  //sp 29
  sp_rates[28] -= (fwd_rates[250] - rev_rates[244]);
  //sp 16
  sp_rates[15] -= (fwd_rates[250] - rev_rates[244]);

  //rxn 251
  //sp 1
  sp_rates[1] += (fwd_rates[251] - rev_rates[245]);
  //sp 27
  sp_rates[26] += (fwd_rates[251] - rev_rates[245]);
  //sp 29
  sp_rates[28] -= (fwd_rates[251] - rev_rates[245]);
  //sp 16
  sp_rates[15] -= (fwd_rates[251] - rev_rates[245]);

  //rxn 252
  //sp 17
  sp_rates[16] -= (fwd_rates[252] - rev_rates[246]);
  //sp 28
  sp_rates[27] += (fwd_rates[252] - rev_rates[246]);
  //sp 29
  sp_rates[28] -= (fwd_rates[252] - rev_rates[246]);
  //sp 16
  sp_rates[15] += (fwd_rates[252] - rev_rates[246]);

  //rxn 253
  //sp 17
  sp_rates[16] -= (fwd_rates[253] - rev_rates[247]);
  //sp 27
  sp_rates[26] += (fwd_rates[253] - rev_rates[247]);
  //sp 29
  sp_rates[28] -= (fwd_rates[253] - rev_rates[247]);
  //sp 16
  sp_rates[15] += (fwd_rates[253] - rev_rates[247]);

  //rxn 254
  //sp 36
  sp_rates[35] = -(fwd_rates[254] - rev_rates[248]);
  //sp 37
  sp_rates[36] = (fwd_rates[254] - rev_rates[248]);

  //rxn 255
  //sp 36
  sp_rates[35] -= (fwd_rates[255] - rev_rates[249]);
  //sp 38
  sp_rates[37] = (fwd_rates[255] - rev_rates[249]);

  //rxn 256
  //sp 36
  sp_rates[35] -= (fwd_rates[256] - rev_rates[250]);
  //sp 7
  sp_rates[6] += (fwd_rates[256] - rev_rates[250]);
  //sp 16
  sp_rates[15] += (fwd_rates[256] - rev_rates[250]);

  //rxn 257
  //sp 0
  sp_rates[0] += (fwd_rates[257] - rev_rates[251]);
  //sp 3
  sp_rates[3] += (fwd_rates[257] - rev_rates[251]);
  //sp 36
  sp_rates[35] -= (fwd_rates[257] - rev_rates[251]);

  //rxn 258
  //sp 37
  sp_rates[36] -= (fwd_rates[258] - rev_rates[252]);
  //sp 38
  sp_rates[37] += (fwd_rates[258] - rev_rates[252]);

  //rxn 259
  //sp 37
  sp_rates[36] -= (fwd_rates[259] - rev_rates[253]);
  //sp 7
  sp_rates[6] += (fwd_rates[259] - rev_rates[253]);
  //sp 16
  sp_rates[15] += (fwd_rates[259] - rev_rates[253]);

  //rxn 260
  //sp 0
  sp_rates[0] += (fwd_rates[260] - rev_rates[254]);
  //sp 3
  sp_rates[3] += (fwd_rates[260] - rev_rates[254]);
  //sp 37
  sp_rates[36] -= (fwd_rates[260] - rev_rates[254]);

  //rxn 261
  //sp 38
  sp_rates[37] -= (fwd_rates[261] - rev_rates[255]);
  //sp 7
  sp_rates[6] += (fwd_rates[261] - rev_rates[255]);
  //sp 16
  sp_rates[15] += (fwd_rates[261] - rev_rates[255]);

  //rxn 262
  //sp 0
  sp_rates[0] += (fwd_rates[262] - rev_rates[256]);
  //sp 3
  sp_rates[3] += (fwd_rates[262] - rev_rates[256]);
  //sp 38
  sp_rates[37] -= (fwd_rates[262] - rev_rates[256]);

  //rxn 263
  //sp 33
  sp_rates[32] = (fwd_rates[263] - rev_rates[257]) * pres_mod[30];
  //sp 10
  sp_rates[9] -= (fwd_rates[263] - rev_rates[257]) * pres_mod[30];
  //sp 3
  sp_rates[3] -= (fwd_rates[263] - rev_rates[257]) * pres_mod[30];

  //rxn 264
  //sp 10
  sp_rates[9] -= (fwd_rates[264] - rev_rates[258]) * pres_mod[31];
  //sp 3
  sp_rates[3] -= (fwd_rates[264] - rev_rates[258]) * pres_mod[31];
  //sp 34
  sp_rates[33] = (fwd_rates[264] - rev_rates[258]) * pres_mod[31];

  //rxn 265
  //sp 10
  sp_rates[9] -= (fwd_rates[265] - rev_rates[259]) * pres_mod[32];
  //sp 35
  sp_rates[34] = (fwd_rates[265] - rev_rates[259]) * pres_mod[32];
  //sp 3
  sp_rates[3] -= (fwd_rates[265] - rev_rates[259]) * pres_mod[32];

  //rxn 266
  //sp 33
  sp_rates[32] -= (fwd_rates[266] - rev_rates[260]) * pres_mod[33];
  //sp 34
  sp_rates[33] += (fwd_rates[266] - rev_rates[260]) * pres_mod[33];

  //rxn 267
  //sp 33
  sp_rates[32] -= (fwd_rates[267] - rev_rates[261]);
  //sp 4
  (*dy_N) += (fwd_rates[267] - rev_rates[261]);
  //sp 12
  sp_rates[11] += (fwd_rates[267] - rev_rates[261]);

  //rxn 268
  //sp 34
  sp_rates[33] -= (fwd_rates[268] - rev_rates[262]);
  //sp 4
  (*dy_N) += (fwd_rates[268] - rev_rates[262]);
  //sp 12
  sp_rates[11] += (fwd_rates[268] - rev_rates[262]);

  //rxn 269
  //sp 35
  sp_rates[34] -= (fwd_rates[269] - rev_rates[263]);
  //sp 12
  sp_rates[11] += (fwd_rates[269] - rev_rates[263]);
  //sp 4
  (*dy_N) += (fwd_rates[269] - rev_rates[263]);

  //rxn 270
  //sp 34
  sp_rates[33] -= (fwd_rates[270] - rev_rates[264]);
  //sp 12
  sp_rates[11] -= (fwd_rates[270] - rev_rates[264]);
  //sp 36
  sp_rates[35] += (fwd_rates[270] - rev_rates[264]);

  //rxn 271
  //sp 34
  sp_rates[33] -= (fwd_rates[271] - rev_rates[265]);
  //sp 12
  sp_rates[11] -= (fwd_rates[271] - rev_rates[265]);
  //sp 37
  sp_rates[36] += (fwd_rates[271] - rev_rates[265]);

  //rxn 272
  //sp 34
  sp_rates[33] -= (fwd_rates[272] - rev_rates[266]);
  //sp 12
  sp_rates[11] -= (fwd_rates[272] - rev_rates[266]);
  //sp 38
  sp_rates[37] += (fwd_rates[272] - rev_rates[266]);

  //rxn 273
  //sp 34
  sp_rates[33] -= (fwd_rates[273] - rev_rates[267]);
  //sp 12
  sp_rates[11] -= (fwd_rates[273] - rev_rates[267]);
  //sp 7
  sp_rates[6] += (fwd_rates[273] - rev_rates[267]);
  //sp 16
  sp_rates[15] += (fwd_rates[273] - rev_rates[267]);

  //rxn 274
  //sp 33
  sp_rates[32] += (fwd_rates[274] - rev_rates[268]);
  //sp 34
  sp_rates[33] -= (fwd_rates[274] - rev_rates[268]);

  //rxn 275
  //sp 0
  sp_rates[0] += (fwd_rates[275] - rev_rates[269]);
  //sp 34
  sp_rates[33] -= (fwd_rates[275] - rev_rates[269]);
  //sp 12
  sp_rates[11] -= (fwd_rates[275] - rev_rates[269]);
  //sp 3
  sp_rates[3] += (fwd_rates[275] - rev_rates[269]);

  //rxn 276
  //sp 33
  sp_rates[32] -= (fwd_rates[276] - rev_rates[270]);
  //sp 36
  sp_rates[35] += (fwd_rates[276] - rev_rates[270]);
  //sp 12
  sp_rates[11] -= (fwd_rates[276] - rev_rates[270]);

  //rxn 277
  //sp 33
  sp_rates[32] -= (fwd_rates[277] - rev_rates[271]);
  //sp 12
  sp_rates[11] -= (fwd_rates[277] - rev_rates[271]);
  //sp 37
  sp_rates[36] += (fwd_rates[277] - rev_rates[271]);

  //rxn 278
  //sp 33
  sp_rates[32] -= (fwd_rates[278] - rev_rates[272]);
  //sp 12
  sp_rates[11] -= (fwd_rates[278] - rev_rates[272]);
  //sp 38
  sp_rates[37] += (fwd_rates[278] - rev_rates[272]);

  //rxn 279
  //sp 33
  sp_rates[32] -= (fwd_rates[279] - rev_rates[273]);
  //sp 12
  sp_rates[11] -= (fwd_rates[279] - rev_rates[273]);
  //sp 7
  sp_rates[6] += (fwd_rates[279] - rev_rates[273]);
  //sp 16
  sp_rates[15] += (fwd_rates[279] - rev_rates[273]);

  //rxn 280
  //sp 33
  sp_rates[32] -= (fwd_rates[280] - rev_rates[274]);
  //sp 0
  sp_rates[0] += (fwd_rates[280] - rev_rates[274]);
  //sp 12
  sp_rates[11] -= (fwd_rates[280] - rev_rates[274]);
  //sp 3
  sp_rates[3] += (fwd_rates[280] - rev_rates[274]);

  //rxn 281
  //sp 33
  sp_rates[32] -= (fwd_rates[281] - rev_rates[275]);
  //sp 3
  sp_rates[3] += (fwd_rates[281] - rev_rates[275]);
  //sp 5
  sp_rates[4] -= (fwd_rates[281] - rev_rates[275]);
  //sp 13
  sp_rates[12] += (fwd_rates[281] - rev_rates[275]);

  //rxn 282
  //sp 33
  sp_rates[32] -= (fwd_rates[282] - rev_rates[276]);
  //sp 2
  sp_rates[2] += (fwd_rates[282] - rev_rates[276]);
  //sp 5
  sp_rates[4] -= (fwd_rates[282] - rev_rates[276]);
  //sp 31
  sp_rates[30] += (fwd_rates[282] - rev_rates[276]);

  //rxn 283
  //sp 33
  sp_rates[32] -= (fwd_rates[283] - rev_rates[277]);
  //sp 34
  sp_rates[33] += (fwd_rates[283] - rev_rates[277]);

  //rxn 284
  //sp 34
  sp_rates[33] -= (fwd_rates[284] - rev_rates[278]);
  //sp 3
  sp_rates[3] += (fwd_rates[284] - rev_rates[278]);
  //sp 5
  sp_rates[4] -= (fwd_rates[284] - rev_rates[278]);
  //sp 13
  sp_rates[12] += (fwd_rates[284] - rev_rates[278]);

  //rxn 285
  //sp 34
  sp_rates[33] -= (fwd_rates[285] - rev_rates[279]);
  //sp 2
  sp_rates[2] += (fwd_rates[285] - rev_rates[279]);
  //sp 5
  sp_rates[4] -= (fwd_rates[285] - rev_rates[279]);
  //sp 31
  sp_rates[30] += (fwd_rates[285] - rev_rates[279]);

  //rxn 286
  //sp 25
  sp_rates[24] += (fwd_rates[286] - rev_rates[280]);
  //sp 35
  sp_rates[34] -= (fwd_rates[286] - rev_rates[280]);
  //sp 5
  sp_rates[4] -= (fwd_rates[286] - rev_rates[280]);
  //sp 7
  sp_rates[6] += (fwd_rates[286] - rev_rates[280]);

  //rxn 287
  //sp 35
  sp_rates[34] -= (fwd_rates[287] - rev_rates[281]);
  //sp 3
  sp_rates[3] += (fwd_rates[287] - rev_rates[281]);
  //sp 5
  sp_rates[4] -= (fwd_rates[287] - rev_rates[281]);
  //sp 13
  sp_rates[12] += (fwd_rates[287] - rev_rates[281]);

  //rxn 288
  //sp 33
  sp_rates[32] -= (fwd_rates[288] - rev_rates[282]);
  //sp 10
  sp_rates[9] -= (fwd_rates[288] - rev_rates[282]);
  //sp 2
  sp_rates[2] += (fwd_rates[288] - rev_rates[282]);
  //sp 16
  sp_rates[15] += (fwd_rates[288] - rev_rates[282]);

  //rxn 289
  //sp 33
  sp_rates[32] -= (fwd_rates[289] - rev_rates[283]) * pres_mod[34];
  //sp 4
  (*dy_N) += (fwd_rates[289] - rev_rates[283]) * pres_mod[34];
  //sp 12
  sp_rates[11] += (fwd_rates[289] - rev_rates[283]) * pres_mod[34];

  //rxn 290
  //sp 33
  sp_rates[32] -= (fwd_rates[290] - rev_rates[284]);
  //sp 10
  sp_rates[9] -= (fwd_rates[290] - rev_rates[284]);
  //sp 4
  (*dy_N) += (fwd_rates[290] - rev_rates[284]);
  //sp 0
  sp_rates[0] += (fwd_rates[290] - rev_rates[284]);

  //rxn 291
  //sp 33
  sp_rates[32] += (fwd_rates[291] - rev_rates[285]);
  //sp 34
  sp_rates[33] -= (fwd_rates[291] - rev_rates[285]);

  //rxn 292
  //sp 10
  sp_rates[9] -= (fwd_rates[292] - rev_rates[286]);
  //sp 2
  sp_rates[2] += (fwd_rates[292] - rev_rates[286]);
  //sp 16
  sp_rates[15] += (fwd_rates[292] - rev_rates[286]);
  //sp 34
  sp_rates[33] -= (fwd_rates[292] - rev_rates[286]);

  //rxn 293
  //sp 34
  sp_rates[33] -= (fwd_rates[293] - rev_rates[287]) * pres_mod[35];
  //sp 4
  (*dy_N) += (fwd_rates[293] - rev_rates[287]) * pres_mod[35];
  //sp 12
  sp_rates[11] += (fwd_rates[293] - rev_rates[287]) * pres_mod[35];

  //rxn 294
  //sp 0
  sp_rates[0] += (fwd_rates[294] - rev_rates[288]);
  //sp 10
  sp_rates[9] -= (fwd_rates[294] - rev_rates[288]);
  //sp 4
  (*dy_N) += (fwd_rates[294] - rev_rates[288]);
  //sp 34
  sp_rates[33] -= (fwd_rates[294] - rev_rates[288]);

  //rxn 295
  //sp 10
  sp_rates[9] -= (fwd_rates[295] - rev_rates[289]);
  //sp 35
  sp_rates[34] -= (fwd_rates[295] - rev_rates[289]);
  //sp 2
  sp_rates[2] += (fwd_rates[295] - rev_rates[289]);
  //sp 16
  sp_rates[15] += (fwd_rates[295] - rev_rates[289]);

  //rxn 296
  //sp 35
  sp_rates[34] -= (fwd_rates[296] - rev_rates[290]) * pres_mod[36];
  //sp 12
  sp_rates[11] += (fwd_rates[296] - rev_rates[290]) * pres_mod[36];
  //sp 4
  (*dy_N) += (fwd_rates[296] - rev_rates[290]) * pres_mod[36];

  //rxn 297
  //sp 0
  sp_rates[0] += (fwd_rates[297] - rev_rates[291]);
  //sp 10
  sp_rates[9] -= (fwd_rates[297] - rev_rates[291]);
  //sp 35
  sp_rates[34] -= (fwd_rates[297] - rev_rates[291]);
  //sp 4
  (*dy_N) += (fwd_rates[297] - rev_rates[291]);

  //rxn 298
  //sp 33
  sp_rates[32] -= (fwd_rates[298] - rev_rates[292]);
  //sp 28
  sp_rates[27] += (fwd_rates[298] - rev_rates[292]);
  //sp 2
  sp_rates[2] += (fwd_rates[298] - rev_rates[292]);
  //sp 12
  sp_rates[11] -= (fwd_rates[298] - rev_rates[292]);

  //rxn 299
  //sp 28
  sp_rates[27] += (fwd_rates[299] - rev_rates[293]);
  //sp 34
  sp_rates[33] -= (fwd_rates[299] - rev_rates[293]);
  //sp 2
  sp_rates[2] += (fwd_rates[299] - rev_rates[293]);
  //sp 12
  sp_rates[11] -= (fwd_rates[299] - rev_rates[293]);

  //rxn 300
  //sp 28
  sp_rates[27] += (fwd_rates[300] - rev_rates[294]);
  //sp 2
  sp_rates[2] += (fwd_rates[300] - rev_rates[294]);
  //sp 35
  sp_rates[34] -= (fwd_rates[300] - rev_rates[294]);
  //sp 12
  sp_rates[11] -= (fwd_rates[300] - rev_rates[294]);

  //rxn 301
  //sp 33
  sp_rates[32] -= (fwd_rates[301] - rev_rates[295]);
  //sp 0
  sp_rates[0] += (fwd_rates[301] - rev_rates[295]);
  //sp 12
  sp_rates[11] -= (fwd_rates[301] - rev_rates[295]);
  //sp 3
  sp_rates[3] += (fwd_rates[301] - rev_rates[295]);

  //rxn 302
  //sp 0
  sp_rates[0] += (fwd_rates[302] - rev_rates[296]);
  //sp 34
  sp_rates[33] -= (fwd_rates[302] - rev_rates[296]);
  //sp 12
  sp_rates[11] -= (fwd_rates[302] - rev_rates[296]);
  //sp 3
  sp_rates[3] += (fwd_rates[302] - rev_rates[296]);

  //rxn 303
  //sp 0
  sp_rates[0] += (fwd_rates[303] - rev_rates[297]);
  //sp 3
  sp_rates[3] += (fwd_rates[303] - rev_rates[297]);
  //sp 35
  sp_rates[34] -= (fwd_rates[303] - rev_rates[297]);
  //sp 12
  sp_rates[11] -= (fwd_rates[303] - rev_rates[297]);

  //rxn 304
  //sp 33
  sp_rates[32] -= (fwd_rates[304] - rev_rates[298]);
  //sp 25
  sp_rates[24] += (fwd_rates[304] - rev_rates[298]);
  //sp 11
  sp_rates[10] -= (fwd_rates[304] - rev_rates[298]);
  //sp 2
  sp_rates[2] += (fwd_rates[304] - rev_rates[298]);

  //rxn 305
  //sp 25
  sp_rates[24] += (fwd_rates[305] - rev_rates[299]);
  //sp 34
  sp_rates[33] -= (fwd_rates[305] - rev_rates[299]);
  //sp 11
  sp_rates[10] -= (fwd_rates[305] - rev_rates[299]);
  //sp 2
  sp_rates[2] += (fwd_rates[305] - rev_rates[299]);

  //rxn 306
  //sp 25
  sp_rates[24] += (fwd_rates[306] - rev_rates[300]);
  //sp 11
  sp_rates[10] -= (fwd_rates[306] - rev_rates[300]);
  //sp 2
  sp_rates[2] += (fwd_rates[306] - rev_rates[300]);
  //sp 35
  sp_rates[34] -= (fwd_rates[306] - rev_rates[300]);

  //rxn 307
  //sp 33
  sp_rates[32] -= (fwd_rates[307] - rev_rates[301]);
  //sp 12
  sp_rates[11] += (fwd_rates[307] - rev_rates[301]);
  //sp 11
  sp_rates[10] -= (fwd_rates[307] - rev_rates[301]);
  //sp 3
  sp_rates[3] += (fwd_rates[307] - rev_rates[301]);

  //rxn 308
  //sp 12
  sp_rates[11] += (fwd_rates[308] - rev_rates[302]);
  //sp 34
  sp_rates[33] -= (fwd_rates[308] - rev_rates[302]);
  //sp 11
  sp_rates[10] -= (fwd_rates[308] - rev_rates[302]);
  //sp 3
  sp_rates[3] += (fwd_rates[308] - rev_rates[302]);

  //rxn 309
  //sp 12
  sp_rates[11] += (fwd_rates[309] - rev_rates[303]);
  //sp 11
  sp_rates[10] -= (fwd_rates[309] - rev_rates[303]);
  //sp 3
  sp_rates[3] += (fwd_rates[309] - rev_rates[303]);
  //sp 35
  sp_rates[34] -= (fwd_rates[309] - rev_rates[303]);

  //rxn 310
  //sp 33
  sp_rates[32] -= (fwd_rates[310] - rev_rates[304]);
  //sp 35
  sp_rates[34] += (fwd_rates[310] - rev_rates[304]);

  //rxn 311
  //sp 34
  sp_rates[33] -= (fwd_rates[311] - rev_rates[305]);
  //sp 35
  sp_rates[34] += (fwd_rates[311] - rev_rates[305]);

  //sp 8
  sp_rates[7] = 0.0;
  //sp 9
  sp_rates[8] = 0.0;
  //sp 39
  sp_rates[38] = 0.0;
  //sp 40
  sp_rates[39] = 0.0;
  //sp 41
  sp_rates[40] = 0.0;
  //sp 42
  sp_rates[41] = 0.0;
} // end eval_spec_rates

