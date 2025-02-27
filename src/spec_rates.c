#include "header.h"
#include "rates.h"

void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  //rxn 0
  //sp 0
  sp_rates[0] = -fwd_rates[0];
  //sp 2
  sp_rates[2] = fwd_rates[0];
  //sp 3
  sp_rates[3] = -fwd_rates[0];
  //sp 4
  sp_rates[4] = fwd_rates[0];

  //rxn 1
  //sp 0
  sp_rates[0] += fwd_rates[1];
  //sp 2
  sp_rates[2] -= fwd_rates[1];
  //sp 3
  sp_rates[3] += fwd_rates[1];
  //sp 4
  sp_rates[4] -= fwd_rates[1];

  //rxn 2
  //sp 0
  sp_rates[0] += fwd_rates[2];
  //sp 1
  sp_rates[1] = -fwd_rates[2];
  //sp 2
  sp_rates[2] -= fwd_rates[2];
  //sp 4
  sp_rates[4] += fwd_rates[2];

  //rxn 3
  //sp 0
  sp_rates[0] -= fwd_rates[3];
  //sp 1
  sp_rates[1] += fwd_rates[3];
  //sp 2
  sp_rates[2] += fwd_rates[3];
  //sp 4
  sp_rates[4] -= fwd_rates[3];

  //rxn 4
  //sp 0
  sp_rates[0] += fwd_rates[4];
  //sp 1
  sp_rates[1] -= fwd_rates[4];
  //sp 4
  sp_rates[4] -= fwd_rates[4];
  //sp 5
  sp_rates[5] = fwd_rates[4];

  //rxn 5
  //sp 0
  sp_rates[0] -= fwd_rates[5];
  //sp 1
  sp_rates[1] += fwd_rates[5];
  //sp 4
  sp_rates[4] += fwd_rates[5];
  //sp 5
  sp_rates[5] -= fwd_rates[5];

  //rxn 6
  //sp 2
  sp_rates[2] -= fwd_rates[6];
  //sp 4
  sp_rates[4] += 2.0 * fwd_rates[6];
  //sp 5
  sp_rates[5] -= fwd_rates[6];

  //rxn 7
  //sp 2
  sp_rates[2] += fwd_rates[7];
  //sp 4
  sp_rates[4] -= 2.0 * fwd_rates[7];
  //sp 5
  sp_rates[5] += fwd_rates[7];

  //rxn 8
  //sp 0
  sp_rates[0] += 2.0 * fwd_rates[8] * pres_mod[0];
  //sp 1
  sp_rates[1] -= fwd_rates[8] * pres_mod[0];

  //rxn 9
  //sp 0
  sp_rates[0] -= 2.0 * fwd_rates[9] * pres_mod[1];
  //sp 1
  sp_rates[1] += fwd_rates[9] * pres_mod[1];

  //rxn 10
  //sp 2
  sp_rates[2] += 2.0 * fwd_rates[10] * pres_mod[2];
  //sp 3
  sp_rates[3] -= fwd_rates[10] * pres_mod[2];

  //rxn 11
  //sp 2
  sp_rates[2] -= 2.0 * fwd_rates[11] * pres_mod[3];
  //sp 3
  sp_rates[3] += fwd_rates[11] * pres_mod[3];

  //rxn 12
  //sp 0
  sp_rates[0] += fwd_rates[12] * pres_mod[4];
  //sp 2
  sp_rates[2] += fwd_rates[12] * pres_mod[4];
  //sp 4
  sp_rates[4] -= fwd_rates[12] * pres_mod[4];

  //rxn 13
  //sp 0
  sp_rates[0] -= fwd_rates[13] * pres_mod[5];
  //sp 2
  sp_rates[2] -= fwd_rates[13] * pres_mod[5];
  //sp 4
  sp_rates[4] += fwd_rates[13] * pres_mod[5];

  //rxn 14
  //sp 0
  sp_rates[0] += fwd_rates[14] * pres_mod[6];
  //sp 4
  sp_rates[4] += fwd_rates[14] * pres_mod[6];
  //sp 5
  sp_rates[5] -= fwd_rates[14] * pres_mod[6];

  //rxn 15
  //sp 0
  sp_rates[0] -= fwd_rates[15] * pres_mod[7];
  //sp 4
  sp_rates[4] -= fwd_rates[15] * pres_mod[7];
  //sp 5
  sp_rates[5] += fwd_rates[15] * pres_mod[7];

  //rxn 16
  //sp 0
  sp_rates[0] -= (fwd_rates[16] - rev_rates[0]) * pres_mod[8];
  //sp 3
  sp_rates[3] -= (fwd_rates[16] - rev_rates[0]) * pres_mod[8];
  //sp 7
  sp_rates[6] = (fwd_rates[16] - rev_rates[0]) * pres_mod[8];

  //rxn 17
  //sp 0
  sp_rates[0] -= fwd_rates[17];
  //sp 1
  sp_rates[1] += fwd_rates[17];
  //sp 3
  sp_rates[3] += fwd_rates[17];
  //sp 7
  sp_rates[6] -= fwd_rates[17];

  //rxn 18
  //sp 0
  sp_rates[0] += fwd_rates[18];
  //sp 1
  sp_rates[1] -= fwd_rates[18];
  //sp 3
  sp_rates[3] -= fwd_rates[18];
  //sp 7
  sp_rates[6] += fwd_rates[18];

  //rxn 19
  //sp 0
  sp_rates[0] -= fwd_rates[19];
  //sp 4
  sp_rates[4] += 2.0 * fwd_rates[19];
  //sp 7
  sp_rates[6] -= fwd_rates[19];

  //rxn 20
  //sp 0
  sp_rates[0] += fwd_rates[20];
  //sp 4
  sp_rates[4] -= 2.0 * fwd_rates[20];
  //sp 7
  sp_rates[6] += fwd_rates[20];

  //rxn 21
  //sp 2
  sp_rates[2] -= fwd_rates[21];
  //sp 3
  sp_rates[3] += fwd_rates[21];
  //sp 4
  sp_rates[4] += fwd_rates[21];
  //sp 7
  sp_rates[6] -= fwd_rates[21];

  //rxn 22
  //sp 2
  sp_rates[2] += fwd_rates[22];
  //sp 3
  sp_rates[3] -= fwd_rates[22];
  //sp 4
  sp_rates[4] -= fwd_rates[22];
  //sp 7
  sp_rates[6] += fwd_rates[22];

  //rxn 23
  //sp 3
  sp_rates[3] += fwd_rates[23];
  //sp 4
  sp_rates[4] -= fwd_rates[23];
  //sp 5
  sp_rates[5] += fwd_rates[23];
  //sp 7
  sp_rates[6] -= fwd_rates[23];

  //rxn 24
  //sp 3
  sp_rates[3] -= fwd_rates[24];
  //sp 4
  sp_rates[4] += fwd_rates[24];
  //sp 5
  sp_rates[5] -= fwd_rates[24];
  //sp 7
  sp_rates[6] += fwd_rates[24];

  //rxn 25
  //sp 3
  sp_rates[3] -= fwd_rates[25];
  //sp 7
  sp_rates[6] += 2.0 * fwd_rates[25];
  //sp 8
  sp_rates[7] = -fwd_rates[25];

  //rxn 26
  //sp 3
  sp_rates[3] += fwd_rates[26];
  //sp 7
  sp_rates[6] -= 2.0 * fwd_rates[26];
  //sp 8
  sp_rates[7] += fwd_rates[26];

  //rxn 27
  //sp 3
  sp_rates[3] -= fwd_rates[27];
  //sp 7
  sp_rates[6] += 2.0 * fwd_rates[27];
  //sp 8
  sp_rates[7] -= fwd_rates[27];

  //rxn 28
  //sp 3
  sp_rates[3] += fwd_rates[28];
  //sp 7
  sp_rates[6] -= 2.0 * fwd_rates[28];
  //sp 8
  sp_rates[7] += fwd_rates[28];

  //rxn 29
  //sp 4
  sp_rates[4] += 2.0 * (fwd_rates[29] - rev_rates[1]) * pres_mod[9];
  //sp 8
  sp_rates[7] -= (fwd_rates[29] - rev_rates[1]) * pres_mod[9];

  //rxn 30
  //sp 0
  sp_rates[0] -= fwd_rates[30];
  //sp 4
  sp_rates[4] += fwd_rates[30];
  //sp 5
  sp_rates[5] += fwd_rates[30];
  //sp 8
  sp_rates[7] -= fwd_rates[30];

  //rxn 31
  //sp 0
  sp_rates[0] += fwd_rates[31];
  //sp 4
  sp_rates[4] -= fwd_rates[31];
  //sp 5
  sp_rates[5] -= fwd_rates[31];
  //sp 8
  sp_rates[7] += fwd_rates[31];

  //rxn 32
  //sp 0
  sp_rates[0] -= fwd_rates[32];
  //sp 1
  sp_rates[1] += fwd_rates[32];
  //sp 7
  sp_rates[6] += fwd_rates[32];
  //sp 8
  sp_rates[7] -= fwd_rates[32];

  //rxn 33
  //sp 0
  sp_rates[0] += fwd_rates[33];
  //sp 1
  sp_rates[1] -= fwd_rates[33];
  //sp 7
  sp_rates[6] -= fwd_rates[33];
  //sp 8
  sp_rates[7] += fwd_rates[33];

  //rxn 34
  //sp 2
  sp_rates[2] -= fwd_rates[34];
  //sp 4
  sp_rates[4] += fwd_rates[34];
  //sp 7
  sp_rates[6] += fwd_rates[34];
  //sp 8
  sp_rates[7] -= fwd_rates[34];

  //rxn 35
  //sp 2
  sp_rates[2] += fwd_rates[35];
  //sp 4
  sp_rates[4] -= fwd_rates[35];
  //sp 7
  sp_rates[6] -= fwd_rates[35];
  //sp 8
  sp_rates[7] += fwd_rates[35];

  //rxn 36
  //sp 4
  sp_rates[4] -= fwd_rates[36];
  //sp 5
  sp_rates[5] += fwd_rates[36];
  //sp 7
  sp_rates[6] += fwd_rates[36];
  //sp 8
  sp_rates[7] -= fwd_rates[36];

  //rxn 37
  //sp 4
  sp_rates[4] += fwd_rates[37];
  //sp 5
  sp_rates[5] -= fwd_rates[37];
  //sp 7
  sp_rates[6] -= fwd_rates[37];
  //sp 8
  sp_rates[7] += fwd_rates[37];

  //rxn 38
  //sp 4
  sp_rates[4] -= fwd_rates[38];
  //sp 5
  sp_rates[5] += fwd_rates[38];
  //sp 7
  sp_rates[6] += fwd_rates[38];
  //sp 8
  sp_rates[7] -= fwd_rates[38];

  //rxn 39
  //sp 4
  sp_rates[4] += fwd_rates[39];
  //sp 5
  sp_rates[5] -= fwd_rates[39];
  //sp 7
  sp_rates[6] -= fwd_rates[39];
  //sp 8
  sp_rates[7] += fwd_rates[39];

  //sp 9
  sp_rates[8] = 0.0;
  //sp 6
  (*dy_N) = 0.0;
} // end eval_spec_rates

