#include "header.cuh"
#include "rates.cuh"

__device__ void eval_spec_rates (const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, double * __restrict__ sp_rates, double * __restrict__ dy_N) {
  extern volatile __shared__ double shared_temp[];
  //rxn 0
  //sp 4
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(0)] - rev_rates[INDEX(0)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] = -(fwd_rates[INDEX(0)] - rev_rates[INDEX(0)]);
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] = (fwd_rates[INDEX(0)] - rev_rates[INDEX(0)]);
  //sp 7
  shared_temp[threadIdx.x + 3 * blockDim.x] = (fwd_rates[INDEX(0)] - rev_rates[INDEX(0)]);

  //rxn 1
  //sp 3
  sp_rates[INDEX(2)] = -(fwd_rates[INDEX(1)] - rev_rates[INDEX(1)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(1)] - rev_rates[INDEX(1)]);
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(1)] - rev_rates[INDEX(1)]);
  //sp 7
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(1)] - rev_rates[INDEX(1)]);

  //rxn 2
  //sp 3
  sp_rates[INDEX(2)] -= (fwd_rates[INDEX(2)] - rev_rates[INDEX(2)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(2)] - rev_rates[INDEX(2)]);
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(2)] - rev_rates[INDEX(2)]);
  //sp 7
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(2)] - rev_rates[INDEX(2)]);

  //rxn 3
  sp_rates[INDEX(3)] = shared_temp[threadIdx.x];
  //sp 9
  shared_temp[threadIdx.x] = (fwd_rates[INDEX(3)] - rev_rates[INDEX(3)]);
  //sp 3
  sp_rates[INDEX(2)] -= (fwd_rates[INDEX(3)] - rev_rates[INDEX(3)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(3)] - rev_rates[INDEX(3)]);
  //sp 7
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(3)] - rev_rates[INDEX(3)]);

  //rxn 4
  //sp 9
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(4)] - rev_rates[INDEX(4)]);
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(4)] - rev_rates[INDEX(4)]);
  //sp 7
  shared_temp[threadIdx.x + 3 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(4)] - rev_rates[INDEX(4)]);

  //rxn 5
  //sp 9
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(5)] - rev_rates[INDEX(5)]);
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(5)] - rev_rates[INDEX(5)]);
  //sp 7
  shared_temp[threadIdx.x + 3 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(5)] - rev_rates[INDEX(5)]);

  //rxn 6
  //sp 3
  sp_rates[INDEX(2)] -= (fwd_rates[INDEX(6)] - rev_rates[INDEX(6)]) * pres_mod[INDEX(0)];
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += 2.0 * (fwd_rates[INDEX(6)] - rev_rates[INDEX(6)]) * pres_mod[INDEX(0)];

  //rxn 7
  //sp 3
  sp_rates[INDEX(2)] -= (fwd_rates[INDEX(7)] - rev_rates[INDEX(7)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += 2.0 * (fwd_rates[INDEX(7)] - rev_rates[INDEX(7)]);

  //rxn 8
  //sp 3
  sp_rates[INDEX(2)] -= (fwd_rates[INDEX(8)] - rev_rates[INDEX(8)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += 2.0 * (fwd_rates[INDEX(8)] - rev_rates[INDEX(8)]);

  //rxn 9
  sp_rates[INDEX(6)] = shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] = (fwd_rates[INDEX(9)] - rev_rates[INDEX(9)]) * pres_mod[INDEX(1)];
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(9)] - rev_rates[INDEX(9)]) * pres_mod[INDEX(1)];

  //rxn 10
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(10)] - rev_rates[INDEX(10)]);
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(10)] - rev_rates[INDEX(10)]);

  //rxn 11
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(11)] - rev_rates[INDEX(11)]);
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(11)] - rev_rates[INDEX(11)]);

  //rxn 12
  sp_rates[INDEX(8)] = shared_temp[threadIdx.x];
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(12)] - rev_rates[INDEX(12)]) * pres_mod[INDEX(2)];
  //sp 6
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(12)] - rev_rates[INDEX(12)]) * pres_mod[INDEX(2)];
  //sp 7
  shared_temp[threadIdx.x] = (fwd_rates[INDEX(12)] - rev_rates[INDEX(12)]) * pres_mod[INDEX(2)];

  //rxn 13
  //sp 9
  sp_rates[INDEX(8)] -= (fwd_rates[INDEX(13)] - rev_rates[INDEX(13)]) * pres_mod[INDEX(3)];
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(13)] - rev_rates[INDEX(13)]) * pres_mod[INDEX(3)];
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(13)] - rev_rates[INDEX(13)]) * pres_mod[INDEX(3)];

  //rxn 14
  //sp 9
  sp_rates[INDEX(8)] -= (fwd_rates[INDEX(14)] - rev_rates[INDEX(14)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(14)] - rev_rates[INDEX(14)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(14)] - rev_rates[INDEX(14)]);

  //rxn 15
  sp_rates[INDEX(5)] = shared_temp[threadIdx.x + 1 * blockDim.x];
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(15)] - rev_rates[INDEX(15)]) * pres_mod[INDEX(4)];
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(15)] - rev_rates[INDEX(15)]) * pres_mod[INDEX(4)];
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] = (fwd_rates[INDEX(15)] - rev_rates[INDEX(15)]) * pres_mod[INDEX(4)];

  //rxn 16
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(16)] - rev_rates[INDEX(16)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(16)] - rev_rates[INDEX(16)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(16)] - rev_rates[INDEX(16)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(16)] - rev_rates[INDEX(16)]);

  //rxn 17
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(17)] - rev_rates[INDEX(17)]);
  //sp 7
  shared_temp[threadIdx.x] += 2.0 * (fwd_rates[INDEX(17)] - rev_rates[INDEX(17)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(17)] - rev_rates[INDEX(17)]);

  //rxn 18
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(18)] - rev_rates[INDEX(18)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(18)] - rev_rates[INDEX(18)]);
  //sp 6
  sp_rates[INDEX(5)] += (fwd_rates[INDEX(18)] - rev_rates[INDEX(18)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(18)] - rev_rates[INDEX(18)]);

  //rxn 19
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(19)] - rev_rates[INDEX(19)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(19)] - rev_rates[INDEX(19)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(19)] - rev_rates[INDEX(19)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(19)] - rev_rates[INDEX(19)]);

  //rxn 20
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(20)] - rev_rates[INDEX(20)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(20)] - rev_rates[INDEX(20)]);
  //sp 7
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(20)] - rev_rates[INDEX(20)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(20)] - rev_rates[INDEX(20)]);

  //rxn 21
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(21)] - rev_rates[INDEX(21)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(21)] - rev_rates[INDEX(21)]);
  //sp 7
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(21)] - rev_rates[INDEX(21)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(21)] - rev_rates[INDEX(21)]);

  //rxn 22
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(22)] - rev_rates[INDEX(22)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(22)] - rev_rates[INDEX(22)]);

  //rxn 23
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(23)] - rev_rates[INDEX(23)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(23)] - rev_rates[INDEX(23)]);
  //sp 7
  shared_temp[threadIdx.x] += 2.0 * (fwd_rates[INDEX(23)] - rev_rates[INDEX(23)]);

  //rxn 24
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(24)] - rev_rates[INDEX(24)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(24)] - rev_rates[INDEX(24)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(24)] - rev_rates[INDEX(24)]);

  //rxn 25
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(25)] - rev_rates[INDEX(25)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(25)] - rev_rates[INDEX(25)]);
  //sp 7
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(25)] - rev_rates[INDEX(25)]);

  //rxn 26
  //sp 10
  sp_rates[INDEX(9)] = (fwd_rates[INDEX(26)] - rev_rates[INDEX(26)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(26)] - rev_rates[INDEX(26)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(26)] - rev_rates[INDEX(26)]);

  //rxn 27
  sp_rates[INDEX(4)] = shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 10
  shared_temp[threadIdx.x + 2 * blockDim.x] = (fwd_rates[INDEX(27)] - rev_rates[INDEX(27)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(27)] - rev_rates[INDEX(27)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(27)] - rev_rates[INDEX(27)]);

  //rxn 28
  //sp 10
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(28)] - rev_rates[INDEX(28)]) * pres_mod[INDEX(5)];
  //sp 7
  shared_temp[threadIdx.x] += 2.0 * (fwd_rates[INDEX(28)] - rev_rates[INDEX(28)]) * pres_mod[INDEX(5)];

  //rxn 29
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(29)] - rev_rates[INDEX(29)]);
  //sp 10
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(29)] - rev_rates[INDEX(29)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(29)] - rev_rates[INDEX(29)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(29)] - rev_rates[INDEX(29)]);

  //rxn 30
  //sp 10
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(30)] - rev_rates[INDEX(30)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(30)] - rev_rates[INDEX(30)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(30)] - rev_rates[INDEX(30)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(30)] - rev_rates[INDEX(30)]);

  //rxn 31
  //sp 10
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(31)] - rev_rates[INDEX(31)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(31)] - rev_rates[INDEX(31)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(31)] - rev_rates[INDEX(31)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(31)] - rev_rates[INDEX(31)]);

  //rxn 32
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(32)] - rev_rates[INDEX(32)]);
  //sp 10
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(32)] - rev_rates[INDEX(32)]);
  //sp 7
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(32)] - rev_rates[INDEX(32)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(32)] - rev_rates[INDEX(32)]);

  //rxn 33
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(33)] - rev_rates[INDEX(33)]);
  //sp 10
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(33)] - rev_rates[INDEX(33)]);
  //sp 7
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(33)] - rev_rates[INDEX(33)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(33)] - rev_rates[INDEX(33)]);

  //rxn 34
  sp_rates[INDEX(3)] += shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(34)] - rev_rates[INDEX(34)]) * pres_mod[INDEX(6)];
  //sp 13
  sp_rates[INDEX(12)] = -(fwd_rates[INDEX(34)] - rev_rates[INDEX(34)]) * pres_mod[INDEX(6)];
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] = (fwd_rates[INDEX(34)] - rev_rates[INDEX(34)]) * pres_mod[INDEX(6)];

  //rxn 35
  sp_rates[INDEX(6)] += shared_temp[threadIdx.x];
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(35)] - rev_rates[INDEX(35)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(35)] - rev_rates[INDEX(35)]);
  //sp 13
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(35)] - rev_rates[INDEX(35)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(35)] - rev_rates[INDEX(35)]);

  //rxn 36
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(36)] - rev_rates[INDEX(36)]);
  //sp 13
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(36)] - rev_rates[INDEX(36)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(36)] - rev_rates[INDEX(36)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(36)] - rev_rates[INDEX(36)]);

  //rxn 37
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(37)] - rev_rates[INDEX(37)]);
  //sp 13
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(37)] - rev_rates[INDEX(37)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(37)] - rev_rates[INDEX(37)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(37)] - rev_rates[INDEX(37)]);

  //rxn 38
  //sp 10
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(38)] - rev_rates[INDEX(38)]);
  //sp 13
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(38)] - rev_rates[INDEX(38)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(38)] - rev_rates[INDEX(38)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(38)] - rev_rates[INDEX(38)]);

  //rxn 39
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(39)] - rev_rates[INDEX(39)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(39)] - rev_rates[INDEX(39)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(39)] - rev_rates[INDEX(39)]);
  //sp 15
  sp_rates[INDEX(14)] = (fwd_rates[INDEX(39)] - rev_rates[INDEX(39)]);

  //rxn 40
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(40)] - rev_rates[INDEX(40)]);
  //sp 19
  sp_rates[INDEX(18)] = (fwd_rates[INDEX(40)] - rev_rates[INDEX(40)]);
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(40)] - rev_rates[INDEX(40)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(40)] - rev_rates[INDEX(40)]);

  //rxn 41
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(41)] - rev_rates[INDEX(41)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(41)] - rev_rates[INDEX(41)]);
  //sp 15
  sp_rates[INDEX(14)] += (fwd_rates[INDEX(41)] - rev_rates[INDEX(41)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(41)] - rev_rates[INDEX(41)]);

  //rxn 42
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(42)] - rev_rates[INDEX(42)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(42)] - rev_rates[INDEX(42)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(42)] - rev_rates[INDEX(42)]);
  //sp 15
  sp_rates[INDEX(14)] += (fwd_rates[INDEX(42)] - rev_rates[INDEX(42)]);

  //rxn 43
  //sp 4
  sp_rates[INDEX(3)] += (fwd_rates[INDEX(43)] - rev_rates[INDEX(43)]);
  //sp 13
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(43)] - rev_rates[INDEX(43)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(43)] - rev_rates[INDEX(43)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(43)] - rev_rates[INDEX(43)]);

  //rxn 44
  //sp 4
  sp_rates[INDEX(3)] += (fwd_rates[INDEX(44)] - rev_rates[INDEX(44)]);
  //sp 13
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(44)] - rev_rates[INDEX(44)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(44)] - rev_rates[INDEX(44)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(44)] - rev_rates[INDEX(44)]);

  //rxn 45
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(45)] - rev_rates[INDEX(45)]);
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(45)] - rev_rates[INDEX(45)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(45)] - rev_rates[INDEX(45)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(45)] - rev_rates[INDEX(45)]);

  //rxn 46
  //sp 21
  sp_rates[INDEX(20)] = (fwd_rates[INDEX(46)] - rev_rates[INDEX(46)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(46)] - rev_rates[INDEX(46)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(46)] - rev_rates[INDEX(46)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(46)] - rev_rates[INDEX(46)]);

  //rxn 47
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(47)] - rev_rates[INDEX(47)]);
  //sp 20
  sp_rates[INDEX(19)] = (fwd_rates[INDEX(47)] - rev_rates[INDEX(47)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(47)] - rev_rates[INDEX(47)]);
  //sp 8
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(47)] - rev_rates[INDEX(47)]);

  //rxn 48
  //sp 6
  sp_rates[INDEX(5)] += (fwd_rates[INDEX(48)] - rev_rates[INDEX(48)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(48)] - rev_rates[INDEX(48)]);
  //sp 21
  sp_rates[INDEX(20)] += (fwd_rates[INDEX(48)] - rev_rates[INDEX(48)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(48)] - rev_rates[INDEX(48)]);

  //rxn 49
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(49)] - rev_rates[INDEX(49)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(49)] - rev_rates[INDEX(49)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(49)] - rev_rates[INDEX(49)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(49)] - rev_rates[INDEX(49)]);

  //rxn 50
  //sp 13
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(50)] - rev_rates[INDEX(50)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(50)] - rev_rates[INDEX(50)]);
  //sp 15
  sp_rates[INDEX(14)] += (fwd_rates[INDEX(50)] - rev_rates[INDEX(50)]);

  //rxn 51
  //sp 12
  sp_rates[INDEX(11)] = (fwd_rates[INDEX(51)] - rev_rates[INDEX(51)]);
  //sp 13
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(51)] - rev_rates[INDEX(51)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(51)] - rev_rates[INDEX(51)]);
  //sp 15
  sp_rates[INDEX(14)] -= (fwd_rates[INDEX(51)] - rev_rates[INDEX(51)]);

  //rxn 52
  //sp 2
  (*dy_N) = (fwd_rates[INDEX(52)] - rev_rates[INDEX(52)]);
  //sp 12
  sp_rates[INDEX(11)] -= (fwd_rates[INDEX(52)] - rev_rates[INDEX(52)]);
  //sp 5
  sp_rates[INDEX(4)] += 2.0 * (fwd_rates[INDEX(52)] - rev_rates[INDEX(52)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(52)] - rev_rates[INDEX(52)]);

  //rxn 53
  sp_rates[INDEX(9)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] = (fwd_rates[INDEX(53)] - rev_rates[INDEX(53)]);
  //sp 19
  sp_rates[INDEX(18)] -= (fwd_rates[INDEX(53)] - rev_rates[INDEX(53)]);
  //sp 13
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(53)] - rev_rates[INDEX(53)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(53)] - rev_rates[INDEX(53)]);

  //rxn 54
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(54)] - rev_rates[INDEX(54)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(54)] - rev_rates[INDEX(54)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(54)] - rev_rates[INDEX(54)]);
  //sp 16
  sp_rates[INDEX(15)] = (fwd_rates[INDEX(54)] - rev_rates[INDEX(54)]);

  //rxn 55
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(55)] - rev_rates[INDEX(55)]);
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(55)] - rev_rates[INDEX(55)]);
  //sp 2
  (*dy_N) += (fwd_rates[INDEX(55)] - rev_rates[INDEX(55)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(55)] - rev_rates[INDEX(55)]);

  //rxn 56
  sp_rates[INDEX(7)] = shared_temp[threadIdx.x + 1 * blockDim.x];
  //sp 25
  sp_rates[INDEX(24)] = -(fwd_rates[INDEX(56)] - rev_rates[INDEX(56)]);
  //sp 13
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(56)] - rev_rates[INDEX(56)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(56)] - rev_rates[INDEX(56)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] = (fwd_rates[INDEX(56)] - rev_rates[INDEX(56)]);

  //rxn 57
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(57)] - rev_rates[INDEX(57)]);
  //sp 18
  sp_rates[INDEX(17)] = (fwd_rates[INDEX(57)] - rev_rates[INDEX(57)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(57)] - rev_rates[INDEX(57)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(57)] - rev_rates[INDEX(57)]);

  //rxn 58
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(58)] - rev_rates[INDEX(58)]);
  //sp 21
  sp_rates[INDEX(20)] += (fwd_rates[INDEX(58)] - rev_rates[INDEX(58)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(58)] - rev_rates[INDEX(58)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(58)] - rev_rates[INDEX(58)]);

  //rxn 59
  sp_rates[INDEX(12)] += shared_temp[threadIdx.x];
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(59)] - rev_rates[INDEX(59)]) * pres_mod[INDEX(7)];
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(59)] - rev_rates[INDEX(59)]) * pres_mod[INDEX(7)];
  //sp 15
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(59)] - rev_rates[INDEX(59)]) * pres_mod[INDEX(7)];

  //rxn 60
  sp_rates[INDEX(16)] = shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(60)] - rev_rates[INDEX(60)]);
  //sp 12
  sp_rates[INDEX(11)] += (fwd_rates[INDEX(60)] - rev_rates[INDEX(60)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] = -(fwd_rates[INDEX(60)] - rev_rates[INDEX(60)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(60)] - rev_rates[INDEX(60)]);

  //rxn 61
  //sp 17
  sp_rates[INDEX(16)] += (fwd_rates[INDEX(61)] - rev_rates[INDEX(61)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(61)] - rev_rates[INDEX(61)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(61)] - rev_rates[INDEX(61)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(61)] - rev_rates[INDEX(61)]);

  //rxn 62
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(62)] - rev_rates[INDEX(62)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(62)] - rev_rates[INDEX(62)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(62)] - rev_rates[INDEX(62)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(62)] - rev_rates[INDEX(62)]);

  //rxn 63
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(63)] - rev_rates[INDEX(63)]);
  //sp 12
  sp_rates[INDEX(11)] += (fwd_rates[INDEX(63)] - rev_rates[INDEX(63)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(63)] - rev_rates[INDEX(63)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(63)] - rev_rates[INDEX(63)]);

  //rxn 64
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(64)] - rev_rates[INDEX(64)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(64)] - rev_rates[INDEX(64)]);
  //sp 6
  sp_rates[INDEX(5)] += (fwd_rates[INDEX(64)] - rev_rates[INDEX(64)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(64)] - rev_rates[INDEX(64)]);

  //rxn 65
  //sp 17
  sp_rates[INDEX(16)] += (fwd_rates[INDEX(65)] - rev_rates[INDEX(65)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(65)] - rev_rates[INDEX(65)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(65)] - rev_rates[INDEX(65)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(65)] - rev_rates[INDEX(65)]);

  //rxn 66
  //sp 12
  sp_rates[INDEX(11)] += (fwd_rates[INDEX(66)] - rev_rates[INDEX(66)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(66)] - rev_rates[INDEX(66)]);
  //sp 15
  shared_temp[threadIdx.x] -= 2.0 * (fwd_rates[INDEX(66)] - rev_rates[INDEX(66)]);

  //rxn 67
  //sp 2
  (*dy_N) += (fwd_rates[INDEX(67)] - rev_rates[INDEX(67)]);
  //sp 12
  sp_rates[INDEX(11)] -= (fwd_rates[INDEX(67)] - rev_rates[INDEX(67)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(67)] - rev_rates[INDEX(67)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(67)] - rev_rates[INDEX(67)]);

  //rxn 68
  //sp 17
  sp_rates[INDEX(16)] -= (fwd_rates[INDEX(68)] - rev_rates[INDEX(68)]);
  //sp 18
  sp_rates[INDEX(17)] += (fwd_rates[INDEX(68)] - rev_rates[INDEX(68)]);
  //sp 5
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(68)] - rev_rates[INDEX(68)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(68)] - rev_rates[INDEX(68)]);

  //rxn 69
  //sp 17
  sp_rates[INDEX(16)] -= (fwd_rates[INDEX(69)] - rev_rates[INDEX(69)]);
  //sp 2
  (*dy_N) += (fwd_rates[INDEX(69)] - rev_rates[INDEX(69)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(69)] - rev_rates[INDEX(69)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(69)] - rev_rates[INDEX(69)]);

  //rxn 70
  //sp 25
  sp_rates[INDEX(24)] -= (fwd_rates[INDEX(70)] - rev_rates[INDEX(70)]);
  //sp 14
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(70)] - rev_rates[INDEX(70)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(70)] - rev_rates[INDEX(70)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(70)] - rev_rates[INDEX(70)]);

  //rxn 71
  //sp 18
  sp_rates[INDEX(17)] += (fwd_rates[INDEX(71)] - rev_rates[INDEX(71)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(71)] - rev_rates[INDEX(71)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(71)] - rev_rates[INDEX(71)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(71)] - rev_rates[INDEX(71)]);

  //rxn 72
  sp_rates[INDEX(4)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] = (fwd_rates[INDEX(72)] - rev_rates[INDEX(72)]);
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(72)] - rev_rates[INDEX(72)]);
  //sp 15
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(72)] - rev_rates[INDEX(72)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(72)] - rev_rates[INDEX(72)]);

  //rxn 73
  sp_rates[INDEX(13)] = shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(73)] - rev_rates[INDEX(73)]);
  //sp 12
  shared_temp[threadIdx.x + 3 * blockDim.x] = -(fwd_rates[INDEX(73)] - rev_rates[INDEX(73)]);
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(73)] - rev_rates[INDEX(73)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(73)] - rev_rates[INDEX(73)]);

  //rxn 74
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(74)] - rev_rates[INDEX(74)]);
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(74)] - rev_rates[INDEX(74)]);
  //sp 12
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(74)] - rev_rates[INDEX(74)]);
  //sp 6
  sp_rates[INDEX(5)] += (fwd_rates[INDEX(74)] - rev_rates[INDEX(74)]);

  //rxn 75
  sp_rates[INDEX(23)] = shared_temp[threadIdx.x + 1 * blockDim.x];
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(75)] - rev_rates[INDEX(75)]);
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] = (fwd_rates[INDEX(75)] - rev_rates[INDEX(75)]);
  //sp 12
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(75)] - rev_rates[INDEX(75)]);
  //sp 6
  sp_rates[INDEX(5)] += (fwd_rates[INDEX(75)] - rev_rates[INDEX(75)]);

  //rxn 76
  sp_rates[INDEX(14)] += shared_temp[threadIdx.x];
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(76)] - rev_rates[INDEX(76)]);
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(76)] - rev_rates[INDEX(76)]);
  //sp 16
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(76)] - rev_rates[INDEX(76)]);

  //rxn 77
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(77)] - rev_rates[INDEX(77)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(77)] - rev_rates[INDEX(77)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(77)] - rev_rates[INDEX(77)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(77)] - rev_rates[INDEX(77)]);

  //rxn 78
  sp_rates[INDEX(16)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 18
  sp_rates[INDEX(17)] += (fwd_rates[INDEX(78)] - rev_rates[INDEX(78)]);
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(78)] - rev_rates[INDEX(78)]);
  //sp 6
  shared_temp[threadIdx.x + 2 * blockDim.x] = -(fwd_rates[INDEX(78)] - rev_rates[INDEX(78)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(78)] - rev_rates[INDEX(78)]);

  //rxn 79
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(79)] - rev_rates[INDEX(79)]);
  //sp 6
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(79)] - rev_rates[INDEX(79)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(79)] - rev_rates[INDEX(79)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(79)] - rev_rates[INDEX(79)]);

  //rxn 80
  //sp 17
  sp_rates[INDEX(16)] += (fwd_rates[INDEX(80)] - rev_rates[INDEX(80)]);
  //sp 6
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(80)] - rev_rates[INDEX(80)]);
  //sp 15
  sp_rates[INDEX(14)] += (fwd_rates[INDEX(80)] - rev_rates[INDEX(80)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(80)] - rev_rates[INDEX(80)]);

  //rxn 81
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(81)] - rev_rates[INDEX(81)]);
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(81)] - rev_rates[INDEX(81)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(81)] - rev_rates[INDEX(81)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(81)] - rev_rates[INDEX(81)]);

  //rxn 82
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(82)] - rev_rates[INDEX(82)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(82)] - rev_rates[INDEX(82)]);
  //sp 8
  sp_rates[INDEX(7)] += (fwd_rates[INDEX(82)] - rev_rates[INDEX(82)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(82)] - rev_rates[INDEX(82)]);

  //rxn 83
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(83)] - rev_rates[INDEX(83)]);
  //sp 14
  sp_rates[INDEX(13)] += (fwd_rates[INDEX(83)] - rev_rates[INDEX(83)]);
  //sp 15
  sp_rates[INDEX(14)] -= (fwd_rates[INDEX(83)] - rev_rates[INDEX(83)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(83)] - rev_rates[INDEX(83)]);

  //rxn 84
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(84)] - rev_rates[INDEX(84)]);
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(84)] - rev_rates[INDEX(84)]);
  //sp 14
  sp_rates[INDEX(13)] -= (fwd_rates[INDEX(84)] - rev_rates[INDEX(84)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(84)] - rev_rates[INDEX(84)]);

  //rxn 85
  //sp 17
  sp_rates[INDEX(16)] -= (fwd_rates[INDEX(85)] - rev_rates[INDEX(85)]);
  //sp 2
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(85)] - rev_rates[INDEX(85)]);
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(85)] - rev_rates[INDEX(85)]);
  //sp 16
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(85)] - rev_rates[INDEX(85)]);

  //rxn 86
  sp_rates[INDEX(11)] += shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 14
  sp_rates[INDEX(13)] += (fwd_rates[INDEX(86)] - rev_rates[INDEX(86)]) * pres_mod[INDEX(8)];
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] = -(fwd_rates[INDEX(86)] - rev_rates[INDEX(86)]) * pres_mod[INDEX(8)];
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(86)] - rev_rates[INDEX(86)]) * pres_mod[INDEX(8)];

  //rxn 87
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(87)] - rev_rates[INDEX(87)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(87)] - rev_rates[INDEX(87)]);
  //sp 22
  sp_rates[INDEX(21)] = (fwd_rates[INDEX(87)] - rev_rates[INDEX(87)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(87)] - rev_rates[INDEX(87)]);

  //rxn 88
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(88)] - rev_rates[INDEX(88)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(88)] - rev_rates[INDEX(88)]);
  //sp 21
  sp_rates[INDEX(20)] += (fwd_rates[INDEX(88)] - rev_rates[INDEX(88)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(88)] - rev_rates[INDEX(88)]);

  //rxn 89
  (*dy_N) += shared_temp[threadIdx.x + 1 * blockDim.x];
  //sp 22
  sp_rates[INDEX(21)] += (fwd_rates[INDEX(89)] - rev_rates[INDEX(89)]);
  //sp 6
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(89)] - rev_rates[INDEX(89)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(89)] - rev_rates[INDEX(89)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] = (fwd_rates[INDEX(89)] - rev_rates[INDEX(89)]);

  //rxn 90
  //sp 21
  sp_rates[INDEX(20)] += (fwd_rates[INDEX(90)] - rev_rates[INDEX(90)]);
  //sp 6
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(90)] - rev_rates[INDEX(90)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(90)] - rev_rates[INDEX(90)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(90)] - rev_rates[INDEX(90)]);

  //rxn 91
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(91)] - rev_rates[INDEX(91)]);
  //sp 22
  sp_rates[INDEX(21)] += (fwd_rates[INDEX(91)] - rev_rates[INDEX(91)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(91)] - rev_rates[INDEX(91)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(91)] - rev_rates[INDEX(91)]);

  //rxn 92
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(92)] - rev_rates[INDEX(92)]);
  //sp 21
  sp_rates[INDEX(20)] += (fwd_rates[INDEX(92)] - rev_rates[INDEX(92)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(92)] - rev_rates[INDEX(92)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(92)] - rev_rates[INDEX(92)]);

  //rxn 93
  sp_rates[INDEX(15)] += shared_temp[threadIdx.x];
  //sp 22
  sp_rates[INDEX(21)] += (fwd_rates[INDEX(93)] - rev_rates[INDEX(93)]);
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(93)] - rev_rates[INDEX(93)]);
  //sp 14
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(93)] - rev_rates[INDEX(93)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(93)] - rev_rates[INDEX(93)]);

  //rxn 94
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(94)] - rev_rates[INDEX(94)]);
  //sp 21
  sp_rates[INDEX(20)] += (fwd_rates[INDEX(94)] - rev_rates[INDEX(94)]);
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(94)] - rev_rates[INDEX(94)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(94)] - rev_rates[INDEX(94)]);

  //rxn 95
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(95)] - rev_rates[INDEX(95)]);
  //sp 22
  sp_rates[INDEX(21)] += (fwd_rates[INDEX(95)] - rev_rates[INDEX(95)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(95)] - rev_rates[INDEX(95)]);
  //sp 15
  sp_rates[INDEX(14)] -= (fwd_rates[INDEX(95)] - rev_rates[INDEX(95)]);

  //rxn 96
  //sp 21
  sp_rates[INDEX(20)] += (fwd_rates[INDEX(96)] - rev_rates[INDEX(96)]);
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(96)] - rev_rates[INDEX(96)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(96)] - rev_rates[INDEX(96)]);
  //sp 15
  sp_rates[INDEX(14)] -= (fwd_rates[INDEX(96)] - rev_rates[INDEX(96)]);

  //rxn 97
  //sp 10
  sp_rates[INDEX(9)] += (fwd_rates[INDEX(97)] - rev_rates[INDEX(97)]);
  //sp 22
  sp_rates[INDEX(21)] += (fwd_rates[INDEX(97)] - rev_rates[INDEX(97)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(97)] - rev_rates[INDEX(97)]);
  //sp 8
  sp_rates[INDEX(7)] -= (fwd_rates[INDEX(97)] - rev_rates[INDEX(97)]);

  //rxn 98
  sp_rates[INDEX(5)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 10
  sp_rates[INDEX(9)] += (fwd_rates[INDEX(98)] - rev_rates[INDEX(98)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] = (fwd_rates[INDEX(98)] - rev_rates[INDEX(98)]);
  //sp 23
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(98)] - rev_rates[INDEX(98)]);
  //sp 8
  sp_rates[INDEX(7)] -= (fwd_rates[INDEX(98)] - rev_rates[INDEX(98)]);

  //rxn 99
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(99)] - rev_rates[INDEX(99)]) * pres_mod[INDEX(9)];
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(99)] - rev_rates[INDEX(99)]) * pres_mod[INDEX(9)];
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(99)] - rev_rates[INDEX(99)]) * pres_mod[INDEX(9)];

  //rxn 100
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(100)] - rev_rates[INDEX(100)]) * pres_mod[INDEX(10)];
  //sp 22
  sp_rates[INDEX(21)] += (fwd_rates[INDEX(100)] - rev_rates[INDEX(100)]) * pres_mod[INDEX(10)];

  //rxn 101
  //sp 17
  sp_rates[INDEX(16)] += (fwd_rates[INDEX(101)] - rev_rates[INDEX(101)]) * pres_mod[INDEX(11)];
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(101)] - rev_rates[INDEX(101)]) * pres_mod[INDEX(11)];
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(101)] - rev_rates[INDEX(101)]) * pres_mod[INDEX(11)];

  //rxn 102
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(102)] - rev_rates[INDEX(102)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(102)] - rev_rates[INDEX(102)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(102)] - rev_rates[INDEX(102)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(102)] - rev_rates[INDEX(102)]);

  //rxn 103
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(103)] - rev_rates[INDEX(103)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(103)] - rev_rates[INDEX(103)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(103)] - rev_rates[INDEX(103)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(103)] - rev_rates[INDEX(103)]);

  //rxn 104
  sp_rates[INDEX(22)] = shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] = (fwd_rates[INDEX(104)] - rev_rates[INDEX(104)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(104)] - rev_rates[INDEX(104)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(104)] - rev_rates[INDEX(104)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(104)] - rev_rates[INDEX(104)]);

  //rxn 105
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(105)] - rev_rates[INDEX(105)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(105)] - rev_rates[INDEX(105)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(105)] - rev_rates[INDEX(105)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(105)] - rev_rates[INDEX(105)]);

  //rxn 106
  //sp 10
  sp_rates[INDEX(9)] += (fwd_rates[INDEX(106)] - rev_rates[INDEX(106)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(106)] - rev_rates[INDEX(106)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(106)] - rev_rates[INDEX(106)]);
  //sp 8
  sp_rates[INDEX(7)] -= (fwd_rates[INDEX(106)] - rev_rates[INDEX(106)]);

  //rxn 107
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(107)] - rev_rates[INDEX(107)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(107)] - rev_rates[INDEX(107)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(107)] - rev_rates[INDEX(107)]);
  //sp 8
  sp_rates[INDEX(7)] += (fwd_rates[INDEX(107)] - rev_rates[INDEX(107)]);

  //rxn 108
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(108)] - rev_rates[INDEX(108)]);
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(108)] - rev_rates[INDEX(108)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(108)] - rev_rates[INDEX(108)]);
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(108)] - rev_rates[INDEX(108)]);

  //rxn 109
  //sp 17
  sp_rates[INDEX(16)] -= (fwd_rates[INDEX(109)] - rev_rates[INDEX(109)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += 2.0 * (fwd_rates[INDEX(109)] - rev_rates[INDEX(109)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(109)] - rev_rates[INDEX(109)]);

  //rxn 110
  //sp 25
  sp_rates[INDEX(24)] += (fwd_rates[INDEX(110)] - rev_rates[INDEX(110)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(110)] - rev_rates[INDEX(110)]);
  //sp 21
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(110)] - rev_rates[INDEX(110)]);
  //sp 24
  sp_rates[INDEX(23)] -= (fwd_rates[INDEX(110)] - rev_rates[INDEX(110)]);

  //rxn 111
  sp_rates[INDEX(6)] += shared_temp[threadIdx.x + 1 * blockDim.x];
  sp_rates[INDEX(13)] += shared_temp[threadIdx.x];
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(111)] - rev_rates[INDEX(111)]) * pres_mod[INDEX(12)];
  //sp 5
  shared_temp[threadIdx.x] = (fwd_rates[INDEX(111)] - rev_rates[INDEX(111)]) * pres_mod[INDEX(12)];
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] = -(fwd_rates[INDEX(111)] - rev_rates[INDEX(111)]) * pres_mod[INDEX(12)];

  //rxn 112
  //sp 14
  sp_rates[INDEX(13)] += (fwd_rates[INDEX(112)] - rev_rates[INDEX(112)]);
  //sp 5
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(112)] - rev_rates[INDEX(112)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(112)] - rev_rates[INDEX(112)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(112)] - rev_rates[INDEX(112)]);

  //rxn 113
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(113)] - rev_rates[INDEX(113)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(113)] - rev_rates[INDEX(113)]);
  //sp 5
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(113)] - rev_rates[INDEX(113)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(113)] - rev_rates[INDEX(113)]);

  //rxn 114
  sp_rates[INDEX(20)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(114)] - rev_rates[INDEX(114)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(114)] - rev_rates[INDEX(114)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(114)] - rev_rates[INDEX(114)]);
  //sp 7
  shared_temp[threadIdx.x + 2 * blockDim.x] = (fwd_rates[INDEX(114)] - rev_rates[INDEX(114)]);

  //rxn 115
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(115)] - rev_rates[INDEX(115)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(115)] - rev_rates[INDEX(115)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(115)] - rev_rates[INDEX(115)]);
  //sp 7
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(115)] - rev_rates[INDEX(115)]);

  //rxn 116
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(116)] - rev_rates[INDEX(116)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(116)] - rev_rates[INDEX(116)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(116)] - rev_rates[INDEX(116)]);
  //sp 7
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(116)] - rev_rates[INDEX(116)]);

  //rxn 117
  sp_rates[INDEX(4)] += shared_temp[threadIdx.x];
  //sp 10
  sp_rates[INDEX(9)] += (fwd_rates[INDEX(117)] - rev_rates[INDEX(117)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(117)] - rev_rates[INDEX(117)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(117)] - rev_rates[INDEX(117)]);
  //sp 8
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(117)] - rev_rates[INDEX(117)]);

  //rxn 118
  //sp 4
  sp_rates[INDEX(3)] += (fwd_rates[INDEX(118)] - rev_rates[INDEX(118)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(118)] - rev_rates[INDEX(118)]);
  //sp 23
  sp_rates[INDEX(22)] += (fwd_rates[INDEX(118)] - rev_rates[INDEX(118)]);
  //sp 8
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(118)] - rev_rates[INDEX(118)]);

  //rxn 119
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(119)] - rev_rates[INDEX(119)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(119)] - rev_rates[INDEX(119)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(119)] - rev_rates[INDEX(119)]);
  //sp 8
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(119)] - rev_rates[INDEX(119)]);

  //rxn 120
  //sp 14
  sp_rates[INDEX(13)] -= (fwd_rates[INDEX(120)] - rev_rates[INDEX(120)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(120)] - rev_rates[INDEX(120)]);
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(120)] - rev_rates[INDEX(120)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(120)] - rev_rates[INDEX(120)]);

  //rxn 121
  //sp 25
  sp_rates[INDEX(24)] += (fwd_rates[INDEX(121)] - rev_rates[INDEX(121)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(121)] - rev_rates[INDEX(121)]);
  //sp 22
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(121)] - rev_rates[INDEX(121)]);
  //sp 24
  sp_rates[INDEX(23)] -= (fwd_rates[INDEX(121)] - rev_rates[INDEX(121)]);

  //rxn 122
  sp_rates[INDEX(6)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] = -(fwd_rates[INDEX(122)] - rev_rates[INDEX(122)]) * pres_mod[INDEX(13)];
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(122)] - rev_rates[INDEX(122)]) * pres_mod[INDEX(13)];
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(122)] - rev_rates[INDEX(122)]) * pres_mod[INDEX(13)];

  //rxn 123
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(123)] - rev_rates[INDEX(123)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(123)] - rev_rates[INDEX(123)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(123)] - rev_rates[INDEX(123)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(123)] - rev_rates[INDEX(123)]);

  //rxn 124
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(124)] - rev_rates[INDEX(124)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(124)] - rev_rates[INDEX(124)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(124)] - rev_rates[INDEX(124)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(124)] - rev_rates[INDEX(124)]);

  //rxn 125
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(125)] - rev_rates[INDEX(125)]);
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(125)] - rev_rates[INDEX(125)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(125)] - rev_rates[INDEX(125)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(125)] - rev_rates[INDEX(125)]);

  //rxn 126
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(126)] - rev_rates[INDEX(126)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(126)] - rev_rates[INDEX(126)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(126)] - rev_rates[INDEX(126)]);
  //sp 8
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(126)] - rev_rates[INDEX(126)]);

  //rxn 127
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(127)] - rev_rates[INDEX(127)]);
  //sp 18
  sp_rates[INDEX(17)] += (fwd_rates[INDEX(127)] - rev_rates[INDEX(127)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(127)] - rev_rates[INDEX(127)]);

  //rxn 128
  sp_rates[INDEX(21)] += shared_temp[threadIdx.x + 1 * blockDim.x];
  //sp 25
  sp_rates[INDEX(24)] += (fwd_rates[INDEX(128)] - rev_rates[INDEX(128)]);
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(128)] - rev_rates[INDEX(128)]);
  //sp 19
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(128)] - rev_rates[INDEX(128)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] = -(fwd_rates[INDEX(128)] - rev_rates[INDEX(128)]);

  //rxn 129
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(129)] - rev_rates[INDEX(129)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(129)] - rev_rates[INDEX(129)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(129)] - rev_rates[INDEX(129)]);
  //sp 8
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(129)] - rev_rates[INDEX(129)]);

  //rxn 130
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(130)] - rev_rates[INDEX(130)]) * pres_mod[INDEX(14)];
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(130)] - rev_rates[INDEX(130)]) * pres_mod[INDEX(14)];
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(130)] - rev_rates[INDEX(130)]) * pres_mod[INDEX(14)];

  //rxn 131
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(131)] - rev_rates[INDEX(131)]) * pres_mod[INDEX(15)];
  //sp 28
  sp_rates[INDEX(27)] = (fwd_rates[INDEX(131)] - rev_rates[INDEX(131)]) * pres_mod[INDEX(15)];
  //sp 8
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(131)] - rev_rates[INDEX(131)]) * pres_mod[INDEX(15)];

  //rxn 132
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(132)] - rev_rates[INDEX(132)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(132)] - rev_rates[INDEX(132)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(132)] - rev_rates[INDEX(132)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(132)] - rev_rates[INDEX(132)]);

  //rxn 133
  sp_rates[INDEX(18)] += shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(133)] - rev_rates[INDEX(133)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] = (fwd_rates[INDEX(133)] - rev_rates[INDEX(133)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(133)] - rev_rates[INDEX(133)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(133)] - rev_rates[INDEX(133)]);

  //rxn 134
  //sp 25
  sp_rates[INDEX(24)] += (fwd_rates[INDEX(134)] - rev_rates[INDEX(134)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(134)] - rev_rates[INDEX(134)]);
  //sp 8
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(134)] - rev_rates[INDEX(134)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(134)] - rev_rates[INDEX(134)]);

  //rxn 135
  //sp 26
  sp_rates[INDEX(25)] = (fwd_rates[INDEX(135)] - rev_rates[INDEX(135)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(135)] - rev_rates[INDEX(135)]);
  //sp 8
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(135)] - rev_rates[INDEX(135)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(135)] - rev_rates[INDEX(135)]);

  //rxn 136
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += 2.0 * (fwd_rates[INDEX(136)] - rev_rates[INDEX(136)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(136)] - rev_rates[INDEX(136)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(136)] - rev_rates[INDEX(136)]);

  //rxn 137
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(137)] - rev_rates[INDEX(137)]);
  //sp 27
  sp_rates[INDEX(26)] = (fwd_rates[INDEX(137)] - rev_rates[INDEX(137)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= 2.0 * (fwd_rates[INDEX(137)] - rev_rates[INDEX(137)]);

  //rxn 138
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(138)] - rev_rates[INDEX(138)]);
  //sp 18
  sp_rates[INDEX(17)] += (fwd_rates[INDEX(138)] - rev_rates[INDEX(138)]);
  //sp 4
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(138)] - rev_rates[INDEX(138)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(138)] - rev_rates[INDEX(138)]);

  //rxn 139
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(139)] - rev_rates[INDEX(139)]) * pres_mod[INDEX(16)];
  //sp 25
  sp_rates[INDEX(24)] += (fwd_rates[INDEX(139)] - rev_rates[INDEX(139)]) * pres_mod[INDEX(16)];
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(139)] - rev_rates[INDEX(139)]) * pres_mod[INDEX(16)];

  //rxn 140
  //sp 25
  sp_rates[INDEX(24)] += (fwd_rates[INDEX(140)] - rev_rates[INDEX(140)]);
  //sp 3
  sp_rates[INDEX(2)] -= (fwd_rates[INDEX(140)] - rev_rates[INDEX(140)]);
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(140)] - rev_rates[INDEX(140)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(140)] - rev_rates[INDEX(140)]);

  //rxn 141
  //sp 18
  sp_rates[INDEX(17)] += (fwd_rates[INDEX(141)] - rev_rates[INDEX(141)]);
  //sp 12
  sp_rates[INDEX(11)] -= (fwd_rates[INDEX(141)] - rev_rates[INDEX(141)]);
  //sp 6
  sp_rates[INDEX(5)] += (fwd_rates[INDEX(141)] - rev_rates[INDEX(141)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(141)] - rev_rates[INDEX(141)]);

  //rxn 142
  sp_rates[INDEX(7)] += shared_temp[threadIdx.x];
  //sp 25
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(142)] - rev_rates[INDEX(142)]);
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(142)] - rev_rates[INDEX(142)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(142)] - rev_rates[INDEX(142)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(142)] - rev_rates[INDEX(142)]);

  //rxn 143
  //sp 25
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(143)] - rev_rates[INDEX(143)]);
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(143)] - rev_rates[INDEX(143)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(143)] - rev_rates[INDEX(143)]);
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(143)] - rev_rates[INDEX(143)]);

  //rxn 144
  //sp 25
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(144)] - rev_rates[INDEX(144)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(144)] - rev_rates[INDEX(144)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(144)] - rev_rates[INDEX(144)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(144)] - rev_rates[INDEX(144)]);

  //rxn 145
  //sp 25
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(145)] - rev_rates[INDEX(145)]);
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(145)] - rev_rates[INDEX(145)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(145)] - rev_rates[INDEX(145)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(145)] - rev_rates[INDEX(145)]);

  //rxn 146
  //sp 25
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(146)] - rev_rates[INDEX(146)]);
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(146)] - rev_rates[INDEX(146)]);
  //sp 28
  sp_rates[INDEX(27)] += (fwd_rates[INDEX(146)] - rev_rates[INDEX(146)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(146)] - rev_rates[INDEX(146)]);

  //rxn 147
  //sp 25
  shared_temp[threadIdx.x] -= 2.0 * (fwd_rates[INDEX(147)] - rev_rates[INDEX(147)]);
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(147)] - rev_rates[INDEX(147)]);
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(147)] - rev_rates[INDEX(147)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(147)] - rev_rates[INDEX(147)]);

  //rxn 148
  sp_rates[INDEX(3)] += shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 25
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(148)] - rev_rates[INDEX(148)]) * pres_mod[INDEX(17)];
  //sp 26
  shared_temp[threadIdx.x + 3 * blockDim.x] = -(fwd_rates[INDEX(148)] - rev_rates[INDEX(148)]) * pres_mod[INDEX(17)];

  //rxn 149
  //sp 26
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(149)] - rev_rates[INDEX(149)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(149)] - rev_rates[INDEX(149)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(149)] - rev_rates[INDEX(149)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(149)] - rev_rates[INDEX(149)]);

  //rxn 150
  //sp 26
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(150)] - rev_rates[INDEX(150)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(150)] - rev_rates[INDEX(150)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(150)] - rev_rates[INDEX(150)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(150)] - rev_rates[INDEX(150)]);

  //rxn 151
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(151)] - rev_rates[INDEX(151)]);
  //sp 26
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(151)] - rev_rates[INDEX(151)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(151)] - rev_rates[INDEX(151)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(151)] - rev_rates[INDEX(151)]);

  //rxn 152
  sp_rates[INDEX(16)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 27
  shared_temp[threadIdx.x + 2 * blockDim.x] = (fwd_rates[INDEX(152)] - rev_rates[INDEX(152)]) * pres_mod[INDEX(18)];
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(152)] - rev_rates[INDEX(152)]) * pres_mod[INDEX(18)];
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(152)] - rev_rates[INDEX(152)]) * pres_mod[INDEX(18)];

  //rxn 153
  //sp 27
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(153)] - rev_rates[INDEX(153)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(153)] - rev_rates[INDEX(153)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(153)] - rev_rates[INDEX(153)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(153)] - rev_rates[INDEX(153)]);

  //rxn 154
  //sp 27
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(154)] - rev_rates[INDEX(154)]);
  //sp 4
  sp_rates[INDEX(3)] += (fwd_rates[INDEX(154)] - rev_rates[INDEX(154)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(154)] - rev_rates[INDEX(154)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(154)] - rev_rates[INDEX(154)]);

  //rxn 155
  //sp 27
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(155)] - rev_rates[INDEX(155)]);
  //sp 8
  sp_rates[INDEX(7)] += (fwd_rates[INDEX(155)] - rev_rates[INDEX(155)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(155)] - rev_rates[INDEX(155)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(155)] - rev_rates[INDEX(155)]);

  //rxn 156
  //sp 4
  sp_rates[INDEX(3)] += (fwd_rates[INDEX(156)] - rev_rates[INDEX(156)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(156)] - rev_rates[INDEX(156)]);
  //sp 8
  sp_rates[INDEX(7)] -= (fwd_rates[INDEX(156)] - rev_rates[INDEX(156)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(156)] - rev_rates[INDEX(156)]);
  //sp 27
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(156)] - rev_rates[INDEX(156)]);

  //rxn 157
  //sp 17
  sp_rates[INDEX(16)] += (fwd_rates[INDEX(157)] - rev_rates[INDEX(157)]);
  //sp 27
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(157)] - rev_rates[INDEX(157)]);
  //sp 4
  sp_rates[INDEX(3)] += (fwd_rates[INDEX(157)] - rev_rates[INDEX(157)]);

  //rxn 158
  sp_rates[INDEX(24)] += shared_temp[threadIdx.x];
  //sp 28
  shared_temp[threadIdx.x] = (fwd_rates[INDEX(158)] - rev_rates[INDEX(158)]) * pres_mod[INDEX(19)];
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(158)] - rev_rates[INDEX(158)]) * pres_mod[INDEX(19)];
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(158)] - rev_rates[INDEX(158)]) * pres_mod[INDEX(19)];

  //rxn 159
  sp_rates[INDEX(25)] += shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 27
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(159)] - rev_rates[INDEX(159)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(159)] - rev_rates[INDEX(159)]);
  //sp 28
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(159)] - rev_rates[INDEX(159)]);
  //sp 5
  shared_temp[threadIdx.x + 3 * blockDim.x] = -(fwd_rates[INDEX(159)] - rev_rates[INDEX(159)]);

  //rxn 160
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(160)] - rev_rates[INDEX(160)]);
  //sp 28
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(160)] - rev_rates[INDEX(160)]);
  //sp 5
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(160)] - rev_rates[INDEX(160)]);
  //sp 24
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(160)] - rev_rates[INDEX(160)]);

  //rxn 161
  //sp 25
  sp_rates[INDEX(24)] += (fwd_rates[INDEX(161)] - rev_rates[INDEX(161)]);
  //sp 28
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(161)] - rev_rates[INDEX(161)]);
  //sp 5
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(161)] - rev_rates[INDEX(161)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(161)] - rev_rates[INDEX(161)]);

  //rxn 162
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(162)] - rev_rates[INDEX(162)]);
  //sp 27
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(162)] - rev_rates[INDEX(162)]);
  //sp 28
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(162)] - rev_rates[INDEX(162)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(162)] - rev_rates[INDEX(162)]);

  //rxn 163
  sp_rates[INDEX(23)] += shared_temp[threadIdx.x + 1 * blockDim.x];
  sp_rates[INDEX(4)] += shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 2
  shared_temp[threadIdx.x + 3 * blockDim.x] = (fwd_rates[INDEX(163)] - rev_rates[INDEX(163)]) * pres_mod[INDEX(20)];
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] = -(fwd_rates[INDEX(163)] - rev_rates[INDEX(163)]) * pres_mod[INDEX(20)];
  //sp 6
  sp_rates[INDEX(5)] += (fwd_rates[INDEX(163)] - rev_rates[INDEX(163)]) * pres_mod[INDEX(20)];

  //rxn 164
  //sp 2
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(164)] - rev_rates[INDEX(164)]);
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(164)] - rev_rates[INDEX(164)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(164)] - rev_rates[INDEX(164)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(164)] - rev_rates[INDEX(164)]);

  //rxn 165
  //sp 2
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(165)] - rev_rates[INDEX(165)]);
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(165)] - rev_rates[INDEX(165)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(165)] - rev_rates[INDEX(165)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(165)] - rev_rates[INDEX(165)]);

  //rxn 166
  //sp 17
  sp_rates[INDEX(16)] += 2.0 * (fwd_rates[INDEX(166)] - rev_rates[INDEX(166)]);
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(166)] - rev_rates[INDEX(166)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(166)] - rev_rates[INDEX(166)]);

  //rxn 167
  //sp 2
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(167)] - rev_rates[INDEX(167)]);
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(167)] - rev_rates[INDEX(167)]);
  //sp 4
  sp_rates[INDEX(3)] += (fwd_rates[INDEX(167)] - rev_rates[INDEX(167)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(167)] - rev_rates[INDEX(167)]);

  //rxn 168
  //sp 2
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(168)] - rev_rates[INDEX(168)]);
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(168)] - rev_rates[INDEX(168)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(168)] - rev_rates[INDEX(168)]);
  //sp 8
  sp_rates[INDEX(7)] += (fwd_rates[INDEX(168)] - rev_rates[INDEX(168)]);

  //rxn 169
  sp_rates[INDEX(26)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] = -(fwd_rates[INDEX(169)] - rev_rates[INDEX(169)]);
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(169)] - rev_rates[INDEX(169)]);
  //sp 19
  sp_rates[INDEX(18)] -= (fwd_rates[INDEX(169)] - rev_rates[INDEX(169)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(169)] - rev_rates[INDEX(169)]);

  //rxn 170
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(170)] - rev_rates[INDEX(170)]);
  //sp 19
  sp_rates[INDEX(18)] -= (fwd_rates[INDEX(170)] - rev_rates[INDEX(170)]);
  //sp 12
  sp_rates[INDEX(11)] -= (fwd_rates[INDEX(170)] - rev_rates[INDEX(170)]);
  //sp 15
  sp_rates[INDEX(14)] += (fwd_rates[INDEX(170)] - rev_rates[INDEX(170)]);

  //rxn 171
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(171)] - rev_rates[INDEX(171)]);
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(171)] - rev_rates[INDEX(171)]);
  //sp 2
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(171)] - rev_rates[INDEX(171)]);
  //sp 24
  sp_rates[INDEX(23)] += (fwd_rates[INDEX(171)] - rev_rates[INDEX(171)]);

  //rxn 172
  //sp 2
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(172)] - rev_rates[INDEX(172)]);
  //sp 18
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(172)] - rev_rates[INDEX(172)]);
  //sp 12
  sp_rates[INDEX(11)] -= (fwd_rates[INDEX(172)] - rev_rates[INDEX(172)]);
  //sp 17
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(172)] - rev_rates[INDEX(172)]);

  //rxn 173
  sp_rates[INDEX(27)] += shared_temp[threadIdx.x];
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(173)] - rev_rates[INDEX(173)]);
  //sp 29
  sp_rates[INDEX(28)] = (fwd_rates[INDEX(173)] - rev_rates[INDEX(173)]);
  //sp 14
  shared_temp[threadIdx.x] = -2.0 * (fwd_rates[INDEX(173)] - rev_rates[INDEX(173)]);

  //rxn 174
  //sp 30
  sp_rates[INDEX(29)] = (fwd_rates[INDEX(174)] - rev_rates[INDEX(174)]);
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(174)] - rev_rates[INDEX(174)]);
  //sp 14
  shared_temp[threadIdx.x] -= 2.0 * (fwd_rates[INDEX(174)] - rev_rates[INDEX(174)]);

  //rxn 175
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(175)] - rev_rates[INDEX(175)]);
  //sp 29
  sp_rates[INDEX(28)] += (fwd_rates[INDEX(175)] - rev_rates[INDEX(175)]);
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(175)] - rev_rates[INDEX(175)]);
  //sp 15
  sp_rates[INDEX(14)] -= (fwd_rates[INDEX(175)] - rev_rates[INDEX(175)]);

  //rxn 176
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(176)] - rev_rates[INDEX(176)]);
  //sp 22
  sp_rates[INDEX(21)] -= (fwd_rates[INDEX(176)] - rev_rates[INDEX(176)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(176)] - rev_rates[INDEX(176)]);
  //sp 32
  sp_rates[INDEX(31)] = (fwd_rates[INDEX(176)] - rev_rates[INDEX(176)]);

  //rxn 177
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(177)] - rev_rates[INDEX(177)]);
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(177)] - rev_rates[INDEX(177)]);
  //sp 30
  sp_rates[INDEX(29)] += (fwd_rates[INDEX(177)] - rev_rates[INDEX(177)]);
  //sp 22
  sp_rates[INDEX(21)] -= (fwd_rates[INDEX(177)] - rev_rates[INDEX(177)]);

  //rxn 178
  sp_rates[INDEX(17)] += shared_temp[threadIdx.x + 1 * blockDim.x];
  //sp 14
  shared_temp[threadIdx.x] -= 2.0 * (fwd_rates[INDEX(178)] - rev_rates[INDEX(178)]) * pres_mod[INDEX(21)];
  //sp 31
  shared_temp[threadIdx.x + 1 * blockDim.x] = (fwd_rates[INDEX(178)] - rev_rates[INDEX(178)]) * pres_mod[INDEX(21)];

  //rxn 179
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(179)] - rev_rates[INDEX(179)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(179)] - rev_rates[INDEX(179)]);
  //sp 31
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(179)] - rev_rates[INDEX(179)]);
  //sp 32
  sp_rates[INDEX(31)] += (fwd_rates[INDEX(179)] - rev_rates[INDEX(179)]);

  //rxn 180
  //sp 15
  sp_rates[INDEX(14)] += (fwd_rates[INDEX(180)] - rev_rates[INDEX(180)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(180)] - rev_rates[INDEX(180)]);
  //sp 31
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(180)] - rev_rates[INDEX(180)]);
  //sp 23
  sp_rates[INDEX(22)] += (fwd_rates[INDEX(180)] - rev_rates[INDEX(180)]);

  //rxn 181
  (*dy_N) += shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(181)] - rev_rates[INDEX(181)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(181)] - rev_rates[INDEX(181)]);
  //sp 31
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(181)] - rev_rates[INDEX(181)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] = (fwd_rates[INDEX(181)] - rev_rates[INDEX(181)]);

  //rxn 182
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(182)] - rev_rates[INDEX(182)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(182)] - rev_rates[INDEX(182)]);
  //sp 31
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(182)] - rev_rates[INDEX(182)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(182)] - rev_rates[INDEX(182)]);

  //rxn 183
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(183)] - rev_rates[INDEX(183)]);
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(183)] - rev_rates[INDEX(183)]);
  //sp 31
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(183)] - rev_rates[INDEX(183)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] += (fwd_rates[INDEX(183)] - rev_rates[INDEX(183)]);

  //rxn 184
  sp_rates[INDEX(16)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] = (fwd_rates[INDEX(184)] - rev_rates[INDEX(184)]);
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(184)] - rev_rates[INDEX(184)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(184)] - rev_rates[INDEX(184)]);

  //rxn 185
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(185)] - rev_rates[INDEX(185)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(185)] - rev_rates[INDEX(185)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(185)] - rev_rates[INDEX(185)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(185)] - rev_rates[INDEX(185)]);

  //rxn 186
  sp_rates[INDEX(13)] += shared_temp[threadIdx.x];
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(186)] - rev_rates[INDEX(186)]);
  //sp 6
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(186)] - rev_rates[INDEX(186)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(186)] - rev_rates[INDEX(186)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(186)] - rev_rates[INDEX(186)]);

  //rxn 187
  //sp 14
  sp_rates[INDEX(13)] += (fwd_rates[INDEX(187)] - rev_rates[INDEX(187)]);
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(187)] - rev_rates[INDEX(187)]);
  //sp 6
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(187)] - rev_rates[INDEX(187)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(187)] - rev_rates[INDEX(187)]);

  //rxn 188
  //sp 5
  sp_rates[INDEX(4)] += fwd_rates[INDEX(188)];
  //sp 6
  shared_temp[threadIdx.x] -= fwd_rates[INDEX(188)];
  //sp 14
  sp_rates[INDEX(13)] += fwd_rates[INDEX(188)];
  //sp 17
  sp_rates[INDEX(16)] += fwd_rates[INDEX(188)];
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= fwd_rates[INDEX(188)];

  //rxn 189
  sp_rates[INDEX(30)] = shared_temp[threadIdx.x + 1 * blockDim.x];
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(189)] - rev_rates[INDEX(188)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(189)] - rev_rates[INDEX(188)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] = -(fwd_rates[INDEX(189)] - rev_rates[INDEX(188)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(189)] - rev_rates[INDEX(188)]);

  //rxn 190
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(190)] - rev_rates[INDEX(189)]);
  //sp 30
  sp_rates[INDEX(29)] += (fwd_rates[INDEX(190)] - rev_rates[INDEX(189)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(190)] - rev_rates[INDEX(189)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(190)] - rev_rates[INDEX(189)]);

  //rxn 191
  //sp 19
  sp_rates[INDEX(18)] += (fwd_rates[INDEX(191)] - rev_rates[INDEX(190)]);
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(191)] - rev_rates[INDEX(190)]);
  //sp 7
  shared_temp[threadIdx.x + 1 * blockDim.x] -= (fwd_rates[INDEX(191)] - rev_rates[INDEX(190)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(191)] - rev_rates[INDEX(190)]);

  //rxn 192
  //sp 10
  sp_rates[INDEX(9)] += (fwd_rates[INDEX(192)] - rev_rates[INDEX(191)]);
  //sp 8
  sp_rates[INDEX(7)] -= (fwd_rates[INDEX(192)] - rev_rates[INDEX(191)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(192)] - rev_rates[INDEX(191)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(192)] - rev_rates[INDEX(191)]);

  //rxn 193
  //sp 4
  sp_rates[INDEX(3)] += (fwd_rates[INDEX(193)] - rev_rates[INDEX(192)]);
  //sp 8
  sp_rates[INDEX(7)] -= (fwd_rates[INDEX(193)] - rev_rates[INDEX(192)]);
  //sp 31
  sp_rates[INDEX(30)] += (fwd_rates[INDEX(193)] - rev_rates[INDEX(192)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(193)] - rev_rates[INDEX(192)]);

  //rxn 194
  sp_rates[INDEX(5)] += shared_temp[threadIdx.x];
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(194)] - rev_rates[INDEX(193)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(194)] - rev_rates[INDEX(193)]);
  //sp 14
  shared_temp[threadIdx.x] = -(fwd_rates[INDEX(194)] - rev_rates[INDEX(193)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(194)] - rev_rates[INDEX(193)]);

  //rxn 195
  //sp 30
  sp_rates[INDEX(29)] += (fwd_rates[INDEX(195)] - rev_rates[INDEX(194)]);
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(195)] - rev_rates[INDEX(194)]);
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(195)] - rev_rates[INDEX(194)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(195)] - rev_rates[INDEX(194)]);

  //rxn 196
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(196)] - rev_rates[INDEX(195)]);
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(196)] - rev_rates[INDEX(195)]);
  //sp 15
  sp_rates[INDEX(14)] -= (fwd_rates[INDEX(196)] - rev_rates[INDEX(195)]);
  //sp 32
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(196)] - rev_rates[INDEX(195)]);

  //rxn 197
  sp_rates[INDEX(6)] += shared_temp[threadIdx.x + 1 * blockDim.x];
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(197)] - rev_rates[INDEX(196)]);
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(197)] - rev_rates[INDEX(196)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] = (fwd_rates[INDEX(197)] - rev_rates[INDEX(196)]);

  //rxn 198
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(198)] - rev_rates[INDEX(197)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(198)] - rev_rates[INDEX(197)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(198)] - rev_rates[INDEX(197)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(198)] - rev_rates[INDEX(197)]);

  //rxn 199
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(199)] - rev_rates[INDEX(198)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(199)] - rev_rates[INDEX(198)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(199)] - rev_rates[INDEX(198)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(199)] - rev_rates[INDEX(198)]);

  //rxn 200
  //sp 17
  sp_rates[INDEX(16)] += (fwd_rates[INDEX(200)] - rev_rates[INDEX(199)]);
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(200)] - rev_rates[INDEX(199)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(200)] - rev_rates[INDEX(199)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(200)] - rev_rates[INDEX(199)]);

  //rxn 201
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(201)] - rev_rates[INDEX(200)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(201)] - rev_rates[INDEX(200)]);
  //sp 7
  sp_rates[INDEX(6)] -= (fwd_rates[INDEX(201)] - rev_rates[INDEX(200)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(201)] - rev_rates[INDEX(200)]);

  //rxn 202
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(202)] - rev_rates[INDEX(201)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(202)] - rev_rates[INDEX(201)]);
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(202)] - rev_rates[INDEX(201)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(202)] - rev_rates[INDEX(201)]);

  //rxn 203
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(203)] - rev_rates[INDEX(202)]);
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(203)] - rev_rates[INDEX(202)]);
  //sp 15
  sp_rates[INDEX(14)] -= (fwd_rates[INDEX(203)] - rev_rates[INDEX(202)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(203)] - rev_rates[INDEX(202)]);

  //rxn 204
  //sp 17
  sp_rates[INDEX(16)] -= (fwd_rates[INDEX(204)] - rev_rates[INDEX(203)]);
  //sp 18
  sp_rates[INDEX(17)] += (fwd_rates[INDEX(204)] - rev_rates[INDEX(203)]);
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(204)] - rev_rates[INDEX(203)]);
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(204)] - rev_rates[INDEX(203)]);

  //rxn 205
  sp_rates[INDEX(31)] += shared_temp[threadIdx.x + 3 * blockDim.x];
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(205)] - rev_rates[INDEX(204)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] = -(fwd_rates[INDEX(205)] - rev_rates[INDEX(204)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(205)] - rev_rates[INDEX(204)]);

  //rxn 206
  //sp 5
  sp_rates[INDEX(4)] += (fwd_rates[INDEX(206)] - rev_rates[INDEX(205)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(206)] - rev_rates[INDEX(205)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(206)] - rev_rates[INDEX(205)]);

  //rxn 207
  //sp 3
  sp_rates[INDEX(2)] += (fwd_rates[INDEX(207)] - rev_rates[INDEX(206)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(207)] - rev_rates[INDEX(206)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(207)] - rev_rates[INDEX(206)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(207)] - rev_rates[INDEX(206)]);

  //rxn 208
  //sp 29
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(208)] - rev_rates[INDEX(207)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(208)] - rev_rates[INDEX(207)]);

  //rxn 209
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(209)] - rev_rates[INDEX(208)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(209)] - rev_rates[INDEX(208)]);
  //sp 7
  sp_rates[INDEX(6)] += (fwd_rates[INDEX(209)] - rev_rates[INDEX(208)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(209)] - rev_rates[INDEX(208)]);

  //rxn 210
  //sp 17
  sp_rates[INDEX(16)] += (fwd_rates[INDEX(210)] - rev_rates[INDEX(209)]);
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(210)] - rev_rates[INDEX(209)]);
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(210)] - rev_rates[INDEX(209)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(210)] - rev_rates[INDEX(209)]);

  //rxn 211
  sp_rates[INDEX(28)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 9
  sp_rates[INDEX(8)] += (fwd_rates[INDEX(211)] - rev_rates[INDEX(210)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(211)] - rev_rates[INDEX(210)]);
  //sp 7
  shared_temp[threadIdx.x + 2 * blockDim.x] = -(fwd_rates[INDEX(211)] - rev_rates[INDEX(210)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(211)] - rev_rates[INDEX(210)]);

  //rxn 212
  //sp 5
  sp_rates[INDEX(4)] += fwd_rates[INDEX(212)];
  //sp 7
  shared_temp[threadIdx.x + 2 * blockDim.x] -= fwd_rates[INDEX(212)];
  //sp 14
  shared_temp[threadIdx.x] += fwd_rates[INDEX(212)];
  //sp 17
  sp_rates[INDEX(16)] += fwd_rates[INDEX(212)];
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= fwd_rates[INDEX(212)];

  //rxn 213
  //sp 7
  shared_temp[threadIdx.x + 2 * blockDim.x] += fwd_rates[INDEX(213)];
  //sp 8
  sp_rates[INDEX(7)] -= fwd_rates[INDEX(213)];
  //sp 14
  shared_temp[threadIdx.x] += fwd_rates[INDEX(213)];
  //sp 17
  sp_rates[INDEX(16)] += fwd_rates[INDEX(213)];
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= fwd_rates[INDEX(213)];

  //rxn 214
  //sp 10
  sp_rates[INDEX(9)] += (fwd_rates[INDEX(214)] - rev_rates[INDEX(211)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(214)] - rev_rates[INDEX(211)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(214)] - rev_rates[INDEX(211)]);
  //sp 8
  sp_rates[INDEX(7)] -= (fwd_rates[INDEX(214)] - rev_rates[INDEX(211)]);

  //rxn 215
  //sp 14
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(215)] - rev_rates[INDEX(212)]);
  //sp 4
  sp_rates[INDEX(3)] -= (fwd_rates[INDEX(215)] - rev_rates[INDEX(212)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(215)] - rev_rates[INDEX(212)]);
  //sp 24
  sp_rates[INDEX(23)] += (fwd_rates[INDEX(215)] - rev_rates[INDEX(212)]);

  //rxn 216
  //sp 14
  shared_temp[threadIdx.x] -= (fwd_rates[INDEX(216)] - rev_rates[INDEX(213)]);
  //sp 13
  sp_rates[INDEX(12)] += (fwd_rates[INDEX(216)] - rev_rates[INDEX(213)]);
  //sp 30
  shared_temp[threadIdx.x + 3 * blockDim.x] -= (fwd_rates[INDEX(216)] - rev_rates[INDEX(213)]);
  //sp 16
  shared_temp[threadIdx.x + 1 * blockDim.x] += (fwd_rates[INDEX(216)] - rev_rates[INDEX(213)]);

  //rxn 217
  //sp 2
  (*dy_N) -= (fwd_rates[INDEX(217)] - rev_rates[INDEX(214)]) * pres_mod[INDEX(22)];
  //sp 12
  sp_rates[INDEX(11)] += 2.0 * (fwd_rates[INDEX(217)] - rev_rates[INDEX(214)]) * pres_mod[INDEX(22)];

  //rxn 218
  //sp 17
  sp_rates[INDEX(16)] += (fwd_rates[INDEX(218)] - rev_rates[INDEX(215)]) * pres_mod[INDEX(23)];
  //sp 12
  sp_rates[INDEX(11)] -= (fwd_rates[INDEX(218)] - rev_rates[INDEX(215)]) * pres_mod[INDEX(23)];
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(218)] - rev_rates[INDEX(215)]) * pres_mod[INDEX(23)];

  //rxn 219
  sp_rates[INDEX(6)] += shared_temp[threadIdx.x + 2 * blockDim.x];
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] = (fwd_rates[INDEX(219)] - rev_rates[INDEX(216)]) * pres_mod[INDEX(24)];
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(219)] - rev_rates[INDEX(216)]) * pres_mod[INDEX(24)];
  //sp 6
  sp_rates[INDEX(5)] -= (fwd_rates[INDEX(219)] - rev_rates[INDEX(216)]) * pres_mod[INDEX(24)];

  //rxn 220
  sp_rates[INDEX(13)] += shared_temp[threadIdx.x];
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(220)] - rev_rates[INDEX(217)]);
  //sp 7
  shared_temp[threadIdx.x] = (fwd_rates[INDEX(220)] - rev_rates[INDEX(217)]);

  //rxn 221
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(221)] - rev_rates[INDEX(218)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(221)] - rev_rates[INDEX(218)]);

  //rxn 222
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(222)] - rev_rates[INDEX(219)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(222)] - rev_rates[INDEX(219)]);

  //rxn 223
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(223)] - rev_rates[INDEX(220)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(223)] - rev_rates[INDEX(220)]);

  //rxn 224
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(224)] - rev_rates[INDEX(221)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(224)] - rev_rates[INDEX(221)]);

  //rxn 225
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(225)] - rev_rates[INDEX(222)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(225)] - rev_rates[INDEX(222)]);

  //rxn 226
  //sp 2
  (*dy_N) += (fwd_rates[INDEX(226)] - rev_rates[INDEX(223)]);
  //sp 18
  sp_rates[INDEX(17)] -= (fwd_rates[INDEX(226)] - rev_rates[INDEX(223)]);
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] += (fwd_rates[INDEX(226)] - rev_rates[INDEX(223)]);
  //sp 5
  sp_rates[INDEX(4)] -= (fwd_rates[INDEX(226)] - rev_rates[INDEX(223)]);

  //rxn 227
  //sp 11
  shared_temp[threadIdx.x + 2 * blockDim.x] -= (fwd_rates[INDEX(227)] - rev_rates[INDEX(224)]);
  //sp 7
  shared_temp[threadIdx.x] += (fwd_rates[INDEX(227)] - rev_rates[INDEX(224)]);

  //sp 0
  sp_rates[INDEX(0)] = 0.0;
  //sp 1
  sp_rates[INDEX(1)] = 0.0;
  sp_rates[INDEX(15)] += shared_temp[threadIdx.x + 1 * blockDim.x];
  sp_rates[INDEX(29)] += shared_temp[threadIdx.x + 3 * blockDim.x];
  sp_rates[INDEX(10)] = shared_temp[threadIdx.x + 2 * blockDim.x];
  sp_rates[INDEX(6)] += shared_temp[threadIdx.x];
} // end eval_spec_rates

