#ifndef MECHANISM_cuh
#define MECHANISM_cuh

#ifdef __GNUG__
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "launch_bounds.cuh"
#include "gpu_macros.cuh"
#endif

struct mechanism_memory {
  double * y;
  double * dy;
  double * conc;
  double * fwd_rates;
  double * rev_rates;
  double * spec_rates;
  double * cp;
  double * h;
  double * dBdT;
  double * jac;
  double * var;
  double * J_nplusjplus;
  double * pres_mod;
};

//last_spec 2
/* Species Indexes
0  HE
1  AR
2  H2
3  O2
4  H
5  O
6  OH
7  HO2
8  H2O
9  H2O2
10  OH*
11  N
12  NH3
13  NH2
14  NH
15  NNH
16  NO
17  N2O
18  HNO
19  HON
20  H2NO
21  HNOH
22  NH2OH
23  NO2
24  HONO
25  HNO2
26  NO3
27  HONO2
28  N2H2
29  H2NN
30  N2H4
31  N2H3
32  N2
*/

//Number of species
#define NSP 33
//Number of variables. NN = NSP + 1 (temperature)
#define NN 34
//Number of forward reactions
#define FWD_RATES 228
//Number of reversible reactions
#define REV_RATES 225
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 25

//Must be implemented by user on a per mechanism basis in mechanism.cu
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

