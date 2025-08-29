#include <math.h>
#include "header.cuh"

__device__ void eval_jacob_4 (const double pres, const double * __restrict__ conc, const double * __restrict__ fwd_rates, const double * __restrict__ rev_rates, const double * __restrict__ pres_mod, const double * __restrict__ spec_rates, const double m, const double mw_avg, const double rho, const double * __restrict__ dBdT, const double T, double * __restric