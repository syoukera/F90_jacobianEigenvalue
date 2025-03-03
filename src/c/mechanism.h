#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
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

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

