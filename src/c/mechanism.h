#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 0
/* Species Indexes
0  H2
1  O2
2  H
3  O
4  OH
5  HO2
6  H2O
7  H2O2
8  N
9  NH3
10  NH2
11  NH
12  NNH
13  NO
14  N2O
15  HNO
16  HON
17  H2NO
18  NO2
19  HONO
20  HONO2
21  N2H2
22  H2NN
23  N2H4
24  N2H3
25  N2
*/

//Number of species
#define NSP 26
//Number of variables. NN = NSP + 1 (temperature)
#define NN 27
//Number of forward reactions
#define FWD_RATES 171
//Number of reversible reactions
#define REV_RATES 168
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 19

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

