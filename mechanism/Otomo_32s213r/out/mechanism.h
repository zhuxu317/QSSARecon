#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 30
/* Species Indexes
0  NO
1  NH3
2  H2
3  O2
4  H
5  O
6  OH
7  HO2
8  N2
9  H2O
10  H2O2
11  NH2
12  NH
13  N
14  NNH
15  NH2OH
16  H2NO
17  HNOH
18  HNO
19  HON
20  NO2
21  HONO
22  HNO2
23  NO3
24  HONO2
25  N2O
26  N2H4
27  N2H3
28  N2H2
29  H2NN
30  AR
*/

//Number of species
#define NSP 31
//Number of variables. NN = NSP + 1 (temperature)
#define NN 32
//Number of forward reactions
#define FWD_RATES 213
//Number of reversible reactions
#define REV_RATES 210
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 26

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

