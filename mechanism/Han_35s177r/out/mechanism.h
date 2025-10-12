#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 34
/* Species Indexes
0  H
1  H2
2  O
3  O2
4  OH
5  H2O
6  HO2
7  H2O2
8  CO
9  CO2
10  HE
11  HCO
12  OHEX
13  N
14  NO
15  NO2
16  N2O
17  NH3
18  HNO
19  H2NO
20  NNH
21  NH2
22  NH
23  N2H2
24  N2H3
25  N2H4
26  H2NN
27  HONO
28  HNOH
29  HNO2
30  HONO2
31  NO3
32  HON
33  N2
34  AR
*/

//Number of species
#define NSP 35
//Number of variables. NN = NSP + 1 (temperature)
#define NN 36
//Number of forward reactions
#define FWD_RATES 177
//Number of reversible reactions
#define REV_RATES 174
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 24

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

