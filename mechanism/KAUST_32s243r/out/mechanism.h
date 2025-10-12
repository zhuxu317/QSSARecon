#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 1
/* Species Indexes
0  AR
1  H
2  H2
3  HE
4  O
5  OH
6  H2O
7  HO2
8  H2O2
9  O2
10  CO
11  CO2
12  HCO
13  HOCO
14  N
15  NH
16  NH2
17  NH3
18  NNH
19  N2H2
20  N2H3
21  N2H4
22  H2NN
23  NO
24  N2O
25  NO2
26  HNO
27  HON
28  HONO
29  H2NO
30  HNOH
31  NH2OH
32  HNO2
33  HONO2
34  NO3
35  N2
*/

//Number of species
#define NSP 36
//Number of variables. NN = NSP + 1 (temperature)
#define NN 37
//Number of forward reactions
#define FWD_RATES 274
//Number of reversible reactions
#define REV_RATES 269
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

