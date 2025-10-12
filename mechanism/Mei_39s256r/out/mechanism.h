#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 37
/* Species Indexes
0  H
1  H2
2  HE
3  O
4  OH
5  H2O
6  HO2
7  H2O2
8  O2
9  OH*
10  N
11  NH
12  NH2
13  NH3
14  NNH
15  N2H2
16  N2H3
17  N2H4
18  H2NN
19  NO
20  N2O
21  NO2
22  HNO
23  HON
24  HONO
25  H2NO
26  HNOH
27  NH2OH
28  HNO2
29  HONO2
30  NO3
31  HNO3
32  CO
33  CO2
34  CH4
35  C2H6
36  N2
37  AR
*/

//Number of species
#define NSP 38
//Number of variables. NN = NSP + 1 (temperature)
#define NN 39
//Number of forward reactions
#define FWD_RATES 265
//Number of reversible reactions
#define REV_RATES 262
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 41

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

