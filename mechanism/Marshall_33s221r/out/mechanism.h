#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 38
/* Species Indexes
0  NH3
1  NO2
2  NO
3  N2O
4  H2
5  O2
6  O3
7  H
8  O
9  OH
10  HO2
11  H2O
12  H2O2
13  CO
14  CO2
15  HOCO
16  CH2O
17  HCO
18  NH2
19  NH
20  N
21  NNH
22  N2H4
23  N2H3
24  tHNNH
25  cHNNH
26  H2NN
27  NH2OH
28  H2NO
29  HNOH
30  HNO
31  HON
32  HONO
33  HNO2
34  NO3
35  HONO2
36  H2NCO
37  HNCO
38  HE
39  N2
40  AR
*/

//Number of species
#define NSP 41
//Number of variables. NN = NSP + 1 (temperature)
#define NN 42
//Number of forward reactions
#define FWD_RATES 271
//Number of reversible reactions
#define REV_RATES 268
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 37

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

