#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 4
/* Species Indexes
0  H2O
1  NH3
2  NO
3  N2O
4  O2
5  H2
6  NO2
7  AR
8  HE
9  H
10  O
11  OH
12  HO2
13  H2O2
14  OHV
15  NH2
16  NH
17  N
18  N2H4
19  N2H3
20  N2H2
21  H2NN
22  NNH
23  NO3
24  HNO
25  HON
26  H2NO
27  HNOH
28  NH2OH
29  HONO
30  HNO2
31  HONO2
32  t-ONNH
33  c-ONNH
34  ONHN
35  H2NNO2
36  t-HNN(O)OH
37  c-HNN(O)OH
38  CO
39  CO2
40  CH4
41  C2H6
42  N2
*/

//Number of species
#define NSP 43
//Number of variables. NN = NSP + 1 (temperature)
#define NN 44
//Number of forward reactions
#define FWD_RATES 312
//Number of reversible reactions
#define REV_RATES 306
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

