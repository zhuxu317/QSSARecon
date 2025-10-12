#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 28
/* Species Indexes
0  N2
1  H
2  H2
3  O
4  OH
5  O2
6  H2O2
7  H2O
8  HO2
9  CO
10  CH3
11  CH2O
12  CO2
13  CH4
14  C2H2
15  C2H4
16  CH2CO
17  C2H6
18  C
19  CH
20  HCO
21  TXCH2
22  SXCH2
23  C2H3
24  C2H5
25  HCCO
26  CH3CHO
27  CH2CHO
28  C2H5O
*/

//Number of species
#define NSP 29
//Number of variables. NN = NSP + 1 (temperature)
#define NN 30
//Number of forward reactions
#define FWD_RATES 172
//Number of reversible reactions
#define REV_RATES 158
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

