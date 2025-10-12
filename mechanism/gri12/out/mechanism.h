#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 31
/* Species Indexes
0  H2
1  H
2  O
3  O2
4  OH
5  H2O
6  HO2
7  H2O2
8  C
9  CH
10  CH2
11  CH2(S)
12  CH3
13  CH4
14  CO
15  CO2
16  HCO
17  CH2O
18  CH2OH
19  CH3O
20  CH3OH
21  C2H
22  C2H2
23  C2H3
24  C2H4
25  C2H5
26  C2H6
27  HCCO
28  CH2CO
29  HCCOH
30  N2
31  AR
*/

//Number of species
#define NSP 32
//Number of variables. NN = NSP + 1 (temperature)
#define NN 33
//Number of forward reactions
#define FWD_RATES 177
//Number of reversible reactions
#define REV_RATES 177
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 35

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

