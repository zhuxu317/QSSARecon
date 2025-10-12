#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 15
/* Species Indexes
0  H
1  H2
2  H2O
3  H2O2
4  HO2
5  N
6  N2
7  N2O
8  NH
9  NH2
10  NNH
11  NO
12  O
13  O2
14  OH
15  AR
*/

//Number of species
#define NSP 16
//Number of variables. NN = NSP + 1 (temperature)
#define NN 17
//Number of forward reactions
#define FWD_RATES 47
//Number of reversible reactions
#define REV_RATES 47
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 7

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

