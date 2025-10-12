#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 2
/* Species Indexes
0  AR
1  HE
2  H
3  O
4  OH
5  H2
6  O2
7  HO2
8  H2O2
9  H2O
10  OH*
11  CO
12  CO2
13  HCO
14  CH2O
15  CH
16  C
17  CH2-3
18  CH2-1
19  C2H2
20  CH3
21  CH4
22  HOCHO
23  CH2OH
24  OCHO
25  CH3O
26  C2H4
27  CH2CO
28  C2H6
29  CH3OH
30  CH3O2
31  CH3O2H
32  C2H3
33  C2H5
34  C2H5OH
35  C2H5O
36  CH3CHO
37  CH3CO
38  CH2O2H
39  C2H
40  C2
41  C2O
42  HCCO
43  H2CC
44  C3H3
45  C4H2
46  C3H5
47  CH2CHO
48  C2H3OO
49  HCCOH
50  C4H4
51  C2H2OH
52  C2H3CHO
53  C3H4O
54  CHCHO
55  CHOCHO
56  CHOCO
57  CH3CO2
58  CH3CO3
59  CH3CO3H
60  C3H6
61  C3H8
62  I-C3H7
63  N-C3H7
64  C3H4
65  C3H4P
66  CH2CHOH
67  CH3CHOH
68  CH2CH2OH
69  C2H5O2
70  C2H4O1-2
71  C2H5O2H
72  C2H5CHO
73  C2H5CO
74  C2H4O2H
75  C2H3O1-2
76  CH3OCH3-DME
77  CH3OCH2
78  CH3OCH2O2
79  CH3OCH2O2H
80  CH3OCH2O
81  CH3OCHO
82  CH2OCHO
83  CH3OCO
84  CH2OCH2O2H
85  O2CH2OCH2O2H
86  HO2CH2OCHO
87  OCH2OCHO
88  HOCH2OCO
89  HOCH2O
90  CH*
91  N
92  NO
93  N2O
94  NO2
95  NH
96  NH2
97  NH3
98  HNO
99  HONO
100  H2NO
101  NNH
102  N2H2
103  N2H3
104  N2H4
105  H2NN
106  HNOH
107  NH2OH
108  HNO2
109  NO3
110  HONO2
111  CN
112  HCN
113  HCNO
114  NCN
115  HNC
116  HNCN
117  H2CN
118  C2N2
119  HOCN
120  HNCO
121  NCO
122  HON
123  HCNH
124  N2
*/

//Number of species
#define NSP 125
//Number of variables. NN = NSP + 1 (temperature)
#define NN 126
//Number of forward reactions
#define FWD_RATES 1099
//Number of reversible reactions
#define REV_RATES 1069
//Number of reactions with pressure modified rates
#define PRES_MOD_RATES 106

//Must be implemented by user on a per mechanism basis in mechanism.c
void set_same_initial_conditions(int, double**, double**);

#if defined (RATES_TEST) || defined (PROFILER)
    void write_jacobian_and_rates_output(int NUM);
#endif
//apply masking of ICs for cache optimized mechanisms
void apply_mask(double*);
void apply_reverse_mask(double*);
#endif

