#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 124
/* Species Indexes
0  HE
1  H
2  O
3  OH
4  H2
5  O2
6  HO2
7  H2O2
8  H2O
9  OH*
10  CO
11  CO2
12  HCO
13  CH2O
14  CH
15  C
16  CH2-3
17  CH2-1
18  C2H2
19  CH3
20  CH4
21  HOCHO
22  CH2OH
23  OCHO
24  CH3O
25  C2H4
26  CH2CO
27  C2H6
28  CH3OH
29  CH3O2
30  CH3O2H
31  C2H3
32  C2H5
33  C2H5OH
34  C2H5O
35  CH3CHO
36  CH3CO
37  CH2O2H
38  C2H
39  C2
40  C2O
41  HCCO
42  H2CC
43  C3H3
44  C4H2
45  C3H5
46  CH2CHO
47  C2H3OO
48  HCCOH
49  C4H4
50  C2H2OH
51  C2H3CHO
52  C3H4O
53  CHCHO
54  CHOCHO
55  CHOCO
56  CH3CO2
57  CH3CO3
58  CH3CO3H
59  C3H6
60  C3H8
61  I-C3H7
62  N-C3H7
63  C3H4
64  C3H4P
65  CH2CHOH
66  CH3CHOH
67  CH2CH2OH
68  C2H5O2
69  C2H4O1-2
70  C2H5O2H
71  C2H5CHO
72  C2H5CO
73  C2H4O2H
74  C2H3O1-2
75  CH3OCH3-DME
76  CH3OCH2
77  CH3OCH2O2
78  CH3OCH2O2H
79  CH3OCH2O
80  CH3OCHO
81  CH2OCHO
82  CH3OCO
83  CH2OCH2O2H
84  O2CH2OCH2O2H
85  HO2CH2OCHO
86  OCH2OCHO
87  HOCH2OCO
88  HOCH2O
89  CH*
90  N
91  NO
92  N2O
93  NO2
94  NH
95  NH2
96  NH3
97  HNO
98  HONO
99  H2NO
100  NNH
101  N2H2
102  N2H3
103  N2H4
104  H2NN
105  HNOH
106  NH2OH
107  HNO2
108  NO3
109  HONO2
110  CN
111  HCN
112  HCNO
113  NCN
114  HNC
115  HNCN
116  H2CN
117  C2N2
118  HOCN
119  HNCO
120  NCO
121  HON
122  HCNH
123  N2
124  AR
*/

//Number of species
#define NSP 125
//Number of variables. NN = NSP + 1 (temperature)
#define NN 126
//Number of forward reactions
#define FWD_RATES 1099
//Number of reversible reactions
#define REV_RATES 1071
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

