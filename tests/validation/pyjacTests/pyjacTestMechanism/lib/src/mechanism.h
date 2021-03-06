#ifndef MECHANISM_h
#define MECHANISM_h

#include <string.h>
//last_spec 47
/* Species Indexes
0  CH4
1  CH2O
2  CH3O
3  H
4  O2
5  H2
6  O
7  OH
8  H2O
9  HO2
10  H2O2
11  C
12  CH
13  CH2
14  CH2(S)
15  CH3
16  CO
17  CO2
18  HCO
19  CH2OH
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
30  N
31  NH
32  NH2
33  NH3
34  NNH
35  NO
36  NO2
37  N2O
38  HNO
39  CN
40  HCN
41  H2CN
42  HCNN
43  HCNO
44  HOCN
45  HNCO
46  NCO
47  AR
48  C3H7
49  C3H8
50  CH2CHO
51  CH3CHO
52  N2
*/

//Number of species
#define NSP 53
//Number of variables. NN = NSP + 1 (temperature)
#define NN 54
//Number of forward reactions
#define FWD_RATES 325
//Number of reversible reactions
#define REV_RATES 309
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

