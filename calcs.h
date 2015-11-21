#ifndef CALCS_H_INCLUDED
#define CALCS_H_INCLUDED

/* Function prototypes */
void ansatz(double*, unsigned, double*, double*, double*);
double bin_search(double*, unsigned, double, double, bool);
short bulk(short, double, double, double, double, double, double*, double*, double*);
double *magnetic_field(double*, double, unsigned, double);
double *num_diff(double*, double, unsigned, double);
double lagrange(double*, unsigned, double, double, unsigned);
double greatest_positive_real_root(double, double, double, double, double);
short solve(double*, double*, unsigned, double, double, unsigned, double*, double*, double*);

#endif
