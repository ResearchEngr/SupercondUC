#ifndef CALCS_H_INCLUDED
#define CALCS_H_INCLUDED

/* Prototipos de funciones. */
void ansatz(double*, unsigned, double*, double*, double*);
double bin_search(double*, unsigned, double, double, bool);
short bulk(short, double, double, double, double, double, double*, double*, double*);
double *campo(double*, double, unsigned, double);
double *deriv(double*, double, unsigned, double);
double lagrange(double*, unsigned, double, double, unsigned);
double mayor_raiz_real_positiva(double, double, double, double, double);
short solve(double*, double*, unsigned, double, double, unsigned, double*, double*, double*);

#endif
