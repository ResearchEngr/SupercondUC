#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include "calcs.h"

/*
 * Returns the greatest positive real root for a quartic equation in one unknown
 */
double greatest_positive_real_root(double a, double b, double c, double d, double e)
{
    double complex a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, *root0, *root;
    double A = b/a, B = c/a, C = d/a, D = e/a, real_root = 0, imag_root = 1;
    short i;

    a1 = -A/4;
    a2 = (3*A*A - 8*B)/12;
    a3 = B*B - 3*A*C + 12*D;
    a4 = 2*B*B*B - 9*A*B*C + 27*C*C + 27*A*A*D - 72*B*D;
    a5 = -(A*A*A)/8 + A*B/2 - C;

    b1 = cpow((a4 + csqrt(a4*a4 - 4*a3*a3*a3))/54, 1.0/3);
    b2 = a3/(9*b1);
    b3 = .5*csqrt(a2 + b1 + b2);
    b4 = 2*a2 - b1 - b2;
    b5 = a5/b3;

    root = malloc(4*sizeof*root);
    root0 = root;
    *root = a1 - b3 - .5*csqrt(b4 - b5);
    *++root = a1 - b3 + .5*csqrt(b4 - b5);
    *++root = a1 + b3 - .5*csqrt(b4 + b5);
    *++root = a1 + b3 + .5*csqrt(b4 + b5);

    root = root0;
    for(i = 0; i < 4; i++)
    {
        if ((creal(*root) > real_root) && (fabs(cimag(*root)) <= imag_root))
        {
            real_root = creal(*root);
            imag_root = fabs(cimag(*root));
        }
        root++;
    }

    free(root0);
    return real_root;
}

/*
 * Returns a pointer to an array containing the numerical differentiation of a discrete function, within a order of approximation $\mathcal{O}(h^4)$
 */
double *num_diff(double *f, double factor, unsigned ni, double h)
{
    double *d, *d0;
    unsigned i, n = ni + 2;

    d = malloc(n*sizeof(double));
    d0 = d;
    factor /= 12*h;

    for (i = 0; i < n; i++)
        if (i < 2) *d++ = factor*(-25*f[i] + 48*f[i + 1] - 36*f[i + 2] + 16*f[i + 3] - 3*f[i + 4]);
        else if (i < ni) *d++ = factor*(-f[i + 2] + 8*f[i + 1] - 8*f[i - 1] + f[i - 2]);
        else *d++ = factor*(25*f[i] - 48*f[i - 1] + 36*f[i - 2] - 16*f[i - 3] + 3*f[i - 4]);

    return d0;
}

/*
 * Gets an interpolation by Lagrange polynomial
 */
double lagrange(double *y, unsigned ni, double h, double x, unsigned n)
{
    unsigned i, i0 = x/h, i1, j;
    double y_x = 0, p;

    i1 = i0 + 1;
    /* Error. Stops execution */
    if ((n += i1 + 1) > ni) exit(1);

    y += i1;
    for (i = i1; i < n; i++)
    {
        p = *y++;
        for (j = i1; j < n; j++) if (i != j) p *= ((double)i0 - j)/(i - j);
        y_x += p;
    }
    return y_x;
}

/*
 * Binary search modified algorithm: returns the position of the exact element searched
 * or a linear interpolation between the two positions where it is contained
 */
double bin_search(double *f, unsigned ni, double h, double y_x, bool sort_order)
{
    unsigned x_a = 0, x_m, x_b = ni + 1;

    if (sort_order)
        while (x_a + 1 < x_b)
            if (y_x < f[x_m = (x_a + x_b)/2])
                x_b = x_m;
            else if (y_x > f[x_m])
                x_a = x_m;
            else
                return h*x_m;
    else
        while (x_a + 1 < x_b)
            if (y_x > f[x_m = (x_a + x_b)/2])
                x_b = x_m;
            else if (y_x < f[x_m])
                x_a = x_m;
            else
                return h*x_m;

    return h*(x_a + (y_x - f[x_a])/(f[x_b] - f[x_a]));
}

/*
 * Set the boundary conditions for the unknowns in the bulk. Returns $0$ if the calculations were successful,
 * otherwise $1$, indicating an inconsistency in the entry data
 */
short bulk(short s, double alpha, double beta, double gamma, double tildegamma, double m, double *f1_infty, double *f2_infty, double *lambda_eff)
{
    double etaR = greatest_positive_real_root(beta*gamma, -s*beta - alpha*tildegamma, 0, alpha + s*tildegamma, -gamma),
           aux = (gamma*etaR - s)/(1 - tildegamma*etaR*etaR);

    if ((aux <= 0) || (!isfinite(aux)))
        return 1;

    *f1_infty = sqrt(aux);
    *f2_infty = etaR**f1_infty;
    *lambda_eff = 1/(*f1_infty*sqrt(1 + m*etaR*etaR));
    return 0;
}

/*
 * Set the initial guess for the unknowns
 */
void ansatz(double *bc, unsigned ni, double *f1, double *f2, double *a)
{
    unsigned i;
    double *f10, *f20, *a0, f1_infty, f2_infty, a_infty;

    f10 = f1;
    f20 = f2;
    a0 = a;

    /* Boundary conditions */
    *f1++ = 0;
    *f2++ = 0;
    *a++ = 0;
    f1_infty = *bc;
    f2_infty = bc[1];
    a_infty = bc[2];

    /* Ansatz's */
    i = 0;
    while (i++ < ni + 1)
    {
        *f1++ = f1_infty;
        *f2++ = f2_infty;
        *a++ = a_infty;
    }

    f1 = f10;
    f2 = f20;
    a = a0;
}

/*
 * Calculation of the magnetic field $B$
 */
double *magnetic_field(double *a, double factor, unsigned ni, double h)
{
    double *B, *B0;
    unsigned i, n = ni + 2;

    B = malloc(n*sizeof(double));
    B0 = B++;
    factor /= 12*h*h;

    *B++ = factor*(-25*a[1] + 48*a[2] - 36*a[3] + 16*a[4] - 3*a[5]);

    for (i = 2; i < n; i++)
        if (i < ni) *B++ = factor*(-a[i + 2] + 8*a[i + 1] - 8*a[i - 1] + a[i - 2])/i;
        else *B++ = factor*(25*a[i] - 48*a[i - 1] + 36*a[i - 2] - 16*a[i - 3] + 3*a[i - 4])/i;

    /*
     * Calculating magnetic field in the core
     * Another approach: lagrange(B0, ni, h, 0, 2);
     */
    *B0 = 3*B0[1] - 3*B0[2] + B0[3];

    return B0;
}

/*
 * Get the GL solutions and set the calculated values to pointers f1, f2 and a. Returns $0$ if the calculations were successful,
 * otherwise $1$, indicating an invalid solution (inconsistency on parameter values)
 */
short solve(double* params, double* bc, unsigned ni, double h, double max_err, unsigned max_iter, double *f1, double *f2, double *a)
{
    double *f11, *f21, *a1, *g, *g0, *y, *y0, *J, *J0,
           s, alpha, beta, gamma, tildegamma, m, kappa1, n,
           c1, c2, c3, c4,
           aux, aux1, aux2, aux3, aux4, f1_2, f2_2,
           pivote, sp, termino;
    unsigned short k1, k2, k_max;
    unsigned i, j, k, ni3, f, c;
    size_t tam3ni;

    /* Recovering the parameter values */
    s = *params++;
    alpha = *params++;
    beta = *params++;
    gamma = *params++;
    tildegamma = *params++;
    m = *params++;
    kappa1 = *params++;
    n = *params;

    /* Optimizing calculations */
    c1 = h*h;
    c2 = kappa1*kappa1*c1;
    c3 = c2/m;
    c4 = -2*n*n;

    /* Allocating required memory */
    ni3 = 3*ni;
    tam3ni = ni3*sizeof(double);
    y = malloc(tam3ni);
    g = malloc(tam3ni);
    J = malloc(ni3*ni3*sizeof*J);
    g0 = g;
    y0 = y;
    J0 = J;

    /* Ansatz's */
    ansatz(bc, ni, f1, f2, a);
    f11 = ++f1;
    f21 = ++f2;
    a1 = ++a;

    for (k = 0; k < max_iter; k++)
    {
        for (i = 1; i <= ni; i++)
        {
            aux = 1.0/i;
            aux1 = pow(n*(*a - 1)*aux, 2);
            aux2 = aux*aux;
            aux *= .5;
            aux3 = 1 + aux;
            aux4 = 1 - aux;
            f1_2 = *f1**f1;
            f2_2 = *f2**f2;

            *g++ = aux3**(f1 + 1) + aux4**(f1 - 1) - 2**f1 - aux1**f1 - c2*(s**f1 + *f1*f1_2 - gamma**f2 - tildegamma**f1*f2_2);
            if (i > 1) *(J - 3) = aux4;
            *J = -2 - aux1 - c2*(s + 3*f1_2 - tildegamma*f2_2);
            *(J + 1) = c2*(gamma + 2*tildegamma**f1**f2);
            *(J + 2) = c4*(*a - 1)**f1*aux2;
            if (i < ni) *(J + 3) = aux3;

            J += ni3 + 1;
            *g++ = aux3**(f2 + 1) + aux4**(f2 - 1) - 2**f2 - aux1**f2 - c3*(alpha**f2 + beta**f2*f2_2 - gamma**f1 - tildegamma*f1_2**f2);
            if (i > 1) *(J - 3) = aux4;
            *(J - 1) = c3*(gamma + 2*tildegamma**f1**f2);
            *J = -2 - aux1 - c3*(alpha + 3*beta*f2_2 - tildegamma*f1_2);
            *(J + 1) = c4*(*a - 1)**f2*aux2;
            if (i < ni) *(J + 3) = aux3;

            J += ni3 + 1;
            *g++ = aux4**(a + 1) + aux3**(a - 1) - 2**a - c1*(f1_2 + m*f2_2)*(*a - 1);
            if (i > 1) *(J - 3) = aux3;
            *(J - 2) = -2*c1*(*a - 1)**f1;
            *(J - 1) = -2*m*c1*(*a - 1)**f2;
            *J = -2 - c1*(f1_2 + m*f2_2);
            if (i < ni) *(J + 3) = aux4;

            J += ni3 + 1;
            f1++;
            f2++;
            a++;
        }

        f1 = f11;
        f2 = f21;
        a = a1;
        J = J0;
        g = g0;

        /* Gaussian elimination (optimized for this particular problem) */
        k_max = 3;
        aux1 = ni3 - 1;
        aux2 = aux1 - 3;
        for (i = 0; i < aux1; i++)
            /* Verify that the Jacobian matrix is not singular */
            if (fabs(pivote = J[i]) > 0)
            {
                if (i > aux2) k_max--;
                f = 0;

                for (k1 = 0; k1 < k_max; k1++)
                    if (fabs(J[++f*ni3 + i]) > 0)
                    {
                        termino = J[f*ni3 + i]/pivote;
                        J[f*ni3 + i] = 0;
                        c = i;
                        for (k2 = 0; k2 < k_max; k2++)
                            if (fabs(J[++c]) > 0)
                                J[f*ni3 + c] -= termino*J[c];
                        g[f] -= termino**g;
                    }
                J += ni3;
                g++;
            }
            else return 1;

        /* Back-Substitution algorithm */
        i = ni3;
        while (i--)
        {
            sp = 0;
            for (j = i + 1; j < ni3; j++)
            {
                sp += J[j]*y[j];
                J[j] = 0;
            }
            y[i] = (-*g-- - sp)/J[i];
            J -= ni3;
        }
        J = J0;
        g = g0;

        for (i = 0; i < ni; i++)
        {
            *f1++ += *y++;
            *f2++ += *y++;
            *a++ += *y++;
        }

        f1 = f11;
        f2 = f21;
        a = a1;
        y = y0;

        /* Check whether the solution error is less than the error set by user */
        sp = 0;
        for (i = 0; i < ni3; i++)
        {
            sp += *y**y;
            y++;
        }
        if (sqrt(sp) < max_err) break;
        y = y0;
    }
    free(J0);
    free(y0);
    free(g0);
    f1--;
    f2--;
    a--;
    return 0;
}
