#include <stdbool.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include "calcs.h"

double mayor_raiz_real_positiva(double a, double b, double c, double d, double e)
/*
 * Descripción:
 * Esta función retorna un valor de tipo double, el cual contiene la mayor raíz real positiva de una ecuación de cuarto grado con una incógnita.
 *
 * Parámetros:
 * a : coeficiente del término de cuarto grado.
 * b : coeficiente del término de tercer grado.
 * c : coeficiente del término de segundo grado.
 * d : coeficiente del término de primer grado.
 * e : término independiente.
 */
{
    double complex a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, *raiz0, *raiz;
    double A = b/a, B = c/a, C = d/a, D = e/a, raizR = 0, raizI = 1;
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

    raiz = malloc(4*sizeof*raiz);
    raiz0 = raiz;
    *raiz = a1 - b3 - .5*csqrt(b4 - b5);
    *++raiz = a1 - b3 + .5*csqrt(b4 - b5);
    *++raiz = a1 + b3 - .5*csqrt(b4 + b5);
    *++raiz = a1 + b3 + .5*csqrt(b4 + b5);

    raiz = raiz0;
    for(i = 0; i < 4; i++)
    {
        if ((creal(*raiz) > raizR) && (fabs(cimag(*raiz)) <= raizI))
        {
            raizR = creal(*raiz);
            raizI = fabs(cimag(*raiz));
        }
        raiz++;
    }

    free(raiz0);
    return raizR;
}

double *deriv(double *f, double factor, unsigned ni, double h)
/*
 * Descripción:
 * Esta función retorna un puntero de tipo double, el cual contiene la dirección de memoria para un arreglo de datos que corresponde a las muestras de la primera derivada de una función discretizada, con un orden de precisión $\mathcal{O}(h^4)$.
 *
 * Parámetros:
 * f      : puntero de memoria que contiene la dirección del arreglo de datos de la función discretizada, de la cual se desea obtener la primera derivada.
 * factor : factor que multiplica a cada una de las muestras de la función discretizada.
 * ni     : cantidad de muestras de la función discretizada.
 * h      : paso entre dos muestras consecutivas.
 */
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

double lagrange(double *y, unsigned ni, double h, double x, unsigned n)
{
    unsigned i, i0 = x/h, i1, j;
    double y_x = 0, p;

    i1 = i0 + 1;
    /* Error. Detiene la ejecución. */
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

double bin_search(double *f, unsigned ni, double h, double y_x, bool orden_ascend)
/*
 * Descripción:
 * Esta función retorna un valor de tipo double, el cual contiene la posición exacta de una muestra que se desea ubicar o la posición que resulta de la interpolación entre las dos muestras consecutivas donde está contenido el valor que se desea ubicar.
 *
 * Parámetros:
 * f            : puntero de memoria que contiene la dirección del arreglo de datos de la función discretizada, sobre la cual se realiza la búsqueda binaria.
 * ni           : cantidad de muestras de la función discretizada.
 * h : paso entre dos muestras consecutivas.
 * y_x          : muestra que se desea localizar.
 * orden_ascend : indica si las muestras de la función discretizada se encuentran ordenadas en forma ascendente (true) o en forma descendente (false).
 */
{
    unsigned x_a = 0, x_m, x_b = ni + 1;

    if (orden_ascend)
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

short bulk(short s, double alpha, double beta, double gamma, double tildegamma, double m, double *f1oo, double *f2oo, double *lambda_eff)
/*
 * Descripción:
 * Esta función retorna un valor de tipo short. Si el resultado es $0$, entonces el cálculo se realizó correctamente, de lo contrario, si es $1$, existe una inconsistencia en los datos de entrada.
 *
 * Parámetros:
 * s           : valor del parámetro $s$.
 * alpha       : valor del parámetro $\alpha$.
 * beta        : valor del parámetro $\beta$.
 * gamma       : valor del parámetro $\gamma$.
 * tildegamma  : valor del parámetro $\tilde\gamma$.
 * m           : valor del parámetro $m$.
 * f1oo        : puntero que almacena el valor de la condición de frontera en el bulk para la primera banda.
 * f2oo        : puntero que almacena el valor de la condición de frontera en el bulk para la segunda banda.
 * lambda_eff  : profundidad de penetración efectiva.
 */
{
    double etaR = mayor_raiz_real_positiva(beta*gamma, -s*beta - alpha*tildegamma, 0, alpha + s*tildegamma, -gamma),
           aux = (gamma*etaR - s)/(1 - tildegamma*etaR*etaR);

    if ((aux <= 0) || (!isfinite(aux)))
        return 1;

    *f1oo = sqrt(aux);
    *f2oo = etaR**f1oo;
    *lambda_eff = 1/(*f1oo*sqrt(1 + m*etaR*etaR));
    return 0;
}

void ansatz(double *cf, unsigned ni, double *f1, double *f2, double *a)
/*
 * Descripción:
 * Esta función tiene como propósito definir una solución estimada inicial para las funciones incógnitas $f_1$, $f_2$ y $a$.
 *
 * Parámetros:
 * cf : puntero donde se almacenan las condiciones de frontera en el bulk de las funciones discretizadas.
 * ni : cantidad de muestras de las funciones discretizadas.
 * f1 : puntero que almacena la estimación inicial de la solución numérica para el parámetro de orden $f_1(r)$.
 * f2 : puntero que almacena la estimación inicial de la solución numérica para el parámetro de orden $f_2(r)$.
 * a  : puntero que almacena la estimación inicial de la solución numérica para la función $a(r)$.
 */
{
    unsigned i;
    double *f10, *f20, *a0, f1oo, f2oo, aoo;

    f10 = f1;
    f20 = f2;
    a0 = a;

    /* Condiciones de frontera. */
    *f1++ = 0;
    *f2++ = 0;
    *a++ = 0;
    f1oo = *cf;
    f2oo = cf[1];
    aoo = cf[2];

    /* Ansatz's. */
    i = 0;
    while (i++ < ni + 1)
    {
        *f1++ = f1oo;
        *f2++ = f2oo;
        *a++ = aoo;
    }

    f1 = f10;
    f2 = f20;
    a = a0;
}

double *campo(double *a, double factor, unsigned ni, double h)
/*
 * Descripción:
 * Esta función retorna un puntero de tipo double, el cual contiene la dirección de memoria para un arreglo de datos que corresponde a la inducción de campo magnético, $B$.
 *
 * Parámetros:
 * a      : puntero que contiene la dirección de memoria donde se almacena la solución numérica de la función $a(r)$.
 * factor : factor $n/\kappa_1$ que multiplica a cada una de las muestras de la función $a(r)$.
 * ni     : cantidad de muestras de la función discretizada.
 * h      : paso entre dos muestras consecutivas.
 */
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

    /* Campo magnético aplicado. Otra forma: lagrange(B0, ni, h, 0, 2);*/
    *B0 = 3*B0[1] - 3*B0[2] + B0[3];

    return B0;
}

short solve(double* params, double* cf, unsigned ni, double h, double max_err, unsigned max_iter, double *f1, double *f2, double *a)
/*
 * Descripción:
 * Esta función retorna un valor de tipo short. Si el resultado es $0$, entonces el cálculo se realizó correctamente, de lo contrario, si el valor devuelto es $1$ existe una inconsistencia en los argumentos de entrada. El propósito de esta función es obtener las soluciones de GL, representadas en las variables f1, f2 y a.
 *
 * Parámetros:
 * params   : puntero de memoria que contiene la dirección de un arreglo de datos con los valores de los parámetros $s$, $\alpha$, $\beta$, $\gamma$, $\gt$, $m$, $\kappa_1$ y $n$.
 * cf       : puntero donde se almacenan las condiciones de frontera en el bulk de cada función incógnita.
 * ni       : cantidad de muestras de las funciones discretizadas.
 * h        : paso entre dos muestras consecutivas.
 * max_err  : error máximo de la solución numérica.
 * max_iter : cantidad máxima de iteraciones previstas para el método de Newton-Raphson.
 * f1       : puntero donde se almacena la solución numérica para el parámetro de orden $f_1(r)$.
 * f2       : puntero donde se almacena la solución numérica para el parámetro de orden $f_2(r)$.
 * a        : puntero donde se almacena la solución numérica para la función $a(r)$.
 */
{
    double *f11, *f21, *a1, *g, *g0, *y, *y0, *J, *J0,
           s, alpha, beta, gamma, tildegamma, m, kappa1, n,
           c1, c2, c3, c4,
           aux, aux1, aux2, aux3, aux4, f1_2, f2_2,
           pivote, sp, termino;
    unsigned short k1, k2, k_max;
    unsigned i, j, k, ni3, f, c;
    size_t tam3ni;

    /* Recuperación de parámetros. */
    s = *params++;
    alpha = *params++;
    beta = *params++;
    gamma = *params++;
    tildegamma = *params++;
    m = *params++;
    kappa1 = *params++;
    n = *params;

    /* Reducción de las operaciones. */
    c1 = h*h;
    c2 = kappa1*kappa1*c1;
    c3 = c2/m;
    c4 = -2*n*n;

    /* Reserva de memoria. */
    ni3 = 3*ni;
    tam3ni = ni3*sizeof(double);
    y = malloc(tam3ni);
    g = malloc(tam3ni);
    J = malloc(ni3*ni3*sizeof*J);
    g0 = g;
    y0 = y;
    J0 = J;

    /* Ansatz's. */
    ansatz(cf, ni, f1, f2, a);
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

        /* Eliminación gaussiana (algoritmo optimizado para este problema en particular). */
        k_max = 3;
        aux1 = ni3 - 1;
        aux2 = aux1 - 3;
        for (i = 0; i < aux1; i++)
            /* Verifica que la matriz no sea singular. */
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

        /* Sustitución regresiva. */
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

        /* Verifica si el error de la solución es menor al establecido por el usuario. */
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
