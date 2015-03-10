#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <windows.h>
#include <Lmcons.h>
#endif
#include "calcs.h"
#include "filerw.h"
#include "pipes.h"

int main(int argc, char **argv)
{
    double temp_c, temp_c1, temp_c2, temp, temp_ajuste, ctetemp_c1, ctetemp_c2, cte_temp, gamma_a, r_infty, h, h_bar, *f1, *f2, *a, *B, *df1, *df2, f1oo, f2oo, x0, x, df10, df20, coef[9], cf[3], alpha, beta, m, gamma, tildegamma, kappa_1, w, k, q0, lambda_eff, lambda, lh1, lh2, la, lh_mayor, lh_mayor_ant = 0, factor_lh, factor_r, approx_err = .1, *lambda_eff_var, *lambda_var, *lh1_var, *lh2_var, *la_var, *f1oo_var, *f2oo_var, *df10_var, *df20_var, *B0_var, *gamma_var, *otra_var, max_err, aux, temp0, tildegamma0;
    int i, j, n_bar, n_val, ni;
    short s, o_tit, o_calc, o_graf, o_res, escala, ajuste, n_res, max_iter = 100;
    unsigned short n, decimales, *cod_err, n_err, prctje, prctje_ant = 0;
    char *material, *nom_arch, **nom_arch_aux, *sufijo, **id_arch, *dir, *dir_prefix, form_arch[4], char_coment, ext[3], *nom_soft_graf, *usuario, *var_ind, *valor, **formato, *ident, *ident_dup, *ident_dup0, **msj_err;
    size_t tam_sol;
    bool *es_invalido;
    struct stat st = {0};
    superconductor* sc;
    info_grafico *info_gr;
    FILE **archivo;

    switch(argc)
    {
        {
        /* Inicio de bloque. Interpretación de comandos. */
        case 1:

            error(5, "debe especificar los argumentos de entrada");
            break;

        default:

#ifdef __linux__
usuario = getenv("USER");
#elif _WIN32
usuario = malloc(UNLEN + 1);
DWORD usuario_len = UNLEN + 1;
GetUserName(usuario, &usuario_len);
#endif

            /* Valores predeterminados de los argumentos de entrada. */
            temp = 0;
            tildegamma = .0;
            kappa_1 = 1/sqrt(2);
            n = 1;
            factor_r = 5;
            factor_lh = .95;
            ni = 750;
            o_calc = 0;
            o_tit = 1;
            n_bar = 750;
            decimales = 10;
            strcpy(ext, "dat");
            o_graf = 3;
            o_res = 0;
            escala = 0;
            ajuste = 0;
            strcpy(form_arch, "png");

            cod_err = malloc((argc - 1)*sizeof*cod_err);
            msj_err = malloc((argc - 1)*sizeof msj_err);
            ident = malloc(argc - 1);

            /* Verificación de la sintaxis y validación de los argumentos introducidos por el usuario.
             * Inicializa la cadena de caracteres "ident" en un valor nulo. De lo contrario, la cadena
             * contendrá basura y generará errores.
             */
            strncpy(ident, "", argc);
            ident_dup = malloc(argc - 1);
            ident_dup0 = ident_dup;
            n_err = 0;

            for(i = 1; i < argc; i++)
            {
                /*
                 * Se verifica que el primer caracter sea un identificador válido.
                 * m : material.
                 * T : temperatura.
                 * g : valor de $\tilde\gamma$.
                 * k : valor de $\kappa_1$.
                 * n : vorticidad.
                 * r : factor de $r_\infty$.
                 * l : factor de $L_h$.
                 * N : cantidad de puntos (diferencias finitas).
                 * o : opción de cálculo (ninguno -> 0, temperatura -> 1, $\gamma_a$ -> 2, $\tilde\gamma$ -> 3).
                 * P : cantidad de puntos (barridos de T, $\gamma$ o $\gamma_a$).
                 * d : cantidad de decimales.
                 * i : formato de los archivos imagen de salida.
                 * f : formato de archivo del software de representación gráfica.
                 * G : especifica si se deja abierta la aplicación de representación gráfica.
                 * R : opción de resolución (800x600 -> 0, 1024x768 -> 1, 1920x1080 -> 2, 2560x2048 -> 3).
                 * t : especifica si se coloca título o no al gráfico.
                 * e : especifica si $T$ se escala a $T_c$.
                 * a : especifica si se ajusta la temperatura inicial en un barrido de $T$.
                 */

                if (strchr("mTgknrlNoPdifGRtea", *argv[i]) == NULL)
                {
                    /* Verifica que el caracter inválido no haya aparecido, para no repetir el mensaje de error. */
                    if (strchr(ident, *argv[i]) == NULL)
                    {
                        cod_err[n_err] = 1;
                        msj_err[n_err] = argv[i];
                        n_err++;
                    }
                    ident[i - 1] = *argv[i];
                    /* Avanza a la siguiente iteración para no mostrar ningún otro mensaje de error relativo a un
                     * identificador inválido.
                     */
                    continue;
                }
                /* Se verifica que el segundo caracter sea el signo igual. */
                if (*(argv[i] + 1) != '=')
                {
                    cod_err[n_err] = 2;
                    msj_err[n_err] = argv[i];
                    n_err++;
                }
                /* Se verifica que no existan identificadores duplicados. */
                if ((strchr(ident, *argv[i]) != NULL) && (strchr(ident_dup0, *argv[i]) == NULL))
                {
                    cod_err[n_err] = 3;
                    msj_err[n_err] = argv[i];
                    n_err++;
                    *ident_dup++ = *argv[i];
                }

                valor = argv[i] + 2;

                switch(*argv[i])
                {
                case 'm':

                    if (strchr(ident, 'm') == NULL) material = valor;
                    break;

                case 'T':
                    temp = atof(valor);
                    break;

                case 'g':
                    tildegamma = atof(valor);
                    break;

                case 't':
                    o_tit = atoi(valor);
                    if ((o_tit < 0) || (o_tit > 1))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la opción de título `t' solo puede ser 0 o 1";
                    }
                    break;

                case 'k':
                    kappa_1 = atof(valor);
                    break;

                case 'n':
                    n = atoi(valor);
                    if ((n < 1) || (n > 2))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la vorticidad `n' solo puede ser 1 o 2";
                    }
                    break;

                case 'r':
                    factor_r = atof(valor);
                    if ((factor_r < 4) || (factor_r > 100))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "el factor `r' de la frontera truncada debe ser un número real positivo: 4 <= r <= 100";
                    }
                    break;

                case 'l':
                    factor_lh = atof(valor);
                    if ((factor_lh < .5) || (factor_lh > 1))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "el factor `l' de los parámetros de orden en el \"bulk\" debe ser un número real positivo: 1/2 <= l <= 1";
                    }
                    break;

                case 'N':
                    ni = atoi(valor);
                    if ((ni < 10) || (ni > 1e4))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la cantidad de puntos debe ser un número entero positivo: 10 <= N <= 10000";
                    }
                    break;

                case 'o':
                    o_calc = atoi(valor);
                    if ((o_calc < 0) || (o_calc > 2))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la opción de cálculo `o' solo puede ser 0, 1 o 2";
                    }
                    break;

                case 'P':
                    n_bar = atoi(valor);
                    if ((n_bar < 10) || (n_bar > 1e4))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la cantidad de puntos debe ser un número entero positivo: 10 <= P <= 10000";
                    }
                    break;

                case 'd':
                    decimales = atoi(valor);

                    if ((decimales < 4) || (decimales > 15))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la cantidad de decimales `d' debe ser un número entero positivo: 4 <= d <= 15";
                    }
                    break;

                case 'i':
                    if ((strcmp(valor, "png") != 0) && (strcmp(valor, "eps") != 0))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "debe indicar un formato de imagen correspondiente a alguna de las dos extensiones admitidas para la salida gráfica (``png'' o ``eps'')";
                    }
                    strcpy(form_arch, valor);
                    break;

                case 'f':
                    if ((strcmp(valor, "dat") != 0) && (strcmp(valor, "txt") != 0))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "Valor fuera de rango: debe indicar un formato de archivo correspondiente a algunos de las dos extensiones admitidas (``dat'' o ``txt''). Si no desea graficar la solución, por favor emplee la extensión ``txt''";
                    }
                    strcpy(ext, valor);
                    break;

                case 'G':
                    o_graf = atoi(valor);
                    if ((o_graf < 0) || (o_graf > 1))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la selección de salida gráfica `G' solo puede ser 0 o 1";
                    }
                    break;

                case 'R':
                    o_res = atoi(valor);
                    if ((o_res < 0) || (o_res > 3))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la opción de resolución de salida gráfica `R' solo puede ser 0, 1, 2 o 3";
                    }
                    break;

                case 'e':
                    escala = atoi(valor);
                    if ((escala < 0) || (escala > 1))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la escala `e' solo puede ser 0 o 1";
                    }
                    break;

                case 'a':
                    ajuste = atoi(valor);
                    if ((ajuste < 0) || (ajuste > 1))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "el reajuste `a' solo puede ser 0 o 1";
                    }
                    break;

                default:
                    break;
                }
                ident[i - 1] = *argv[i];
            }

            if ((o_calc < 1) && (strchr(ident, 'P') != NULL))
            {
                cod_err[n_err] = 5;
                msj_err[n_err++] = "el argumento `P' se emplea solo en las opciones de barridos";
            }
            if (strchr(ident, 'T') == NULL)
            {
                cod_err[n_err] = 7;
                msj_err[n_err++] = "El argumento temperatura no es opcional. Debe especificar el valor de la temperatura a través del identificador `T'";
            }

            /* Verifica que el usuario haya introducido el argumento material. */
            if (strchr(ident, 'm') != NULL)
            {
                /** Get superconductor parameters. */
                sc = obtener_parametros(material);

                if (sc == NULL)
                {
                    cod_err[n_err] = 0;
                    msj_err[n_err++] = material;
                }
                else
                    /* Validación de la temperatura. */
                    if ((temp < 0) || (temp > *sc->temp_c) || (escala*temp > 1))
                    {
                        cod_err[n_err] = 4;
                        msj_err[n_err++] = "la temperatura debe ser un número real positivo o cero y no mayor a la temperatura crítica del material: 0 <= T/T_c <= 1";
                    }
            }
            else
            {
                cod_err[n_err] = 7;
                msj_err[n_err++] = "El argumento material no es opcional. Debe especificar el material a través del identificador `m'";
            }

            /* Impresión de mensajes de error e interrupción del programa si corresponde. */
            for(i = 0; i < n_err; i++) error(cod_err[i], msj_err[i]);
            if (n_err > 0)
            {
                printf("\n\n");
                exit(1);
            }
        }
        /* Fin de bloque. Interpretación de comandos. */

        /* Ajustes de los parámetros. */
        if (escala > 0) temp *= *sc->temp_c;
        max_err = pow(10, -decimales);

        /* Variables auxiliares para optimizar los cálculos. */
        temp_c = *sc->temp_c++;
        temp_c1 = *sc->temp_c++;
        temp_c2 = *sc->temp_c++;
        w = sc->w;
        gamma_a = sc->gamma_a;
        cte_temp = pow(temp_c1/temp_c2, 2);
        q0 = temp_c1*temp_c2/(temp_c - temp_c1)/(temp_c - temp_c2);
        k = 1 - 4*gamma_a*gamma_a*q0;

        if ((k < 0) || (fabs(k - 1) < 1e-15))
        {
            error(6, "verifique los parámetros del material superconductor en la biblioteca compuestos.txt");
            exit(1);
        }
        k = sqrt(k);


        /* Determinación de la menor temperatura ($T_0$) a partir de la cual las soluciones son válidas. */
        beta = (1 + k)/(1 - k)*cte_temp; /* Simplificación de (-1 - k)/fabs(-1 + k)*cte_temp. */
        m = beta*w;
        temp_ajuste = temp;
        h_bar = (temp_c - temp)/n_bar;
        for (i = 0; i < n_bar; i++)
        {
            ctetemp_c1 = 1 - temp_ajuste/temp_c1;
            ctetemp_c2 = 1 - temp_ajuste/temp_c2;
            s = copysign(1, -ctetemp_c1);
            aux = (1 - k)*fabs(ctetemp_c1); /* Simplificación de fabs((-1 + k)*ctetemp_c1). */
            alpha = -(1 + k)*ctetemp_c2/aux; /* Simplificación de (-1 - k)*ctetemp_c2/aux. */
            gamma = 2*gamma_a/aux;
            if (bulk(s, alpha, beta, gamma, tildegamma, m, &f1oo, &f2oo, &lambda_eff) < 1) break;
            temp_ajuste += h_bar;
        }
        if (i == n_bar)
        {
            error(6, "no se puede determinar la solución. Verifique que los valores de los parámetros sean correctos");
            exit(1);
        }
        if (temp < temp_ajuste)
        {
            if (ajuste > 0)
            {
                printf("\nReajuste de temperatura: T = %g K ---> T = %g K.\n", temp, temp_ajuste);
                temp = temp_ajuste;
            }
            else if (o_calc < 1)
            {
                printf("\nLos parámetros introducidos requieren una temperatura mínima de %g K, para que exista una solución válida.\n", temp_ajuste);
                exit(1);
            }
        }
        temp0 = temp;

        lh_mayor = lambda_eff;
        if (kappa_1 < 1) lh_mayor /= kappa_1;

        /* Reserva de memoria. */
        tam_sol = (ni + 2)*sizeof(double);
        f1 = malloc(tam_sol);
        f2 = malloc(tam_sol);
        a = malloc(tam_sol);
        /* Para barridos de temperaturas */
        if (o_calc > 0)
        {
            n_res = 6;
            tam_sol = n_bar*sizeof(double);
            lambda_eff_var = malloc(tam_sol);
            lambda_var = malloc(tam_sol);
            lh1_var = malloc(tam_sol);
            lh2_var = malloc(tam_sol);
            la_var = malloc(tam_sol);
            f1oo_var = malloc(tam_sol);
            f2oo_var = malloc(tam_sol);
            df10_var = malloc(tam_sol);
            df20_var = malloc(tam_sol);
            B0_var = malloc(tam_sol);
            gamma_var = malloc(tam_sol);
            otra_var = malloc(tam_sol);
            es_invalido = malloc(n_bar*sizeof(bool));
        }
        else n_res = 3;
        archivo = malloc(n_res*sizeof(*archivo));
        id_arch = malloc(n_res*sizeof(*id_arch));
        info_gr = malloc(n_res*sizeof(*info_gr));
        formato = malloc(n_res*sizeof(*formato));

        ctetemp_c1 = 1 - temp/temp_c1;
        ctetemp_c2 = 1 - temp/temp_c2;
        s = copysign(1, -ctetemp_c1);
        aux = (1 - k)*fabs(ctetemp_c1); /* Simplificación de fabs((-1 + k)*ctetemp_c1). */
        alpha = -(1 + k)*ctetemp_c2/aux; /* Simplificación de (-1 - k)*ctetemp_c2/aux. */
        gamma = 2*gamma_a/aux;

        /* Definición del paso y el número de cálculos realizados. */
        switch (o_calc)
        {
        case 0:
            n_bar = 1;
            break;
        case 1:
            h_bar = (temp_c - temp)/n_bar;
            printf("T = %g, T0 = %g, h_bar = %g, n_bar = %d", temp, temp_ajuste, h_bar, n_bar);
            temp0 = temp -= h_bar;
            break;
        default:
            h_bar = (3 - tildegamma)/n_bar;
            tildegamma0 = tildegamma -= h_bar;
            break;
        }

        n_val = 0;

        if (o_calc > 0) printf("\nProgreso:\n0 %%\n");

        /** Resolución del problema. */
        for (i = 0; i < n_bar; i++)
        {
            if (o_calc > 0)
            {
                if ((prctje = 100*(i + 1)/n_bar) > prctje_ant)
                {
                    printf("%d %%\n", prctje);
                    prctje_ant = prctje;
                }
                if (o_calc < 2)
                {
                    temp += h_bar;
                    ctetemp_c1 = 1 - temp/temp_c1;
                    ctetemp_c2 = 1 - temp/temp_c2;
                    s = copysign(1, -ctetemp_c1);
                    aux = (1 - k)*fabs(ctetemp_c1);
                    alpha = (-1 - k)*ctetemp_c2/aux;
                    gamma = 2*gamma_a/aux;
                }
                else if (o_calc < 3) tildegamma += h_bar;
            }

            /* Parámetros constantes de la ecuación diferencial. */
            coef[0] = s;
            coef[1] = alpha;
            coef[2] = beta;
            coef[3] = gamma;
            coef[4] = tildegamma;
            coef[5] = m;
            coef[6] = kappa_1;
            coef[7] = n;

            /* Condiciones de frontera en el límite de la fase volumétrica. */
            if (bulk(s, alpha, beta, gamma, tildegamma, m, &f1oo, &f2oo, &lambda_eff) > 0)
            {
                if (o_calc > 0)
                {
                    es_invalido[i] = true;
                    continue;
                }
                else
                {
                    printf("\nT = %g\n", temp);
                    error(6, "no se puede determinar la solución. Verifique que los valores de los parámetros sean correctos");
                    exit(1);
                }
            }

            cf[0] = f1oo;
            cf[1] = f2oo;
            cf[2] = 1;

            /* Ajuste de $r_\infty$ y $h$ (punto fijo). */
            for (j = 0; j < max_iter; j++)
            {
                h = factor_r*lh_mayor/(ni + 1);

                solve(coef, cf, ni, h, max_err, max_iter, f1, f2, a);

                lh1 = bin_search(f1, ni, h, factor_lh*f1oo, true);
                lh2 = bin_search(f2, ni, h, factor_lh*f2oo, true);
                la = bin_search(a, ni, h, factor_lh, true);

                if (lh1 > lh2)
                    if (lh1 > la)
                        lh_mayor = lh1;
                    else
                        lh_mayor = la;
                else if (lh2 > la)
                    lh_mayor = lh2;
                else
                    lh_mayor = la;

                if (fabs(1 - lh_mayor_ant/lh_mayor)*100 < approx_err) break;

                lh_mayor_ant = lh_mayor;
            }
            r_infty = factor_r*lh_mayor;
            h = r_infty/(ni + 1);

            if (solve(coef, cf, ni, h, max_err, max_iter, f1, f2, a) > 0)
            {
                if (o_calc > 0)
                {
                    es_invalido[i] = true;
                    continue;
                }
                else
                {
                    error(6, "no se puede determinar la solución. Verifique que los valores de los parámetros sean correctos");
                    exit(1);
                }
            }

            /* Cálculo de $B_0$ y $\lambda$. */
            B = campo(a, n/kappa_1, ni, h);
            lambda = bin_search(B, ni, h, *B/exp(1), false);

            aux = 12*h;
            df10 = (-25*f1[0] + 48*f1[1] - 36*f1[2] + 16*f1[3] - 3*f1[4])/aux;
            df20 = (-25*f2[0] + 48*f2[1] - 36*f2[2] + 16*f2[3] - 3*f2[4])/aux;

            /* Cálculos de $L_{h_1}$ y $L_{h_2}$. */
            lh1 = bin_search(f1, ni, h, factor_lh*f1oo, true);
            lh2 = bin_search(f2, ni, h, factor_lh*f2oo, true);
            la = bin_search(a, ni, h, factor_lh, true);

            if (o_calc > 0)
            {
                lambda_eff_var[i] = lambda_eff;
                lambda_var[i] = lambda;
                lh1_var[i] = lh1;
                lh2_var[i] = lh2;
                la_var[i] = la;
                B0_var[i] = *B;
                f1oo_var[i] = f1oo;
                f2oo_var[i] = f2oo;
                df10_var[i] = df10;
                df20_var[i] = df20;
                gamma_var[i] = gamma;
                otra_var[i] = lambda/lambda_eff*(1 - temp/temp_c);
                free(B);
            }
            n_val++;
        }

        /** Directorio y archivos de salida. */
        if (o_calc < 1)
        {
            dir_prefix = "Soluciones";
            sufijo = "GL";
        }
        else
        {
            dir_prefix = "Barridos";
            if (o_calc < 2) sufijo = "temp";
            else sufijo = "tildegamma";
        }

        /* Reserva de memoria para el los nombres de los archivos y el directorio. */
        dir = malloc(strlen(dir_prefix) + strlen(sufijo) + 2);
        nom_arch = malloc(strlen(material) + strlen(sufijo) + 2);

        /* Nombre base del directorio. */
        strcpy(dir, dir_prefix);
        strcat(strcat(dir, "_"), sufijo);
        if (stat(dir, &st) == -1) mkdir(dir, 0700);

        /* Nombre base del archivo. */
        strcpy(nom_arch, material);
        strcpy(nom_arch, material);
        strcat(strcat(nom_arch, "_"), sufijo);

        /* Formato de archivo de salida. */
        if (strcmp(ext, "dat") == 0)
        {
            nom_soft_graf = "Grace";
            char_coment = '#';
        }
        else
        {
            nom_soft_graf = "Texto";
            char_coment = '*';
        }

        if (o_calc < 1)
        {
            info_gr[0].n = 4;
            info_gr[0].etq_x = "r/\\lambda_1(0)";
            info_gr[0].etq_y = "f_1, f_2, B";
            id_arch[0] = "_f1,f2,B";
            info_gr[0].leyenda = malloc((info_gr[0].n - 1)*sizeof(*info_gr[0].leyenda));
            info_gr[0].leyenda[0] = "f_1";
            info_gr[0].leyenda[1] = "f_2";
            info_gr[0].leyenda[2] = "B";

            info_gr[1].n = 2;
            info_gr[1].etq_x = "r/\\lambda_1(0)";
            info_gr[1].etq_y = "f_2/f_1";
            id_arch[1] = "_f2_div_f1";
            info_gr[1].leyenda = malloc((info_gr[1].n - 1)*sizeof(*info_gr[1].leyenda));
            info_gr[1].leyenda[0] = "f_2/f_1";

            info_gr[2].n = 3;
            info_gr[2].etq_x = "r/\\lambda_1(0)";
            info_gr[2].etq_y = "b_1, b_2";
            id_arch[2] = "_b1,b2";
            info_gr[2].leyenda = malloc((info_gr[2].n - 1)*sizeof(*info_gr[2].leyenda));
            info_gr[2].leyenda[0] = "b_1";
            info_gr[2].leyenda[1] = "b_2";
        }
        else
        {
            if (o_calc < 2)
            {
                if (escala < 1) var_ind = "T";
                else var_ind = "T / T_c";
            }
            else var_ind = "\\tilde\\gamma";

            info_gr[0].n = 4;
            info_gr[0].etq_x = var_ind;
            info_gr[0].etq_y = "\\lamda, L_{h_1}, L_{h_2}";
            id_arch[0] = "_l,lh1,lh2";
            info_gr[0].leyenda = malloc((info_gr[0].n - 1)*sizeof(*info_gr[0].leyenda));
            info_gr[0].leyenda[0] = "\\lambda";
            info_gr[0].leyenda[1] = "L_{h_1}";
            info_gr[0].leyenda[2] = "L_{h_2}";

            info_gr[1].n = 2;
            info_gr[1].etq_x = var_ind;
            info_gr[1].etq_y = "B_0";
            id_arch[1] = "_B0";
            info_gr[1].leyenda = malloc((info_gr[1].n - 1)*sizeof(*info_gr[1].leyenda));
            info_gr[1].leyenda[0] = "B_0";

            info_gr[2].n = 2;
            info_gr[2].etq_x = var_ind;
            info_gr[2].etq_y = "L_{h_1}/L_{h_2}";
            id_arch[2] = "_lh1_div_lh2";
            info_gr[2].leyenda = malloc((info_gr[2].n - 1)*sizeof(*info_gr[2].leyenda));
            info_gr[2].leyenda[0] = "L_{h_1}/L_{h_2}";

            info_gr[3].n = 3;
            info_gr[3].etq_x = var_ind;
            info_gr[3].etq_y = "f_{1_\\infty}, f_{2_\\infty}";
            id_arch[3] = "_f1oo,f2oo";
            info_gr[3].leyenda = malloc((info_gr[3].n - 1)*sizeof(*info_gr[3].leyenda));
            info_gr[3].leyenda[0] = "f_{1_\\infty}";
            info_gr[3].leyenda[1] = "f_{2_\\infty}";

            info_gr[4].n = 2;
            info_gr[4].etq_x = var_ind;
            info_gr[4].etq_y = "\\eta = f_{2_\\infty}/f_{1_\\infty}";
            id_arch[4] = "_f2oo_div_f1oo";
            info_gr[4].leyenda = malloc((info_gr[4].n - 1)*sizeof(*info_gr[4].leyenda));
            info_gr[4].leyenda[0] = "\\eta = f_{2_\\infty}/f_{1_\\infty}";

            info_gr[5].n = 3;
            info_gr[5].etq_x = var_ind;
            info_gr[5].etq_y = "b_{1_0}, b_{2_0}";
            id_arch[5] = "_b1,b2";
            info_gr[5].leyenda = malloc((info_gr[5].n - 1)*sizeof(*info_gr[5].leyenda));
            info_gr[5].leyenda[0] = "b_{1_0}";
            info_gr[5].leyenda[1] = "b_{2_0}";
        }


        nom_arch_aux = malloc(n_res*sizeof*nom_arch_aux);

        for (i = 0; i < n_res; i++)
        {
            if (i < n_res) formato[i] = formato_datos(info_gr[i].n, decimales);
            nom_arch_aux[i] = malloc(strlen(nom_arch) + strlen(id_arch[i]) + 1);
            strcat(strcpy(nom_arch_aux[i], nom_arch), id_arch[i]);
            archivo[i] = crear_archivo(dir, nom_arch_aux[i], ext, char_coment, material, nom_soft_graf, usuario, o_calc, ni, n_bar, n_val);
        }

        if (o_calc < 1)
        {
            df1 = deriv(f1, 1, ni, h);
            df2 = deriv(f2, 1, ni, h);

            fprintf(*archivo, "%c\n%c Parámetros de interés:\n", char_coment, char_coment);
            fprintf(*archivo, formato_result("$f_{1_\\infty}$",decimales, char_coment), f1oo);
            fprintf(*archivo, formato_result("$f_{2_\\infty}$",decimales, char_coment), f2oo);
            fprintf(*archivo, formato_result("$L_{h_1}$",decimales, char_coment), lh1);
            fprintf(*archivo, formato_result("$L_{h_2}$",decimales, char_coment), lh2);
            fprintf(*archivo, formato_result("$\\lambda$",decimales, char_coment), lambda);
            fprintf(*archivo, formato_result("$B_0$",decimales, char_coment), *B);

            /* $r>0$. */
            for (j = 0; j < ni + 2; j++)
            {
                fprintf(*archivo, *formato, j*h, f1[j], f2[j], B[j]);

                if (isfinite(aux = f2[j]/f1[j])) fprintf(archivo[1], formato[1], j*h, f2[j]/f1[j]);
                else fprintf(archivo[1], "\n%c Error de cálculo: valor indeterminado.", char_coment);

                fprintf(archivo[2], formato[2], j*h, df1[j], df2[j]);
            }
        }

        else
        {
            if (o_calc < 2) x0 = (temp0 + h_bar)/(aux = pow(temp_c, escala));
            else
            {
                aux = 1;
                x0 = tildegamma0 + h_bar;
            }

            for (j = 0; j < n_bar; j++)
            {
                x = x0 + j*h_bar/aux;

                if (!es_invalido[j]) fprintf(*archivo, *formato, x, lambda_var[j], lh1_var[j], lh2_var[j]);
                else fprintf(*archivo, "\n%c Error de cálculo: valor indeterminado.", char_coment);

                if (!es_invalido[j]) fprintf(archivo[1], formato[1], x, B0_var[j]);
                else fprintf(archivo[1], "\n%c Error de cálculo: valor indeterminado.", char_coment);

                if (!es_invalido[j]) fprintf(archivo[2], formato[2], x, lh1_var[j]/lh2_var[j]);
                else fprintf(archivo[2], "\n%c Error de cálculo: valor indeterminado.", char_coment);

                if (!es_invalido[j]) fprintf(archivo[3], formato[3], x, f1oo_var[j], f2oo_var[j]);
                else fprintf(archivo[4], "\n%c Error de cálculo: valor indeterminado.", char_coment);

                if (!es_invalido[j]) fprintf(archivo[4], formato[4], x, f2oo_var[j]/f1oo_var[j]);
                else fprintf(archivo[5], "\n%c Error de cálculo: valor indeterminado.", char_coment);

                if (!es_invalido[j]) fprintf(archivo[5], formato[5], x, df10_var[j], df20_var[j]);
                else fprintf(archivo[5], "\n%c Error de cálculo: valor indeterminado.", char_coment);
            }
        }

        for (i = 0; i < n_res; i++)
        {
            /* Indicador de fin del archivo de datos. */
            if (strcmp(ext, "dat") == 0) fprintf(archivo[i], "\n&");
            fclose(archivo[i]);

            /* Generación del gráfico. */
            info_gr[i].titulo = malloc(50 + strlen(tipo(material)));
            info_gr[i].subtitulo = malloc(255);

            if (o_calc == 1)
            {
                sprintf(info_gr[i].titulo, "Superconductor: %s", tipo(material));
                sprintf(info_gr[i].subtitulo, "%s = %g,  %s = %g,  m = %g,  %s = %g,  n = %d", tipo("\\beta"), beta, tipo("\\tilde\\gamma"), tildegamma, m, tipo("\\kappa_1"), kappa_1, n);
            }
            else
            {
                sprintf(info_gr[i].titulo, "Superconductor: %s @ T = %g K (%g%s)", tipo(material), temp, temp/temp_c, tipo("T_c"));

                if (o_calc < 1)
                    sprintf(info_gr[i].subtitulo, "s = %d,  %s = %g,  %s = %g,  %s = %g,  %s = %g,  m = %g,  %s = %g,  n = %d", s, tipo("\\alpha"), alpha, tipo("\\beta"), beta, tipo("\\gamma"), gamma, tipo("\\tilde\\gamma"), tildegamma, m, tipo("\\kappa_1"), kappa_1, n);
                else
                    sprintf(info_gr[i].subtitulo, "s = %d,  %s = %g,  %s = %g,  %s = %g,  m = %g,  %s = %g,  n = %d", s, tipo("\\alpha"), alpha, tipo("\\beta"), beta, tipo("\\gamma"), gamma, m, tipo("\\kappa_1"), kappa_1, n);

            }
            if (strcmp(ext, "dat") == 0)
            {
                if (strcmp(info_gr[i].etq_y, "f_1, f_2, B") == 0) crear_grafico(info_gr[i], dir, nom_arch_aux[i], o_graf, o_res, o_tit, form_arch, ni);
                else crear_grafico(info_gr[i], dir, nom_arch_aux[i], o_graf, o_res, o_tit, form_arch, o_calc > 0 ? n_val : ni);
            }
        }

        if (o_calc > 0)
        {
            free(B0_var);
            free(lambda_eff_var);
            free(lambda_var);
            free(lh1_var);
            free(lh2_var);
            free(f1oo_var);
            free(f2oo_var);
            free(df10_var);
            free(df20_var);
            free(gamma_var);
            free(otra_var);
            free(es_invalido);
            for (i = 1; i < n_res; i++)
                free(nom_arch_aux[i]);
            free(*nom_arch_aux);
        }
        else
        {
            free(f1);
            free(f2);
            free(a);
            free(B);
        }

        free(cod_err);
        free(ident);
        free(msj_err);

        free(info_gr);
        printf("\n¡Los cálculos fueron realizados con éxito!\n");
    }
    printf("\n");
    return 0;
}
