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
    double temp_c, temp_c1, temp_c2, temp, temp_tune, const_temp_c1, const_temp_c2, const_temp, gamma_a, r_infty, h, h_iter, *f1, *f2, *a, *B, *df1, *df2, f1_infty, f2_infty, x0, x, df10, df20, coef[9], cf[3], alpha, beta, m, gamma, tildegamma, kappa_1, w, k, q0, lambda_eff, lambda, lh1, lh2, la, lh_mayor, lh_mayor_ant = 0, factor_lh, factor_r, approx_err = .1, *lambda_eff_var, *lambda_var, *lh1_var, *lh2_var, *la_var, *f1_infty_var, *f2_infty_var, *df10_var, *df20_var, *B0_var, *gamma_var, *another_param, max_err, aux, temp0, tildegamma0;
    int i, j, n_iter, n_val, ni;
    short s, show_title, calc_opt, plot_soft_open, resolution_opt, scale, temp0_fix, num_of_results, max_iter = 100;
    unsigned short n,  dec_digits, *err_code, num_of_errors, pct, pct_ant = 0;
    char *material, *file_name, **file_name_aux, *sufix, **id_file, *dir, *dir_prefix, plot_file_format[4], char_com, ext[3], *plot_soft_name, *usuario, *var_ind, *valor, **formato, *ident, *ident_dup, *ident_dup0, **err_msg;
    size_t solution_size;
    bool *is_invalid;
    struct stat st = {0};
    superconductor* sc;
    plot_info *p_info;
    FILE **file;

    switch(argc)
    {
        {
        /* Block start - Command Interpreter */
        case 1:

            display_err_msg(5, "debe especificar los argumentos de entrada");
            break;

        default:

#ifdef __linux__
            usuario = getenv("USER");
#elif _WIN32
            usuario = malloc(UNLEN + 1);
            DWORD usuario_len = UNLEN + 1;
            GetUserName(usuario, &usuario_len);
#endif

            /* Setting default values for entry arguments */
            temp = 0;
            tildegamma = .0;
            kappa_1 = 1/sqrt(2);
            n = 1;
            factor_r = 5;
            factor_lh = .95;
            ni = 750;
            calc_opt = 0;
            show_title = 1;
            n_iter = 750;
            dec_digits = 10;
            strcpy(ext, "dat");
            plot_soft_open = 3;
            resolution_opt = 0;
            scale = 0;
            temp0_fix = 0;
            strcpy(plot_file_format, "png");

            err_code = malloc((argc - 1)*sizeof*err_code);
            err_msg = malloc((argc - 1)*sizeof err_msg);
            ident = malloc(argc - 1);

            /* Syntax checking and validation for parameters entered by user.
             * The string "ident" is set to null, otherwise the string will contain
             * 'garbage' data, causing bugs.
             */
            strncpy(ident, "", argc);
            ident_dup = malloc(argc - 1);
            ident_dup0 = ident_dup;
            num_of_errors = 0;

            for(i = 1; i < argc; i++)
            {
                /*
                 * First character is checked to be a valid identifier
                 * m : material
                 * T : temperature
                 * g : value of $\tilde\gamma$
                 * k : value of $\kappa_1$
                 * n : vorticity
                 * r : multiple of $r_\infty$.
                 * l : multiple of $L_h$.
                 * N : number of points (finite difference method).
                 * o : calculation option (GL solution -> 0, temperature variation -> 1, $\tilde\gamma$ variation -> 2).
                 * P : number of points (temperature or $\tilde\gamma$ variation)
                 * d : number of decimal digits
                 * i : output image file format
                 * f : plotting software file format
                 * G : specifies whether the plot software process remains after getting the results
                 * R : resolution option (800x600 -> 0, 1024x768 -> 1, 1920x1080 -> 2, 2560x2048 -> 3).
                 * t : specifies whether a title is shown at top of the graphs.
                 * e : specifies whether $T$ is scaled to $T_c$.
                 * a : specifies whether initial temperature is adjusted to the first valid solution.
                 */

                if (strchr("mTgknrlNoPdifGRtea", *argv[i]) == NULL)
                {
                    /* Checking whether an invalid character is repeated, in order to not displaying the same error message */
                    if (strchr(ident, *argv[i]) == NULL)
                    {
                        err_code[num_of_errors] = 1;
                        err_msg[num_of_errors++] = argv[i];
                    }
                    ident[i - 1] = *argv[i];
                    /* Move to next iteration so it will not show any other error message related to an invalid identifier */
                    continue;
                }
                /* Checking whether the second character (after an identifier) is the equal "=" symbol */
                if (*(argv[i] + 1) != '=')
                {
                    err_code[num_of_errors] = 2;
                    err_msg[num_of_errors++] = argv[i];
                }
                /* It is verified that there are no duplicate identifiers */
                if ((strchr(ident, *argv[i]) != NULL) && (strchr(ident_dup0, *argv[i]) == NULL))
                {
                    err_code[num_of_errors] = 3;
                    err_msg[num_of_errors++] = argv[i];
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
                    show_title = atoi(valor);
                    if ((show_title < 0) || (show_title > 1))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la opción de título `t' solo puede ser 0 o 1";
                    }
                    break;

                case 'k':
                    kappa_1 = atof(valor);
                    break;

                case 'n':
                    n = atoi(valor);
                    if ((n < 1) || (n > 2))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la vorticidad `n' solo puede ser 1 o 2";
                    }
                    break;

                case 'r':
                    factor_r = atof(valor);
                    if ((factor_r < 4) || (factor_r > 100))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "el factor `r' de la frontera truncada debe ser un número real positivo: 4 <= r <= 100";
                    }
                    break;

                case 'l':
                    factor_lh = atof(valor);
                    if ((factor_lh < .5) || (factor_lh > 1))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "el factor `l' de los parámetros de orden en el \"bulk\" debe ser un número real positivo: 1/2 <= l <= 1";
                    }
                    break;

                case 'N':
                    ni = atoi(valor);
                    if ((ni < 10) || (ni > 1e4))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la cantidad de puntos debe ser un número entero positivo: 10 <= N <= 10000";
                    }
                    break;

                case 'o':
                    calc_opt = atoi(valor);
                    if ((calc_opt < 0) || (calc_opt > 2))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la opción de cálculo `o' solo puede ser 0, 1 o 2";
                    }
                    break;

                case 'P':
                    n_iter = atoi(valor);
                    if ((n_iter < 10) || (n_iter > 1e4))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la cantidad de puntos debe ser un número entero positivo: 10 <= P <= 10000";
                    }
                    break;

                case 'd':
                    dec_digits = atoi(valor);

                    if (( dec_digits < 4) || ( dec_digits > 15))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la cantidad de decimales `d' debe ser un número entero positivo: 4 <= d <= 15";
                    }
                    break;

                case 'i':
                    if ((strcmp(valor, "png") != 0) && (strcmp(valor, "eps") != 0))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "debe indicar un formato de imagen correspondiente a alguna de las dos extensiones admitidas para la salida gráfica (``png'' o ``eps'')";
                    }
                    strcpy(plot_file_format, valor);
                    break;

                case 'f':
                    if ((strcmp(valor, "dat") != 0) && (strcmp(valor, "txt") != 0))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "Valor fuera de rango: debe indicar un formato de archivo correspondiente a algunos de las dos extensiones admitidas (``dat'' o ``txt''). Si no desea graficar la solución, por favor emplee la extensión ``txt''";
                    }
                    strcpy(ext, valor);
                    break;

                case 'G':
                    plot_soft_open = atoi(valor);
                    if ((plot_soft_open < 0) || (plot_soft_open > 1))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la selección de salida gráfica `G' solo puede ser 0 o 1";
                    }
                    break;

                case 'R':
                    resolution_opt = atoi(valor);
                    if ((resolution_opt < 0) || (resolution_opt > 3))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la opción de resolución de salida gráfica `R' solo puede ser 0, 1, 2 o 3";
                    }
                    break;

                case 'e':
                    scale = atoi(valor);
                    if ((scale < 0) || (scale > 1))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la escala `e' solo puede ser 0 o 1";
                    }
                    break;

                case 'a':
                    temp0_fix = atoi(valor);
                    if ((temp0_fix < 0) || (temp0_fix > 1))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "el reajuste `a' solo puede ser 0 o 1";
                    }
                    break;

                default:
                    break;
                }
                ident[i - 1] = *argv[i];
            }

            if ((calc_opt < 1) && (strchr(ident, 'P') != NULL))
            {
                err_code[num_of_errors] = 5;
                err_msg[num_of_errors++] = "el argumento `P' se emplea solo en las opciones de barridos";
            }
            if (strchr(ident, 'T') == NULL)
            {
                err_code[num_of_errors] = 7;
                err_msg[num_of_errors++] = "El argumento temperatura no es opcional. Debe especificar el valor de la temperatura a través del identificador `T'";
            }

            /* Check whether material has been entered */
            if (strchr(ident, 'm') != NULL)
            {
                /** Get superconductor parameters. */
                sc = get_sc_parameters(material);

                if (sc == NULL)
                {
                    err_code[num_of_errors] = 0;
                    err_msg[num_of_errors++] = material;
                }
                else
                    /* Temperature validation. */
                    if ((temp < 0) || (temp > *sc->temp_c) || (scale*temp > 1))
                    {
                        err_code[num_of_errors] = 4;
                        err_msg[num_of_errors++] = "la temperatura debe ser un número real positivo o cero y no mayor a la temperatura crítica del material: 0 <= T/T_c <= 1";
                    }
            }
            else
            {
                err_code[num_of_errors] = 7;
                err_msg[num_of_errors++] = "El argumento material no es opcional. Debe especificar el material a través del identificador `m'";
            }

            /* Print error message and finish execution whether something is wrong */
            for(i = 0; i < num_of_errors; i++) display_err_msg(err_code[i], err_msg[i]);
            if (num_of_errors > 0)
            {
                printf("\n\n");
                exit(1);
            }
        }
        /* End of block. Command interpreter */

        /* Setting parameters */
        if (scale > 0) temp *= *sc->temp_c;
        max_err = pow(10, - dec_digits);

        /* Auxiliary variables to optimize performance */
        temp_c = *sc->temp_c++;
        temp_c1 = *sc->temp_c++;
        temp_c2 = *sc->temp_c++;
        w = sc->w;
        gamma_a = sc->gamma_a;
        const_temp = pow(temp_c1/temp_c2, 2);
        q0 = temp_c1*temp_c2/(temp_c - temp_c1)/(temp_c - temp_c2);
        k = 1 - 4*gamma_a*gamma_a*q0;

        if ((k < 0) || (fabs(k - 1) < 1e-15))
        {
            display_err_msg(6, "verifique los parámetros del material superconductor en la biblioteca sc_materials.txt");
            exit(1);
        }
        k = sqrt(k);


        /* Lowest temperature ($T_0$) is calculated, from which the solutions are valid */
        beta = (1 + k)/(1 - k)*const_temp; /* Simplification of (-1 - k)/fabs(-1 + k)*const_temp. */
        m = beta*w;
        temp_tune = temp;
        h_iter = (temp_c - temp)/n_iter;
        for (i = 0; i < n_iter; i++)
        {
            const_temp_c1 = 1 - temp_tune/temp_c1;
            const_temp_c2 = 1 - temp_tune/temp_c2;
            s = copysign(1, -const_temp_c1);
            aux = (1 - k)*fabs(const_temp_c1); /* Simplification of fabs((-1 + k)*const_temp_c1). */
            alpha = -(1 + k)*const_temp_c2/aux; /* Simplification of (-1 - k)*const_temp_c2/aux. */
            gamma = 2*gamma_a/aux;
            if (bulk(s, alpha, beta, gamma, tildegamma, m, &f1_infty, &f2_infty, &lambda_eff) < 1) break;
            temp_tune += h_iter;
        }
        if (i == n_iter)
        {
            display_err_msg(6, "no se puede determinar la solución. Verifique que los valores de los parámetros sean correctos");
            exit(1);
        }
        if (temp < temp_tune)
        {
            if (temp0_fix > 0)
            {
                printf("\nReajuste de temperatura: T = %g K ---> T = %g K.\n", temp, temp_tune);
                temp = temp_tune;
            }
            else if (calc_opt < 1)
            {
                printf("\nLos parámetros introducidos requieren una temperatura mínima de %g K, para que exista una solución válida.\n", temp_tune);
                exit(1);
            }
        }
        temp0 = temp;

        lh_mayor = lambda_eff;
        if (kappa_1 < 1) lh_mayor /= kappa_1;

        /* Allocating memory */
        solution_size = (ni + 2)*sizeof(double);
        f1 = malloc(solution_size);
        f2 = malloc(solution_size);
        a = malloc(solution_size);

        /* Temperature iteration */
        if (calc_opt > 0)
        {
            num_of_results = 6;
            solution_size = n_iter*sizeof(double);
            lambda_eff_var = malloc(solution_size);
            lambda_var = malloc(solution_size);
            lh1_var = malloc(solution_size);
            lh2_var = malloc(solution_size);
            la_var = malloc(solution_size);
            f1_infty_var = malloc(solution_size);
            f2_infty_var = malloc(solution_size);
            df10_var = malloc(solution_size);
            df20_var = malloc(solution_size);
            B0_var = malloc(solution_size);
            gamma_var = malloc(solution_size);
            another_param = malloc(solution_size);
            is_invalid = malloc(n_iter*sizeof(bool));
        }
        else num_of_results = 3;

        file = malloc(num_of_results*sizeof(*file));
        id_file = malloc(num_of_results*sizeof(*id_file));
        p_info = malloc(num_of_results*sizeof(*p_info));
        formato = malloc(num_of_results*sizeof(*formato));

        const_temp_c1 = 1 - temp/temp_c1;
        const_temp_c2 = 1 - temp/temp_c2;
        s = copysign(1, -const_temp_c1);
        aux = (1 - k)*fabs(const_temp_c1); /* Simplification of fabs((-1 + k)*const_temp_c1). */
        alpha = -(1 + k)*const_temp_c2/aux; /* Simplification of (-1 - k)*const_temp_c2/aux. */
        gamma = 2*gamma_a/aux;

        /* Setting the step and the number of iterations required */
        switch (calc_opt)
        {
        case 0:
            n_iter = 1;
            break;
        case 1:
            h_iter = (temp_c - temp)/n_iter;
            printf("T = %g, T0 = %g, h_iter = %g, n_iter = %d", temp, temp_tune, h_iter, n_iter);
            temp0 = temp -= h_iter;
            break;
        default:
            h_iter = (3 - tildegamma)/n_iter;
            tildegamma0 = tildegamma -= h_iter;
            break;
        }

        n_val = 0;

        if (calc_opt > 0) printf("\nProgreso:\n0 %%\n");

        /** Obtaining the numerical solution */
        for (i = 0; i < n_iter; i++)
        {
            if (calc_opt > 0)
            {
                if ((pct = 100*(i + 1)/n_iter) > pct_ant)
                {
                    printf("%d %%\n", pct);
                    pct_ant = pct;
                }
                if (calc_opt < 2)
                {
                    temp += h_iter;
                    const_temp_c1 = 1 - temp/temp_c1;
                    const_temp_c2 = 1 - temp/temp_c2;
                    s = copysign(1, -const_temp_c1);
                    aux = (1 - k)*fabs(const_temp_c1);
                    alpha = (-1 - k)*const_temp_c2/aux;
                    gamma = 2*gamma_a/aux;
                }
                else if (calc_opt < 3) tildegamma += h_iter;
            }

            /* Constant parameters of the coupled differential equations */
            coef[0] = s;
            coef[1] = alpha;
            coef[2] = beta;
            coef[3] = gamma;
            coef[4] = tildegamma;
            coef[5] = m;
            coef[6] = kappa_1;
            coef[7] = n;

            /* Boundary conditions in the bulk */
            if (bulk(s, alpha, beta, gamma, tildegamma, m, &f1_infty, &f2_infty, &lambda_eff) > 0)
            {
                if (calc_opt > 0)
                {
                    is_invalid[i] = true;
                    continue;
                }
                else
                {
                    printf("\nT = %g\n", temp);
                    display_err_msg(6, "no se puede determinar la solución. Verifique que los valores de los parámetros sean correctos");
                    exit(1);
                }
            }

            cf[0] = f1_infty;
            cf[1] = f2_infty;
            cf[2] = 1;

            /* Setting $r_\infty$ and $h$ (fixed point iteration). */
            for (j = 0; j < max_iter; j++)
            {
                h = factor_r*lh_mayor/(ni + 1);

                solve(coef, cf, ni, h, max_err, max_iter, f1, f2, a);

                lh1 = bin_search(f1, ni, h, factor_lh*f1_infty, true);
                lh2 = bin_search(f2, ni, h, factor_lh*f2_infty, true);
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
                if (calc_opt > 0)
                {
                    is_invalid[i] = true;
                    continue;
                }
                else
                {
                    display_err_msg(6, "no se puede determinar la solución. Verifique que los valores de los parámetros sean correctos");
                    exit(1);
                }
            }

            /* Calculating magnetic field ($B_0$) and penetration depth ($\lambda$) */
            B = magnetic_field(a, n/kappa_1, ni, h);
            lambda = bin_search(B, ni, h, *B/exp(1), false);

            aux = 12*h;
            df10 = (-25*f1[0] + 48*f1[1] - 36*f1[2] + 16*f1[3] - 3*f1[4])/aux;
            df20 = (-25*f2[0] + 48*f2[1] - 36*f2[2] + 16*f2[3] - 3*f2[4])/aux;

            /* Calculating healing lengths $L_{h_1}$ and $L_{h_2}$. */
            lh1 = bin_search(f1, ni, h, factor_lh*f1_infty, true);
            lh2 = bin_search(f2, ni, h, factor_lh*f2_infty, true);
            la = bin_search(a, ni, h, factor_lh, true);

            if (calc_opt > 0)
            {
                lambda_eff_var[i] = lambda_eff;
                lambda_var[i] = lambda;
                lh1_var[i] = lh1;
                lh2_var[i] = lh2;
                la_var[i] = la;
                B0_var[i] = *B;
                f1_infty_var[i] = f1_infty;
                f2_infty_var[i] = f2_infty;
                df10_var[i] = df10;
                df20_var[i] = df20;
                gamma_var[i] = gamma;
                another_param[i] = lambda/lambda_eff*(1 - temp/temp_c);
                free(B);
            }
            n_val++;
        }

        /** Output directories and files */
        if (calc_opt < 1)
        {
            dir_prefix = "Soluciones";
            sufix = "GL";
        }
        else
        {
            dir_prefix = "Barridos";
            if (calc_opt < 2) sufix = "temp";
            else sufix = "tildegamma";
        }

        /* Allocating memory for directory and output file names */
        dir = malloc(strlen(dir_prefix) + strlen(sufix) + 2);
        file_name = malloc(strlen(material) + strlen(sufix) + 2);

        /* Output directory name */
        strcpy(dir, dir_prefix);
        strcat(strcat(dir, "_"), sufix);
        if (stat(dir, &st) == -1) mkdir(dir, 0700);

        /* Output file name */
        strcpy(file_name, material);
        strcpy(file_name, material);
        strcat(strcat(file_name, "_"), sufix);

        /* Output file format */
        if (strcmp(ext, "dat") == 0)
        {
            plot_soft_name = "Grace";
            char_com = '#';
        }
        else
        {
            plot_soft_name = "Texto";
            char_com = '*';
        }

        if (calc_opt < 1)
        {
            p_info[0].n = 4;
            p_info[0].lbl_x = "r/\\lambda_1(0)";
            p_info[0].lbl_y = "f_1, f_2, B";
            id_file[0] = "_f1,f2,B";
            p_info[0].key = malloc((p_info[0].n - 1)*sizeof(*p_info[0].key));
            p_info[0].key[0] = "f_1";
            p_info[0].key[1] = "f_2";
            p_info[0].key[2] = "B";

            p_info[1].n = 2;
            p_info[1].lbl_x = "r/\\lambda_1(0)";
            p_info[1].lbl_y = "f_2/f_1";
            id_file[1] = "_f2_div_f1";
            p_info[1].key = malloc((p_info[1].n - 1)*sizeof(*p_info[1].key));
            p_info[1].key[0] = "f_2/f_1";

            p_info[2].n = 3;
            p_info[2].lbl_x = "r/\\lambda_1(0)";
            p_info[2].lbl_y = "b_1, b_2";
            id_file[2] = "_b1,b2";
            p_info[2].key = malloc((p_info[2].n - 1)*sizeof(*p_info[2].key));
            p_info[2].key[0] = "b_1";
            p_info[2].key[1] = "b_2";
        }
        else
        {
            if (calc_opt < 2)
            {
                if (scale < 1) var_ind = "T";
                else var_ind = "T / T_c";
            }
            else var_ind = "\\tilde\\gamma";

            p_info[0].n = 4;
            p_info[0].lbl_x = var_ind;
            p_info[0].lbl_y = "\\lambda, L_{h_1}, L_{h_2}";
            id_file[0] = "_l,lh1,lh2";
            p_info[0].key = malloc((p_info[0].n - 1)*sizeof(*p_info[0].key));
            p_info[0].key[0] = "\\lambda";
            p_info[0].key[1] = "L_{h_1}";
            p_info[0].key[2] = "L_{h_2}";

            p_info[1].n = 2;
            p_info[1].lbl_x = var_ind;
            p_info[1].lbl_y = "B_0";
            id_file[1] = "_B0";
            p_info[1].key = malloc((p_info[1].n - 1)*sizeof(*p_info[1].key));
            p_info[1].key[0] = "B_0";

            p_info[2].n = 2;
            p_info[2].lbl_x = var_ind;
            p_info[2].lbl_y = "L_{h_1}/L_{h_2}";
            id_file[2] = "_lh1_div_lh2";
            p_info[2].key = malloc((p_info[2].n - 1)*sizeof(*p_info[2].key));
            p_info[2].key[0] = "L_{h_1}/L_{h_2}";

            p_info[3].n = 3;
            p_info[3].lbl_x = var_ind;
            p_info[3].lbl_y = "f_{1_\\infty}, f_{2_\\infty}";
            id_file[3] = "_f1_infty,f2_infty";
            p_info[3].key = malloc((p_info[3].n - 1)*sizeof(*p_info[3].key));
            p_info[3].key[0] = "f_{1_\\infty}";
            p_info[3].key[1] = "f_{2_\\infty}";

            p_info[4].n = 2;
            p_info[4].lbl_x = var_ind;
            p_info[4].lbl_y = "\\eta = f_{2_\\infty}/f_{1_\\infty}";
            id_file[4] = "_f2_infty_div_f1_infty";
            p_info[4].key = malloc((p_info[4].n - 1)*sizeof(*p_info[4].key));
            p_info[4].key[0] = "\\eta = f_{2_\\infty}/f_{1_\\infty}";

            p_info[5].n = 3;
            p_info[5].lbl_x = var_ind;
            p_info[5].lbl_y = "b_{1_0}, b_{2_0}";
            id_file[5] = "_b1,b2";
            p_info[5].key = malloc((p_info[5].n - 1)*sizeof(*p_info[5].key));
            p_info[5].key[0] = "b_{1_0}";
            p_info[5].key[1] = "b_{2_0}";
        }

        file_name_aux = malloc(num_of_results*sizeof*file_name_aux);

        for (i = 0; i < num_of_results; i++)
        {
            if (i < num_of_results) formato[i] = data_format(p_info[i].n, dec_digits);
            file_name_aux[i] = malloc(strlen(file_name) + strlen(id_file[i]) + 1);
            strcat(strcpy(file_name_aux[i], file_name), id_file[i]);
            file[i] = create_file(dir, file_name_aux[i], ext, char_com, material, plot_soft_name, usuario, calc_opt, ni, n_iter, n_val);
        }

        if (calc_opt < 1)
        {
            df1 = num_diff(f1, 1, ni, h);
            df2 = num_diff(f2, 1, ni, h);

            fprintf(*file, "%c\n%c Parámetros de interés:\n", char_com, char_com);
            fprintf(*file, solution_data_format("$f_{1_\\infty}$", dec_digits, char_com), f1_infty);
            fprintf(*file, solution_data_format("$f_{2_\\infty}$", dec_digits, char_com), f2_infty);
            fprintf(*file, solution_data_format("$L_{h_1}$", dec_digits, char_com), lh1);
            fprintf(*file, solution_data_format("$L_{h_2}$", dec_digits, char_com), lh2);
            fprintf(*file, solution_data_format("$\\lambda$", dec_digits, char_com), lambda);
            fprintf(*file, solution_data_format("$B_0$", dec_digits, char_com), *B);

            /* $r>0$. */
            for (j = 0; j < ni + 2; j++)
            {
                fprintf(*file, *formato, j*h, f1[j], f2[j], B[j]);

                if (isfinite(aux = f2[j]/f1[j])) fprintf(file[1], formato[1], j*h, f2[j]/f1[j]);
                else fprintf(file[1], "\n%c Error de cálculo: valor indeterminado.", char_com);

                fprintf(file[2], formato[2], j*h, df1[j], df2[j]);
            }
        }

        else
        {
            if (calc_opt < 2) x0 = (temp0 + h_iter)/(aux = pow(temp_c, scale));
            else
            {
                aux = 1;
                x0 = tildegamma0 + h_iter;
            }

            for (j = 0; j < n_iter; j++)
            {
                x = x0 + j*h_iter/aux;

                if (!is_invalid[j]) fprintf(*file, *formato, x, lambda_var[j], lh1_var[j], lh2_var[j]);
                else fprintf(*file, "\n%c Error de cálculo: valor indeterminado.", char_com);

                if (!is_invalid[j]) fprintf(file[1], formato[1], x, B0_var[j]);
                else fprintf(file[1], "\n%c Error de cálculo: valor indeterminado.", char_com);

                if (!is_invalid[j]) fprintf(file[2], formato[2], x, lh1_var[j]/lh2_var[j]);
                else fprintf(file[2], "\n%c Error de cálculo: valor indeterminado.", char_com);

                if (!is_invalid[j]) fprintf(file[3], formato[3], x, f1_infty_var[j], f2_infty_var[j]);
                else fprintf(file[4], "\n%c Error de cálculo: valor indeterminado.", char_com);

                if (!is_invalid[j]) fprintf(file[4], formato[4], x, f2_infty_var[j]/f1_infty_var[j]);
                else fprintf(file[5], "\n%c Error de cálculo: valor indeterminado.", char_com);

                if (!is_invalid[j]) fprintf(file[5], formato[5], x, df10_var[j], df20_var[j]);
                else fprintf(file[5], "\n%c Error de cálculo: valor indeterminado.", char_com);
            }
        }

        for (i = 0; i < num_of_results; i++)
        {
            /* End of file marker */
            if (strcmp(ext, "dat") == 0) fprintf(file[i], "\n&");
            fclose(file[i]);

            /* Generating plot */
            p_info[i].title = malloc(50 + strlen(typo_convert(material)));
            p_info[i].subtitle = malloc(255);

            if (calc_opt == 1)
            {
                sprintf(p_info[i].title, "Superconductor: %s", typo_convert(material));
                sprintf(p_info[i].subtitle, "%s = %g,  %s = %g,  m = %g,  %s = %g,  n = %d", typo_convert("\\beta"), beta, typo_convert("\\tilde\\gamma"), tildegamma, m, typo_convert("\\kappa_1"), kappa_1, n);
            }
            else
            {
                sprintf(p_info[i].title, "Superconductor: %s @ T = %g K (%g%s)", typo_convert(material), temp, temp/temp_c, typo_convert("T_c"));

                if (calc_opt < 1)
                    sprintf(p_info[i].subtitle, "s = %d,  %s = %g,  %s = %g,  %s = %g,  %s = %g,  m = %g,  %s = %g,  n = %d", s, typo_convert("\\alpha"), alpha, typo_convert("\\beta"), beta, typo_convert("\\gamma"), gamma, typo_convert("\\tilde\\gamma"), tildegamma, m, typo_convert("\\kappa_1"), kappa_1, n);
                else
                    sprintf(p_info[i].subtitle, "s = %d,  %s = %g,  %s = %g,  %s = %g,  m = %g,  %s = %g,  n = %d", s, typo_convert("\\alpha"), alpha, typo_convert("\\beta"), beta, typo_convert("\\gamma"), gamma, m, typo_convert("\\kappa_1"), kappa_1, n);

            }
            if (strcmp(ext, "dat") == 0)
            {
                if (strcmp(p_info[i].lbl_y, "f_1, f_2, B") == 0) create_plot(p_info[i], dir, file_name_aux[i], plot_soft_open, resolution_opt, show_title, plot_file_format, ni);
                else create_plot(p_info[i], dir, file_name_aux[i], plot_soft_open, resolution_opt, show_title, plot_file_format, calc_opt > 0 ? n_val : ni);
            }
        }

        if (calc_opt > 0)
        {
            free(B0_var);
            free(lambda_eff_var);
            free(lambda_var);
            free(lh1_var);
            free(lh2_var);
            free(f1_infty_var);
            free(f2_infty_var);
            free(df10_var);
            free(df20_var);
            free(gamma_var);
            free(another_param);
            free(is_invalid);
            for (i = 1; i < num_of_results; i++)
                free(file_name_aux[i]);
            free(*file_name_aux);
        }
        else
        {
            free(f1);
            free(f2);
            free(a);
            free(B);
        }

        free(err_code);
        free(ident);
        free(err_msg);

        free(p_info);
        printf("\n¡Los cálculos fueron realizados con éxito!\n");
    }
    printf("\n");
    return 0;
}
