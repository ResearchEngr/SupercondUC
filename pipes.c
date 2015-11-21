#ifdef __linux__
#include <grace_np.h>
#endif // __linux__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pipes.h"

/*
 * Convert a \TeX typo contained in *text parameter into a Grace typo, in order to be represented correctly on the plot
 */
char* typo_convert(char *text)
{
    char *t;

    if (strcmp(text, "MgB2") == 0)
        t = "MgB\\s2\\N";
    else if (strcmp(text, "OsB2") == 0)
        t = "OsB\\s2\\N";
    else if (strcmp(text, "FeSe1-x") == 0)
        t = "FeSe\\s1-x\\N";
    else if (strcmp(text, "V3Si") == 0)
        t = "V\\s3\\NSi";
    else if (strcmp(text, "b_{1_0}") == 0)
        t = "b\\s1\\s0\\N";
    else if (strcmp(text, "b_{2_0}") == 0)
        t = "b\\s2\\s0\\N";
    else if (strcmp(text, "f_\\xi") == 0)
        t = "f\\s\\xx\\N\\f{}";
    else if (strcmp(text, "\\alpha") == 0)
        t = "\\xa\\f{}";
    else if (strcmp(text, "\\beta") == 0)
        t = "\\xb\\f{}";
    else if (strcmp(text, "\\gamma") == 0)
        t = "\\xg\\N\\f{}";
    else if (strcmp(text, "\\gamma_a") == 0)
        t = "\\xg\\f{}\\sa\\N";
    else if (strcmp(text, "\\gamma(T)") == 0)
        t = "\\xg\\N\\f{}(T)";
    else if (strcmp(text, "\\kappa_1") == 0)
        t = "\\xk\\s1\\N\\f{}";
    else if (strcmp(text, "\\lambda") == 0)
        t = "\\xl\\f{}";
    else if (strcmp(text, "\\lambda_{eff}") == 0)
        t = "\\xl\\f{}\\seff\\N";
    else if (strcmp(text, "L_{h_1}") == 0)
        t = "L\\sh\\s1\\N\\f{}";
    else if (strcmp(text, "L_{h_2}") == 0)
        t = "L\\sh\\s2\\N\\f{}";
    else if (strcmp(text, "L_{h_a}") == 0)
        t = "L\\sh\\sa\\N\\f{}";
    else if (strcmp(text, "\\eta") == 0)
        t = "\\xh\\f{}";
    else if (strcmp(text, "f_1, f_2, B") == 0)
        t = "f\\s1\\N, f\\s2\\N, B";
    else if (strcmp(text, "r/\\lambda_1(0)") == 0)
        t = "r / \\xl\\s1\\N\\f{}(0)";
    else if (strcmp(text, "T_c") == 0)
        t = "T\\sc\\N";
    else if (strcmp(text, "T / T_c") == 0)
        t = "T / T\\sc\\N";
    else if (strcmp(text, "f_1") == 0)
        t = "f\\s1\\N";
    else if (strcmp(text, "f_2") == 0)
        t = "f\\s2\\N";
    else if (strcmp(text, "f_2/f_1") == 0)
        t = "f\\s2\\N/f\\s1\\N";
    else if (strcmp(text, "b_{1_0}") == 0)
        t = "b\\s1\\s0\\N";
    else if (strcmp(text, "b_{2_0}") == 0)
        t = "b\\s2\\s0\\N";
    else if (strcmp(text, "B_0") == 0)
        t = "B\\s0\\N";
    else if (strcmp(text, "\\tilde\\gamma") == 0)
        t = "\\xg\\h{-0.45}\\v{0.45}~\\N\\f{}";
    else if (strcmp(text, "r_\\infty") == 0)
        t = "r\\s\\x\\c%\\C\\N\\f{}";
    else if (strcmp(text, "f_{1_\\infty}") == 0)
        t = "f\\s1\\s\\x\\c%\\C\\N\\f{}";
    else if (strcmp(text, "f_{2_\\infty}") == 0)
        t = "f\\s2\\s\\x\\c%\\C\\N\\f{}";
    else if (strcmp(text, "f_{1_\\infty}/b_{1_0}") == 0)
        t = "(f\\s1\\s\\x\\c%\\C\\f{}\\N\\s/b\\s1\\s0\\N\\s)";
    else if (strcmp(text, "f_{2_\\infty}/b_{2_0}") == 0)
        t = "(f\\s2\\s\\x\\c%\\C\\f{}\\N\\s/b\\s2\\s0\\N\\s)";
    else if (strcmp(text, "\\lambda, L_{h_1}, L_{h_2}") == 0)
        t = "\\xl\\f{}, L\\sh\\s1\\N, L\\sh\\s2\\N\\f{}";
    else if (strcmp(text, "\\lambda/\\xi_1, \\lambda/\\xi_2") == 0)
        t = "\\xl\\N / \\xx\\s1\\N, \\xl\\N / \\xx\\s2\\N\\f{}";
    else if (strcmp(text, "\\lambda/\\xi_1") == 0)
        t = "\\xl\\N / \\xx\\s1\\N\\f{}";
    else if (strcmp(text, "\\lambda/\\xi_2") == 0)
        t = "\\xl\\N / \\xx\\s2\\N\\f{}";
    else if (strcmp(text, "L_{h_1}/L_{h_2}") == 0)
        t = "L\\sh\\s1\\N / L\\sh\\s2\\N\\f{}";
    else if (strcmp(text, "b_{1_0}, b_{2_0}") == 0)
        t = "b\\s1\\s0\\N, b\\s2\\s0\\N";
    else if (strcmp(text, "b_{1_0}/b_{2_0}") == 0)
        t = "b\\s1\\s0\\N/b\\s2\\s0\\N";
    else if (strcmp(text, "b_1") == 0)
        t = "b\\s1\\N";
    else if (strcmp(text, "b_2") == 0)
        t = "b\\s2\\N";
    else if (strcmp(text, "b_1, b_2") == 0)
        t = "b\\s1\\N, b\\s2\\N";
    else if (strcmp(text, "b_1/b_2") == 0)
        t = "b\\s1\\N/b\\s2\\N";
    else if (strcmp(text, "f_{1_\\infty}, f_{2_\\infty}") == 0)
        t = "f\\s1\\s\\x\\c%\\C\\N, \\f{}f\\s2\\s\\x\\c%\\C\\N";
    else if (strcmp(text, "\\eta = f_{2_\\infty}/f_{1_\\infty}") == 0)
        t = "\\xh\\f{} = f\\s2\\s\\x\\c%\\C\\N\\f{} / f\\s1\\s\\x\\c%\\C\\N\\f{}";
    else if (strcmp(text, "B_0") == 0)
        t = "B\\s0\\N";
    else if (strcmp(text, "\\gamma") == 0)
        t = "\\xg\\N";
    else
        t = text;

    return t;
}

/*
 * Get the data bounds commands which Grace will process to represent data on the plot
 */
void plot_bounds(char var, short n, char *min, char *max)
{
    switch (n)
    {
    case 1:
        sprintf(min, "min(s0.%c)", var);
        sprintf(max, "max(s0.%c)", var);
        break;
    case 2:
        sprintf(min, "minof(min(s0.%c),min(s1.%c))", var, var);
        sprintf(max, "maxof(max(s0.%c),max(s1.%c))", var, var);
        break;
    case 3:
        sprintf(min, "minof(minof(min(s0.%c),min(s1.%c)),min(s2.%c))", var, var, var);
        sprintf(max, "maxof(maxof(max(s0.%c),max(s1.%c)),max(s2.%c))", var, var, var);
        break;
    default:
        sprintf(min, "minof(minof(minof(min(s0.%c),min(s1.%c)),min(s2.%c)),min(s3.%c))", var, var, var, var);
        sprintf(max, "maxof(maxof(maxof(max(s0.%c),max(s1.%c)),max(s2.%c)),max(s3.%c))", var, var, var, var);
        break;
    }
}

/* Start Grace and creates a plot */
void create_plot(plot_info p_info, char *dir, char *file_name, short plot_soft_open, short resolution_opt, short show_title, char *plot_file_format, unsigned n_val)
/*
 * Descripción:
 * Esta función ejecuta un proceso Grace y crea un gráfico a partir de la información contenida en los argumentos de entrada.
 *
 * Parámetros:
 * p_info   : estructura compleja de datos que almacena información sobre las propiedades del gráfico.
 * dir       : cadena de caracteres que contiene la ruta relativa donde se almacenan los archivos de salida del sistema.
 * file_name  : cadena de caracteres que contiene el nombre del archivo donde se almacenan los datos que serán representados en la aplicación Grace.
 * plot_soft_open    : opción que permite al usuario decidir si se deja habilitada la aplicación Grace tras finalizar la ejecución del programa.
 * resolution_opt     : opción de resolución (en píxeles) para la representación gráfica de los datos.
 * plot_file_format : cadena de caracteres que define el formato del archivo de imagen de salida (<<png>> o <<eps>>).
 * n_val     : número de datos válidos que se representarán en la aplicación Grace. Este parámetro se emplea para ubicar adecuadamente los símbolos que diferencian las curvas en el gráfico.
 */
{
    short i, n = p_info.n - 1,
             line_style[] = {1, 1, 1, 1},
                            symbol[] = {0, 1, 3, 4},
                                       line_color[] = {2, 15, 4, 1},
                                               symbol_color[] = {2, 15, 4, 1};
    unsigned width, height;
    char *min_x, *max_x, *min_y, *max_y;

    /* Start Grace with a buffer size of 2048 and open the pipe */
#ifdef __linux__
    if (GraceOpen(2048) == -1)
    {
        fprintf(stderr, "No se puede ejecutar Grace.\n");
        exit(1);
    }

    GracePrintf("default font 4");
    GracePrintf("title font 4");
    GracePrintf("subtitle font 4");
    GracePrintf("xaxis label font 4");
    GracePrintf("yaxis label font 4");
    GracePrintf("xaxis ticklabel font 4");
    GracePrintf("yaxis ticklabel font 4");
    GracePrintf("legend font 4");

    GracePrintf("xaxis ticklabel char size .8");
    GracePrintf("yaxis ticklabel char size .8");
    GracePrintf("title size 1.1");
    GracePrintf("subtitle size .8");
    GracePrintf("legend char size .85");

    GracePrintf("xaxis tick default 10");
    GracePrintf("xaxis tick out");
    GracePrintf("yaxis tick default 10");
    GracePrintf("yaxis tick out");

    GracePrintf("xaxis tick major grid on");
    GracePrintf("xaxis tick minor grid on");
    GracePrintf("xaxis tick major linewidth .5");
    GracePrintf("xaxis tick major linestyle 2");
    GracePrintf("xaxis tick minor linewidth .5");
    GracePrintf("xaxis tick minor linestyle 2");
    GracePrintf("xaxis tick place normal");
    GracePrintf("yaxis tick major grid on");
    GracePrintf("yaxis tick minor grid on");
    GracePrintf("yaxis tick major linewidth .5");
    GracePrintf("yaxis tick minor linestyle 2");
    GracePrintf("yaxis tick minor linewidth .5");
    GracePrintf("yaxis tick major linestyle 2");
    GracePrintf("yaxis tick place normal");

    if (show_title > 0)
    {
        GracePrintf("title \"%s\"", p_info.title);
        GracePrintf("subtitle \"%s\"", p_info.subtitle);
    }
    GracePrintf("xaxis label \"%s\"", typo_convert(p_info.lbl_x));
    GracePrintf("yaxis label \"%s\"", typo_convert(p_info.lbl_y));

    for (i = 0; i < n; i++)
    {
        GracePrintf("s%d legend \"%s\"", i, typo_convert(p_info.key[i]));
        GracePrintf("s%d symbol %d", i, symbol[i]);
        GracePrintf("s%d symbol color %d", i, symbol_color[i]);
        GracePrintf("s%d symbol size .85", i);
        GracePrintf("s%d symbol skip %d", i, n_val/50);
        GracePrintf("s%d line linewidth 3", i);
        GracePrintf("s%d line color %d", i, line_color[i]);
        GracePrintf("s%d line linestyle %d", i, line_style[i]);
    }

    switch(resolution_opt)
    {
    case 1:
        width = 1188;
        height = 918;
        break;
    case 2:
        width = 1584;
        height = 1224;
        break;
    case 3:
        width = 2376;
        height = 1836;
        break;
    default:
        width = 792;
        height = 612;
        break;
    }

    min_y = malloc(65);
    max_y = malloc(65);
    min_x = malloc(65);
    max_x = malloc(65);

    plot_bounds('x', n, min_x, max_x);
    plot_bounds('y', n, min_y, max_y);

    GracePrintf("read nxy \"%s/%s.dat\"", dir, file_name);
    GracePrintf("page size %u, %u", width, height);
    GracePrintf("world xmin %s", min_x);
    GracePrintf("world ymin 1.02*%s-.02*%s", min_y, max_y);
    if (strcmp(p_info.lbl_x, "r/\\lambda_1(0)") == 0)
        GracePrintf("world xmax .55*%s", max_x);
    else
        GracePrintf("world xmax %s", max_x);
    GracePrintf("world ymax 1.05*%s-.05*%s", max_y, min_y);
    GracePrintf("legend on");
    GracePrintf("saveall \"%s/%s.agr\"", dir, file_name);

    GracePrintf("print to \"%s/%s.%s\"", dir, file_name, plot_file_format);
    /* Output file format */
    if (strcmp(plot_file_format, "png") == 0)
    {
        GracePrintf("hardcopy device \"PNG\"");
        GracePrintf("device \"PNG\" font antialiasing on");
        GracePrintf("device \"PNG\" op \"compression:9\"");
        /* Optional:
         * GracePrintf("device \"PNG\" op \"transparent:on\"");
         */
    }
    else if (strcmp(plot_file_format, "eps") == 0)
    {
        GracePrintf("hardcopy device \"EPS\"");
        GracePrintf("device \"EPS\" op \"level2\"");
    }
    GracePrintf("print");

    if (plot_soft_open < 1) GraceClose();
    else GraceClosePipe();
    free(min_x);
    free(max_x);
    free(min_y);
    free(max_y);
#endif // __linux__
}
