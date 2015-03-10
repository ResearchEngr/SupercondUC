#ifdef __linux__
#include <grace_np.h>
#endif // __linux__
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pipes.h"

char* tipo(char *texto)
/*
 * Descripción:
 * Esta función retorna un puntero de tipo char, el cual contiene la dirección de memoria de una cadena de caracteres que define el formato tipográfico de un texto, compatible con la aplicación Grace.
 *
 * Parámetro:
 * texto : código en \TeX de un texto que se desea representar en la aplicación Grace.
 */
{
    char *t;

    if (strcmp(texto, "MgB2") == 0)
        t = "MgB\\s2\\N";
    else if (strcmp(texto, "OsB2") == 0)
        t = "OsB\\s2\\N";
    else if (strcmp(texto, "FeSe1-x") == 0)
        t = "FeSe\\s1-x\\N";
    else if (strcmp(texto, "V3Si") == 0)
        t = "V\\s3\\NSi";
    else if (strcmp(texto, "b_{1_0}") == 0)
        t = "b\\s1\\s0\\N";
    else if (strcmp(texto, "b_{2_0}") == 0)
        t = "b\\s2\\s0\\N";
    else if (strcmp(texto, "f_\\xi") == 0)
        t = "f\\s\\xx\\N\\f{}";
    else if (strcmp(texto, "\\alpha") == 0)
        t = "\\xa\\f{}";
    else if (strcmp(texto, "\\beta") == 0)
        t = "\\xb\\f{}";
    else if (strcmp(texto, "\\gamma") == 0)
        t = "\\xg\\N\\f{}";
    else if (strcmp(texto, "\\gamma_a") == 0)
        t = "\\xg\\f{}\\sa\\N";
    else if (strcmp(texto, "\\gamma(T)") == 0)
        t = "\\xg\\N\\f{}(T)";
    else if (strcmp(texto, "\\kappa_1") == 0)
        t = "\\xk\\s1\\N\\f{}";
    else if (strcmp(texto, "\\lambda") == 0)
        t = "\\xl\\f{}";
    else if (strcmp(texto, "\\lambda_{eff}") == 0)
        t = "\\xl\\f{}\\seff\\N";
    else if (strcmp(texto, "L_{h_1}") == 0)
        t = "L\\sh\\s1\\N\\f{}";
    else if (strcmp(texto, "L_{h_2}") == 0)
        t = "L\\sh\\s2\\N\\f{}";
    else if (strcmp(texto, "L_{h_a}") == 0)
        t = "L\\sh\\sa\\N\\f{}";
    else if (strcmp(texto, "\\eta") == 0)
        t = "\\xh\\f{}";
    else if (strcmp(texto, "f_1, f_2, B") == 0)
        t = "f\\s1\\N, f\\s2\\N, B";
    else if (strcmp(texto, "r/\\lambda_1(0)") == 0)
        t = "r / \\xl\\s1\\N\\f{}(0)";
    else if (strcmp(texto, "T_c") == 0)
        t = "T\\sc\\N";
    else if (strcmp(texto, "T / T_c") == 0)
        t = "T / T\\sc\\N";
    else if (strcmp(texto, "f_1") == 0)
        t = "f\\s1\\N";
    else if (strcmp(texto, "f_2") == 0)
        t = "f\\s2\\N";
    else if (strcmp(texto, "f_2/f_1") == 0)
        t = "f\\s2\\N/f\\s1\\N";
    else if (strcmp(texto, "b_{1_0}") == 0)
        t = "b\\s1\\s0\\N";
    else if (strcmp(texto, "b_{2_0}") == 0)
        t = "b\\s2\\s0\\N";
    else if (strcmp(texto, "B_0") == 0)
        t = "B\\s0\\N";
    else if (strcmp(texto, "\\tilde\\gamma") == 0)
        t = "\\xg\\h{-0.45}\\v{0.45}~\\N\\f{}";
    else if (strcmp(texto, "r_\\infty") == 0)
        t = "r\\s\\x\\c%\\C\\N\\f{}";
    else if (strcmp(texto, "f_{1_\\infty}") == 0)
        t = "f\\s1\\s\\x\\c%\\C\\N\\f{}";
    else if (strcmp(texto, "f_{2_\\infty}") == 0)
        t = "f\\s2\\s\\x\\c%\\C\\N\\f{}";
    else if (strcmp(texto, "f_{1_\\infty}/b_{1_0}") == 0)
        t = "(f\\s1\\s\\x\\c%\\C\\f{}\\N\\s/b\\s1\\s0\\N\\s)";
    else if (strcmp(texto, "f_{2_\\infty}/b_{2_0}") == 0)
        t = "(f\\s2\\s\\x\\c%\\C\\f{}\\N\\s/b\\s2\\s0\\N\\s)";
    else if (strcmp(texto, "\\lamda, L_{h_1}, L_{h_2}") == 0)
        t = "\\xl\\f{}, L\\sh\\s1\\N, L\\sh\\s2\\N\\f{}";
    else if (strcmp(texto, "\\lamda/\\xi_1, \\lambda/\\xi_2") == 0)
        t = "\\xl\\N / \\xx\\s1\\N, \\xl\\N / \\xx\\s2\\N\\f{}";
    else if (strcmp(texto, "\\lamda/\\xi_1") == 0)
        t = "\\xl\\N / \\xx\\s1\\N\\f{}";
    else if (strcmp(texto, "\\lamda/\\xi_2") == 0)
        t = "\\xl\\N / \\xx\\s2\\N\\f{}";
    else if (strcmp(texto, "L_{h_1}/L_{h_2}") == 0)
        t = "L\\sh\\s1\\N / L\\sh\\s2\\N\\f{}";
    else if (strcmp(texto, "b_{1_0}, b_{2_0}") == 0)
        t = "b\\s1\\s0\\N, b\\s2\\s0\\N";
    else if (strcmp(texto, "b_{1_0}/b_{2_0}") == 0)
        t = "b\\s1\\s0\\N/b\\s2\\s0\\N";
    else if (strcmp(texto, "b_1") == 0)
        t = "b\\s1\\N";
    else if (strcmp(texto, "b_2") == 0)
        t = "b\\s2\\N";
    else if (strcmp(texto, "b_1, b_2") == 0)
        t = "b\\s1\\N, b\\s2\\N";
    else if (strcmp(texto, "b_1/b_2") == 0)
        t = "b\\s1\\N/b\\s2\\N";
    else if (strcmp(texto, "f_{1_\\infty}, f_{2_\\infty}") == 0)
        t = "f\\s1\\s\\x\\c%\\C\\N, \\f{}f\\s2\\s\\x\\c%\\C\\N";
    else if (strcmp(texto, "\\eta = f_{2_\\infty}/f_{1_\\infty}") == 0)
        t = "\\xh\\f{} = f\\s2\\s\\x\\c%\\C\\N\\f{} / f\\s1\\s\\x\\c%\\C\\N\\f{}";
    else if (strcmp(texto, "B_0") == 0)
        t = "B\\s0\\N";
    else if (strcmp(texto, "\\gamma") == 0)
        t = "\\xg\\N";
    else
        t = texto;

    return t;
}

void limites(char var, short n, char *minimo, char *maximo)
/*
 * Descripción:
 * Esta función tiene como propósito determinar los límites establecidos por el conjunto de datos que definen un gráfico.
 *
 * Parámetros:
 * var    : caracter (`x' o `y') que define si los límites corresponden al eje x o al eje y del gráfico, respectivamente.
 * n      : cantidad de columnas que conforman el conjunto de datos que serán representados en el gráfico. Estos datos incluyen aquellos de la variable independiente y los representados por las curvas del gráfico.
 * minimo : cadena de caracteres que almacena una instrucción para ser enviada al intérprete de órdenes de la aplicación Grace. Esta orden determina el menor de los elementos del conjunto de datos que definen el gráfico.
 * maximo : cadena de caracteres que almacena una instrucción para ser enviada al intérprete de órdenes de la aplicación Grace. Esta orden determina el mayor de los elementos del conjunto de datos que definen el gráfico.
 */
{
    switch (n)
    {
    case 1:
        sprintf(minimo, "min(s0.%c)", var);
        sprintf(maximo, "max(s0.%c)", var);
        break;
    case 2:
        sprintf(minimo, "minof(min(s0.%c),min(s1.%c))", var, var);
        sprintf(maximo, "maxof(max(s0.%c),max(s1.%c))", var, var);
        break;
    case 3:
        sprintf(minimo, "minof(minof(min(s0.%c),min(s1.%c)),min(s2.%c))", var, var, var);
        sprintf(maximo, "maxof(maxof(max(s0.%c),max(s1.%c)),max(s2.%c))", var, var, var);
        break;
    default:
        sprintf(minimo, "minof(minof(minof(min(s0.%c),min(s1.%c)),min(s2.%c)),min(s3.%c))", var, var, var, var);
        sprintf(maximo, "maxof(maxof(maxof(max(s0.%c),max(s1.%c)),max(s2.%c)),max(s3.%c))", var, var, var, var);
        break;
    }
}

void crear_grafico(info_grafico info_gr, char *dir, char *nom_arch, short o_graf, short o_res, short o_tit, char *form_arch, unsigned n_val)
/*
 * Descripción:
 * Esta función ejecuta un proceso Grace y crea un gráfico a partir de la información contenida en los argumentos de entrada.
 *
 * Parámetros:
 * info_gr   : estructura compleja de datos que almacena información sobre las propiedades del gráfico.
 * dir       : cadena de caracteres que contiene la ruta relativa donde se almacenan los archivos de salida del sistema.
 * nom_arch  : cadena de caracteres que contiene el nombre del archivo donde se almacenan los datos que serán representados en la aplicación Grace.
 * o_graf    : opción que permite al usuario decidir si se deja habilitada la aplicación Grace tras finalizar la ejecución del programa.
 * o_res     : opción de resolución (en píxeles) para la representación gráfica de los datos.
 * form_arch : cadena de caracteres que define el formato del archivo de imagen de salida (<<png>> o <<eps>>).
 * n_val     : número de datos válidos que se representarán en la aplicación Grace. Este parámetro se emplea para ubicar adecuadamente los símbolos que diferencian las curvas en el gráfico.
 */
{
    short i, n = info_gr.n - 1,
             estilo_linea[] = {1, 1, 1, 1},
                              simbolo[] = {0, 1, 3, 4},
                                          color_linea[] = {2, 15, 4, 1},
                                                  color_simbolo[] = {2, 15, 4, 1};
    unsigned ancho_pag, largo_pag;
    char *min_x, *max_x, *min_y, *max_y;

    /* Inicia un proceso Grace con un búfer de 2048 bytes y habilita la comunicación con el proceso. */
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

    if (o_tit > 0)
    {
        GracePrintf("title \"%s\"", info_gr.titulo);
        GracePrintf("subtitle \"%s\"", info_gr.subtitulo);
    }
    GracePrintf("xaxis label \"%s\"", tipo(info_gr.etq_x));
    GracePrintf("yaxis label \"%s\"", tipo(info_gr.etq_y));

    for (i = 0; i < n; i++)
    {
        GracePrintf("s%d legend \"%s\"", i, tipo(info_gr.leyenda[i]));
        GracePrintf("s%d symbol %d", i, simbolo[i]);
        GracePrintf("s%d symbol color %d", i, color_simbolo[i]);
        GracePrintf("s%d symbol size .85", i);
        GracePrintf("s%d symbol skip %d", i, n_val/50);
        GracePrintf("s%d line linewidth 3", i);
        GracePrintf("s%d line color %d", i, color_linea[i]);
        GracePrintf("s%d line linestyle %d", i, estilo_linea[i]);
    }

    switch(o_res)
    {
    case 1:
        ancho_pag = 1188;
        largo_pag = 918;
        break;
    case 2:
        ancho_pag = 1584;
        largo_pag = 1224;
        break;
    case 3:
        ancho_pag = 2376;
        largo_pag = 1836;
        break;
    default:
        ancho_pag = 792;
        largo_pag = 612;
        break;
    }

    min_y = malloc(65);
    max_y = malloc(65);
    min_x = malloc(65);
    max_x = malloc(65);

    limites('x', n, min_x, max_x);
    limites('y', n, min_y, max_y);

    GracePrintf("read nxy \"%s/%s.dat\"", dir, nom_arch);
    GracePrintf("page size %u, %u", ancho_pag, largo_pag);
    GracePrintf("world xmin %s", min_x);
    GracePrintf("world ymin 1.02*%s-.02*%s", min_y, max_y);
    if (strcmp(info_gr.etq_x, "r/\\lambda_1(0)") == 0)
        GracePrintf("world xmax .55*%s", max_x);
    else
        GracePrintf("world xmax %s", max_x);
    GracePrintf("world ymax 1.05*%s-.05*%s", max_y, min_y);
    GracePrintf("legend on");
    GracePrintf("saveall \"%s/%s.agr\"", dir, nom_arch);

    GracePrintf("print to \"%s/%s.%s\"", dir, nom_arch, form_arch);
    /* Formato de archivo de salida. */
    if (strcmp(form_arch, "png") == 0)
    {
        GracePrintf("hardcopy device \"PNG\"");
        GracePrintf("device \"PNG\" font antialiasing on");
        GracePrintf("device \"PNG\" op \"compression:9\"");
        /* Opcional: GracePrintf("device \"PNG\" op \"transparent:on\""); */
    }
    else if (strcmp(form_arch, "eps") == 0)
    {
        GracePrintf("hardcopy device \"EPS\"");
        GracePrintf("device \"EPS\" op \"level2\"");
    }
    GracePrintf("print");

    if (o_graf < 1) GraceClose();
    else GraceClosePipe();
    free(min_x);
    free(max_x);
    free(min_y);
    free(max_y);
#endif // __linux__
}
