#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "filerw.h"

static const char *biblioteca = "compuestos.txt";

void error(int codigo, char *msje)
/*
 * Descripción:
 * Esta función muestra un mensaje de error que resulta de la validación de los argumentos de entrada de la aplicación.
 *
 * Parámetros:
 * codigo : número que indica un código de error.
 * msje   : cadena de caracteres que almacena un mensaje, indicando detalles acerca del error.
 */
{
    printf("\n");
    switch(codigo)
    {
    case 0:
        printf("* Error de sintaxis: el material ``%s'' no existe en la biblioteca compuestos.txt", msje);
        break;
    case 1:
        printf("* Error de sintaxis: identificador `%c' inválido", *msje);
        break;
    case 2:
        printf("* Error de sintaxis: se esperaba el signo \"=\" después del identificador `%c'", *msje);
        break;
    case 3:
        printf("* Error de sintaxis: identificador `%c' duplicado", *msje);
        break;
    case 4:
        printf("* Valor fuera de rango: %s", msje);
        break;
    case 5:
        printf("* Inconsistencia: %s", msje);
        break;
    case 6:
        printf("* Error de cálculo: %s", msje);
        break;
    default:
        printf("* %s", msje);
        break;
    }
    printf(".");
}

char *formato_datos(short n_datos, short decimales)
/*
 * Descripción:
 * Esta función retorna un puntero de tipo char, el cual apunta a una cadena de caracteres que almacena el formato adecuado de representación de los datos en los archivos de salida.
 *
 * Parámetros:
 * n_datos   : cantidad de datos contenidos en una línea del archivo.
 * decimales : cantidad de decimales mostrados en la solución.
 */
{
    char *dato = malloc(3 + pow(2, decimales/10)),
          *formato = malloc(n_datos*(4 + pow(2, decimales/10)) + 1); /* longitud(nDatos*(longitud("%.g ") = 4 + n.ro de decimales de precisión de cada dato) - un espacio en blanco que no se coloca en el último término + 2 (debido a "\n")) */
    short i;

    sprintf(dato, "%%.%dg", decimales);
    sprintf(formato, "\n");
    for (i = 0; i < n_datos; i++)
    {
        strcat(formato, dato);
        if (i < n_datos - 1) strcat(formato, " ");
    }

    return formato;
}


FILE* crear_archivo(char *dir, char *nom_arch, char *ext, char char_coment, char *material, char *nom_soft_graf, char *usuario, short o_calc, unsigned ni, unsigned n_bar, unsigned n_val)
/*
 * Descripción:
 * Esta función retorna un puntero de tipo FILE que contiene la dirección de memoria del archivo de datos que se desea crear.
 *
 * Parámetros:
 * dir           : cadena de caracteres que contiene la ruta relativa donde se creará el archivo de datos.
 * nom_arch      : cadena de caracteres que contiene el nombre del archivo de datos.
 * char_coment   : caracter utilizado como delimitador de comentarios. Por defecto, se emplea el caracter `#'.
 * material      : cadena de caracteres que contiene el nombre del material seleccionado por el usuario.
 * nom_soft_graf : cadena de caracteres que contiene el nombre del software de representación gráfica. Por defecto, Grace.
 * usuario       : nombre del usuario del sistema.
 * o_calc        : opción que define el tipo de cálculo a realizar.
 * ni            : cantidad de puntos internos de la retícula, que constituyen las soluciones de GL.
 * n_bar         : cantidad de puntos definidos en el rango de variación de $T$ o de $\gt$. Solo si o_calc$>0$.
 * n_val         : cantidad de puntos solución válidos.
 */
{
    time_t segundos;
    struct tm *tiempo;
    char *ruta = malloc(strlen(dir) + strlen(nom_arch) + 5);
    FILE *archivo;

    /* Encabezado del archivo de datos. */
    sprintf(ruta, "%s/%s.%s", dir, nom_arch, ext);
    archivo = fopen(ruta, "wt");
    fprintf(archivo, "%c Archivo de datos generado con SupercondUC.\n", char_coment);
    if (strcmp(nom_soft_graf, "") != 0) fprintf(archivo, "%c Formato compatible de archivo: %s.\n", char_coment, nom_soft_graf);
    fprintf(archivo, "%c Usuario: %s.\n", char_coment, usuario);

    /* Fecha y hora del sistema. */
    segundos = time(NULL);
    tiempo = localtime(&segundos);
    fprintf(archivo, "%c Fecha (día/mes/año): %02d/%02d/%04d.\n", char_coment, tiempo->tm_mday, tiempo->tm_mon + 1, tiempo->tm_year + 1900);
    fprintf(archivo, "%c Hora: %02dh %02dm %02ds.\n%c\n", char_coment, tiempo->tm_hour, tiempo->tm_min, tiempo->tm_sec, char_coment);

    fprintf(archivo, "%c Compuesto: %s.\n", char_coment, material);
    fprintf(archivo, "%c Método numérico: Diferencias finitas / Newton-Raphson.\n", char_coment);

    if (o_calc > 0) /* Para variaciones de $T$ o $\tilde\gamma$. */
        fprintf(archivo, "%c Cantidad de puntos válidos calculados: %d de %d.\n", char_coment, n_val, n_bar);
    else /* Para soluciones numéricas de las ecuaciones de GL. */
        fprintf(archivo, "%c Cantidad de puntos calculados: %d.\n", char_coment, ni);

    fprintf(archivo, "%c\n%c Observaciones del usuario:\n", char_coment, char_coment);

    return archivo;
}

char *formato_result(char *parametro, short decimales, char char_coment)
/*
 * Descripción:
 * Esta función retorna un puntero de tipo char, que apunta a una cadena de caracteres que almacena el formato adecuado para la representación de los parámetros $f_{1_\infty}$, $f_{2_\infty}$, $L_{h_1}$, $L_{h_2}$, $\lambda$ y $B_0$, en el archivo de datos que contiene las soluciones de GL.
 *
 * Parámetros:
 * parametro   : código en \TeX del parámetro físico.
 * decimales   : cantidad de decimales mostrados en la solución.
 * char_coment : caracter utilizado como delimitador de comentarios. Por defecto, se emplea el caracter `#'.
 */
{
    char *result = malloc(4 + pow(2, decimales/10) + strlen(parametro));
    sprintf(result, "%c %s = %s\n", char_coment, parametro, formato_datos(1, decimales) + 1);
    return result;
}

struct nodo* nuevo_nodo(char *msje, int lin, int col)
/*
 * Descripción:
 * Esta función retorna una estructura compleja de datos de tipo nodo, utilizada para crear un nodo nuevo que se agregue a una lista enlazada simple.
 *
 * Parámetros:
 * msje : cadena de caracteres que contiene un mensaje que será almacenado en el nodo.
 * lin  : número de línea del archivo <<compuestos.txt>> donde se presenta alguna excepción en la lectura del archivo.
 * col  : número de columna del archivo <<compuestos.txt>> donde se presenta alguna excepción en la lectura del archivo.
 */
{
    nodo *nodo_nuevo = malloc(sizeof(nodo));
    if (nodo_nuevo == NULL)
    {
        printf("Error crítico: la memoria libre requerida por el programa es insuficiente.");
        return NULL;
    }
    nodo_nuevo->msje = msje;
    nodo_nuevo->lin = lin;
    nodo_nuevo->col = col;
    nodo_nuevo->sgte = NULL;
    return nodo_nuevo;
}

void agrega(nodo **plista, char *msje, int lin, int col)
/*
 * Descripción:
 * Esta función tiene como propósito añadir un nodo a una lista enlazada simple.
 *
 * Parámetros:
 * plista : arreglo de punteros que contiene las direcciones de memoria de los nodos que componen la lista enlazada.
 * msje   : cadena de caracteres que contiene un mensaje que será almacenado en el nodo nuevo.
 * lin    : número de línea del archivo <<compuestos.txt>> donde se presenta alguna excepción en la lectura del archivo.
 * col    : número de columna del archivo <<compuestos.txt>> donde se presenta alguna excepción en la lectura del archivo.
 */
{
    nodo *nodoAux;
    if ((nodoAux = nuevo_nodo(msje, lin, col)) == NULL) exit(1);
    if (*plista == NULL) *plista = nodoAux;
    else *plista = ((*plista)->sgte = nodoAux);
}

superconductor *obtener_parametros(char *material)
/*
 * Descripción:
 * Esta función recupera del archivo <<compuestos.txt>>, información sobre el material superconductor seleccionado por el usuario. Esta información se almacena en una estructura compleja de datos de tipo <<superconductor>> y es retornada por la función.
 *
 * Parámetro:
 * material : cadena de caracteres que contiene el nombre del material seleccionado por el usuario.
 */
{
    FILE *archivo;
    char *dato, *iniDato, *excep, caracter;
    nodo *nuevaExc, *nodoAux;
    superconductor *sc[10], *scAux = malloc(sizeof(superconductor));
    unsigned lin = 0, col;
    unsigned short i, j, ns = 0;
    bool esComentario, hayDatos, estaEntreParentesis, estaLeyendoUnDato, leyoElDato, encontroElMaterial;

    archivo = fopen(biblioteca,"rt");
    if (archivo == NULL)
    {
        error(7, "No se encuentra el archivo ``compuestos.txt''");
        printf("\n\n");
        exit(1);
    }

    /* Empleando como referencia la fórmula química IUPAC más larga hasta el momento (titina). */
    dato = malloc(31);
    iniDato = dato;

    nuevaExc = NULL;
    /* Se agrega un nodo centinela como cabeza de la lista. */
    agrega(&nuevaExc, NULL, 0, 0);
    nodoAux = nuevaExc;
    scAux->temp_c = NULL;

    while (caracter != EOF)
    {
        /* Comienzo de línea. */
        i = j = col = 0;
        lin++;
        esComentario = hayDatos = estaEntreParentesis = estaLeyendoUnDato = leyoElDato = encontroElMaterial = false;

        do
        {
            caracter = fgetc(archivo);
            col++;
            if ((caracter != ' ') && (caracter != '\r') && (caracter != '\n')) hayDatos = true;

            if ((!esComentario) && (hayDatos))
            {
                excep = NULL;
                switch (caracter)
                {
                case ' ':
                    if (estaLeyendoUnDato)
                    {
                        if (estaEntreParentesis)
                        {
                            if ((i > 2) && (j < scAux->n_bandas - 1)) excep = "Error de sintaxis: datos incompletos.";
                            else if ((i == 2) && (j < scAux->n_bandas)) excep = "Error de sintaxis: por favor, revise la(s) temperatura(s) crítica(s) del compuesto. Recuerde definir la temperatura crítica del sistema.";
                        }
                        estaLeyendoUnDato = false;
                    }
                    break;
                case '(':
                    if ((i > 1) && (scAux->n_bandas < 2)) excep = "Error de sintaxis: no debe emplearse paréntesis en materiales de una sola banda.";
                    else if (i < 2) excep = "Error de sintaxis: uso incorrecto del paréntesis.";
                    else if (estaEntreParentesis) excep = "Error de sintaxis: no se ha cerrado un paréntesis anterior.";
                    estaEntreParentesis = true;
                    estaLeyendoUnDato = false;
                    break;
                case ')':
                    if ((i > 1) && (scAux->n_bandas < 2)) excep = "Error de sintaxis: no debe emplearse paréntesis en materiales de una sola banda.";
                    else if (i < 2) excep = "Error de sintaxis: uso incorrecto del paréntesis.";
                    else if (!estaEntreParentesis) excep = "Error de sintaxis: paréntesis de cierre sin su paréntesis de apertura correspondiente.";
                    estaEntreParentesis = false;
                    estaLeyendoUnDato = false;
                    break;
                case ',':
                    if (!estaLeyendoUnDato) excep = "Error de sintaxis.";
                    else if (!estaEntreParentesis) excep = "Error de sintaxis: emplee el caracter ',' sólo para separar los datos que se encuentran entre paréntesis.";
                    else if ((i > 2) && (j == scAux->n_bandas - 1)) excep = "Error de sintaxis: se definieron más datos de los requeridos.";
                    else if ((i == 2) && (j == scAux->n_bandas)) excep = "Error de sintaxis: se definieron más datos de los requeridos.";
                    estaLeyendoUnDato = false;
                    break;
                case '#':
                    esComentario = true;
                case '\n':
                case EOF:
                    if (((i > 3) && (dato == iniDato)) || ((i < 4) && ((estaLeyendoUnDato) || (leyoElDato)))) excep = "Error de sintaxis: datos incompletos.";
                    estaLeyendoUnDato = false;
                    break;
                default:
                    if (i < 5)
                    {
                        if (!estaLeyendoUnDato)
                        {
                            estaLeyendoUnDato = true;
                            leyoElDato = false;
                        }
                        *dato++ = caracter;
                    }
                    else excep = "Error de sintaxis: se definieron más datos de los requeridos.";
                    break;
                }
                if ((!estaLeyendoUnDato) && (!leyoElDato) && (i < 5))
                {
                    *dato = '\0';
                    dato = iniDato;
                    switch (i)
                    {
                    case 0:
                        if (strcmp(material, dato) == 0)
                        {
                            encontroElMaterial = true;
                            if (ns > 0) excep = "Advertencia: se encontró más de una definición para este compuesto. Se seleccionarán los datos de la primera definición encontrada.";
                        }
                        break;
                    case 1:
                        if (scAux->n_bandas != atoi(dato))
                        {
                            scAux->n_bandas = atoi(dato);
                            /* Se libera memoria reservada con anterioridad para los parámetros del struct. */
                            scAux->temp_c = malloc((scAux->n_bandas + 1)*sizeof(*scAux->temp_c)); /* Se suma 1 espacio de memoria para incluir temp_c. */
                        }
                        break;
                    case 2:
                        scAux->temp_c[j] = atof(dato);
                        break;
                    case 3:
                        scAux->w = atof(dato);
                        break;
                    case 4:
                        scAux->gamma_a = atof(dato);
                        break;
                    default:
                        printf("\n%d\n", i);
                        break;
                    }
                    if ((i == 2) && (scAux->n_bandas > 1) && (j < scAux->n_bandas)) j++;
                    else
                    {
                        j = 0;
                        i++;
                    }
                    leyoElDato = true;
                }
                if (excep != NULL) agrega(&nuevaExc, excep, lin, col);
            }
        }
        while ((caracter != '\n') && (caracter != EOF)); /* Descartando líneas en blanco. */
        if (encontroElMaterial)
        {
            sc[ns++] = scAux;
            scAux = malloc(sizeof(superconductor));
        }
    }

    free(dato);
    free(scAux->temp_c);
    free(scAux);

    /* Impresión de mensajes de error y/o advertencias respecto al archivo componentes.txt. */
    nuevaExc = nodoAux;
    if (nuevaExc->sgte != NULL) printf("\nErrores y/o advertencias:");
    while ((nuevaExc = nuevaExc->sgte) != NULL)
    {
        free(nodoAux);
        printf("\n(lín: %d, col: %d) -> %s", nuevaExc->lin, nuevaExc->col, nuevaExc->msje);
        nodoAux = nuevaExc;
    }

    free(nodoAux);
    if (ns > 0) return *sc;
    else return NULL;
}
