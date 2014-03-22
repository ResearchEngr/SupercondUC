#ifndef SUPERC_H_INCLUDED
#define SUPERC_H_INCLUDED

/* Inclusión de archivos de cabecera de la biblioteca estándar de C. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>

/* Definificiones de tipos. */
typedef struct nodo nodo;
typedef struct superconductor superconductor;
//typedef struct error error;

double raizRealPos(double a, double b, double c, double d, double e)
{
    double complex a1, a2, a3, a4, a5, b1, b2, b3, b4, b5, *raiz0, *raiz;
    double A = b/a, B = c/a, C = d/a, D = e/a, /* Dividiendo por el coeficiente que multiplica al término de grado cuatro. */
           raizR, epsilon;
    unsigned short i;

    A = b/a;
    B = c/a;
    C = d/a;
    D = e/a;

    //printf("\na = %.17g, b = %.17g, c = %.17g, d = %.17g., e = %.17g.", a, b, c, d, e);
    //printf("\nA = %.17g, B = %.17g, C = %.17g, D = %.17g.\n", A, B, C, D);

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

    raiz0 = (double complex *) malloc(4*sizeof(double complex));
    raiz = raiz0;
    *raiz = a1 - b3 - .5*csqrt(b4 - b5);
    *++raiz = a1 - b3 + .5*csqrt(b4 - b5);
    *++raiz = a1 + b3 - .5*csqrt(b4 + b5);
    *++raiz = a1 + b3 + .5*csqrt(b4 + b5);

    epsilon = 1e-6;
    raiz = raiz0;
    for(i = 0; i < 4; i++)
    {
        raizR = creal(*raiz);
        if ((raizR > 0) && (fabs(cimag(*raiz)) < epsilon)) break;
        raiz++;
    }

    printf("\n%.17g + j(%.17g)\n%.17g + j(%.17g)\n%.17g + j(%.17g)\n%.17g + j(%.17g)\n", *raiz0, *(raiz0 + 1), *(raiz0 + 2), *(raiz0 + 3));
    printf("\nSe selecciona: %.17g\n", raizR);

    free(raiz0);
    return raizR;
}

void mostrarAyuda(char* nomArch)
{
    //¡! Esta función debería estar definida en un archivo de cabecera aparte y leer un archivo de texto que contenga la ayuda del programa completo.
    printf("\nSintaxis: param m=<Material> T=<Temperatura_inicial> n=<Cantidad_de_puntos> ...");
    printf("\n\nMaterial: fórmula química del material. Debe estar contenido en la biblioteca compuestos.txt.");
    printf("\n\nTemperatura_inicial: temperatura a partir de la cual se comienzan a realizar los cálculos.");
    printf("\n\nCantidad_de_puntos: cantidad de puntos de temperatura para los cuales se realizan los cálculos pertinentes.");
    printf("\n\n");
    exit(1);
}

void error(int codigo, char* cadena)
{
    printf("\n");
    switch(codigo)
    {
    case 0:
        printf("* Error de sintaxis: el material %s no existe en la biblioteca componentes.txt.\n", cadena);
        break;
    case 1:
        printf("* Error de sintaxis: identificador \"%c\" inválido.", *cadena);
        break;
    case 2:
        printf("* Error de sintaxis: se esperaba el signo \"=\" después del identificador \"%c\".", *cadena);
        break;
    case 3:
        printf("* Error de sintaxis: identificador \"%c\" duplicado.", *cadena);
        break;
    case 4:
        printf("* Valor fuera de rango: %s", cadena);
        break;
    case 5:
        printf("* Inconsistencia: %s", cadena);
        break;
    case 6:
        printf("* Error de cálculo: %s", cadena);
        break;
    case 7:
        printf("* %s", cadena);
        break;
    }
    printf("\n");
}

char *biblioteca = "compuestos.txt";

struct superconductor
{
    /* Estructura utilizada para almacenar los valores de un superconductor de una o más bandas. */
    int nBandas; /* Define el número de bandas a considerar para un compuesto en particular. */
    float *Tc, *lambda0, *xi0, *m, *Hc, *Vf, gammaA;
};

struct nodo
{
    int lin, col;
    char *msje;
    nodo *sgte;
};

struct nodo* NuevoNodo(char *msje, int lin, int col)
{
    nodo *nodoNuevo = (nodo *) malloc(sizeof(nodo));
    if (nodoNuevo == NULL)
    {
        printf("Error crítico: la memoria libre requerida por el programa es insuficiente.");
        return NULL;
    }
    nodoNuevo->msje = msje;
    nodoNuevo->lin = lin;
    nodoNuevo->col = col;
    nodoNuevo->sgte = NULL;
    return nodoNuevo;
};

void Agrega(nodo **plista, char *msje, int lin, int col)
{
    nodo *nodoAux;
    if ((nodoAux = NuevoNodo(msje, lin, col)) == NULL) exit(1);
    if (*plista == NULL) *plista = nodoAux;
    else *plista = ((*plista)->sgte = nodoAux);
};

//struct nodo
//{
//    int lin, col;
//    char *msje;
//    struct nodo *sgte;
//};
//
//struct nodo* NuevoNodo(char *msje, int lin, int col)
//{
//    struct nodo *nodoNuevo = (struct nodo *) malloc(sizeof(struct nodo));
//    if (nodoNuevo == NULL)
//    {
//        printf("Error crítico: la memoria libre requerida por el programa es insuficiente.");
//        return NULL;
//    }
//    nodoNuevo->msje = msje;
//    nodoNuevo->lin = lin;
//    nodoNuevo->col = col;
//    nodoNuevo->sgte = NULL;
//    return nodoNuevo;
//};
//
//void Agrega(struct nodo **plista, char *msje, int lin, int col)
//{
//    struct nodo *nodoAux;
//    if ((nodoAux = NuevoNodo(msje, lin, col)) == NULL) exit(1);
//    if (*plista == NULL) *plista = nodoAux;
//    else *plista = ((*plista)->sgte = nodoAux);
//};

superconductor* obtenerParametros(char *material)
/*
 * Función que devuelve un dato de tipo superconductor de acuerdo a lo indicado en
 * el parámetro material y las propiedades físicas inherentes al material.
 */
{
    FILE* archivo;
    char *dato, *iniDato, *excep, caracter;
    nodo *nuevaExc, *nodoAux;
    superconductor *sc[10], *scAux = (superconductor *) malloc(sizeof(superconductor));
    unsigned lin = 0, col;
    unsigned short i, j, ns = 0;
    bool esComentario, hayDatos, estaEntreParentesis, estaLeyendoUnDato, leyoElDato, encontroElMaterial;

    archivo = fopen(biblioteca,"rt");
    if (archivo == NULL)
    {
        /* Se notifica al usuario que el archivo se perdió o fue borrado y por lo tanto se creará uno predeterminado */
        // Sería necesario verificar que el archivo compuestos~.txt tampoco existe, con el propósito de que se le conceda
        // al usuario la opción de utilizarlo como archivo predeterminado, en caso de que el original haya sido borrado.
        // Este archivo sería una copia de seguridad que se crearía al inicio de la aplicación o cuando se haya cargado co-
        // rrectamente una información al archivo, por lo cual, sería necesario validar la integridad del archivo original
        // antes de realizar una copia.
        printf("\nEl archivo \"compuestos.txt\" no existe y se creará un archivo con la base de datos predeterminada. ");
        printf("Consulte su contenido para más detalles sobre los materiales que puede utilizar.\n");

        /* Se crea un archivo predeterminado en caso de que no exista un archivo de copia de seguridad */
        archivo = fopen(biblioteca,"wt");
        fprintf(archivo,"Este archivo comprende una biblioteca de materiales que permite al usuario añadir nuevos compuestos o editar\n");
        fprintf(archivo,"# información existente. Cada línea del archivo corresponde a un material diferente y el orden en el cual se\n");
        fprintf(archivo,"# definen sus parámetros es como sigue:\n#\n");
        fprintf(archivo,"# - Fórmula química (material):\n");
        fprintf(archivo,"#       En este campo se escribe la fórmula química del compuesto.\n#\n");
        fprintf(archivo,"# - Número de bandas (nBandas):\n");
        fprintf(archivo,"#       Un número entero mayor que cero. Indica la cantidad de bandas a considerar para el material.\n#\n");

        archivo = fopen(biblioteca,"rt");
    }


    // Verificando la integridad del archivo compuestos.txt.
    // Es necesario tener un archivo réplica del original que sirva como modelo y respaldo en caso de que el
    // usuario borre intencional o accidentalmente las instrucciones referidas al formato. Este archivo
    // ("compuestos~.txt") debería crearse automáticamente al ejecutar la aplicación, en caso de que no exista.

    {
        // ...
        // ...
        // ...
    }





    dato = (char *) malloc (31); // Empleando como referencia la fórmula química IUPAC más larga hasta el momento (titina).
    iniDato = dato;

    nuevaExc = NULL;
    /* Se agrega un nodo centinela como cabeza de la lista. */
    Agrega(&nuevaExc, NULL, 0, 0);
    nodoAux = nuevaExc;
    scAux->Tc = scAux->lambda0 = scAux->xi0 = scAux->m = scAux->Hc = scAux->Vf = NULL;

    while (caracter != EOF)
    {
        /* Comienzo de línea. */
        i = j = col = 0;
        lin++;
        esComentario = hayDatos = estaEntreParentesis = estaLeyendoUnDato = leyoElDato = encontroElMaterial = false;

        do
        {
            caracter = fgetc(archivo);
//            printf("%c",caracter);
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
                            if ((i > 2) && (j < scAux->nBandas - 1)) excep = "Error de sintaxis: datos incompletos.";
                            else if ((i == 2) && (j < scAux->nBandas)) excep = "Error de sintaxis: por favor, revise la(s) temperatura(s) crítica(s) del compuesto. Recuerde definir la temperatura crítica del sistema.";
                        estaLeyendoUnDato = false;
                    }
                    break;
                case '(':
                    if ((i > 1) && (scAux->nBandas < 2)) excep = "Error de sintaxis: no debe emplearse paréntesis en materiales de una sola banda.";
                    else if (i < 2) excep = "Error de sintaxis: uso incorrecto del paréntesis.";
                    else if (estaEntreParentesis) excep = "Error de sintaxis: no se ha cerrado un paréntesis anterior.";
                    estaEntreParentesis = true;
                    estaLeyendoUnDato = false;
                    break;
                case ')':
                    if ((i > 1) && (scAux->nBandas < 2)) excep = "Error de sintaxis: no debe emplearse paréntesis en materiales de una sola banda.";
                    else if (i < 2) excep = "Error de sintaxis: uso incorrecto del paréntesis.";
                    else if (!estaEntreParentesis) excep = "Error de sintaxis: paréntesis de cierre sin su paréntesis de apertura correspondiente.";
                    estaEntreParentesis = false;
                    estaLeyendoUnDato = false;
                    break;
                case ',':
                    if (!estaLeyendoUnDato) excep = "Error de sintaxis.";
                    else if (!estaEntreParentesis) excep = "Error de sintaxis: emplee el caracter ',' sólo para separar los datos que se encuentran entre paréntesis.";
                    else if ((i > 2) && (j == scAux->nBandas - 1)) excep = "Error de sintaxis: se definieron más datos de los requeridos.";
                    else if ((i == 2) && (j == scAux->nBandas)) excep = "Error de sintaxis: se definieron más datos de los requeridos.";
                    estaLeyendoUnDato = false;
                    break;
                case '#':
                    esComentario = true;
                case '\n':
                case EOF:
                    if ((i == 8) && (dato == iniDato) || (i < 8) && ((estaLeyendoUnDato) || (leyoElDato))) excep = "Error de sintaxis: datos incompletos.";
                    estaLeyendoUnDato = false;
                    break;
                default:
                    if (i < 9)
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
                if ((!estaLeyendoUnDato) && (!leyoElDato) && (i < 9))
                {
                    *dato = '\0';
                    dato = iniDato;
                    switch (i)
                    {
                    case 0:
                        if (strcmp(material, dato) == 0)
                        {
                            encontroElMaterial = true;
                            if (ns > 0) excep = "Advertencia: se encontró una definición adicional para este compuesto.";
                        }
                        break;
                    case 1:
                        if (scAux->nBandas != atoi(dato))
                        {
                            scAux->nBandas = atoi(dato);
                            /* Se libera memoria reservada con anterioridad para los parámetros del struct. */

                            scAux->Tc = (float *) malloc((scAux->nBandas + 1) * sizeof(float)); /* Se suma 1 espacio de memoria para incluir Tc. */
                            scAux->lambda0 = (float *) malloc(scAux->nBandas * sizeof(float));
                            scAux->xi0 = (float *) malloc(scAux->nBandas * sizeof(float));
                            scAux->m = (float *) malloc(scAux->nBandas * sizeof(float));
                            scAux->Hc = (float *) malloc(scAux->nBandas * sizeof(float));
                            scAux->Vf = (float *) malloc(scAux->nBandas * sizeof(float));

//                            scAux->Tc = (float *) realloc(scAux->Tc, (scAux->nBandas + 1) * sizeof(float)); /* Se suma 1 espacio de memoria para incluir Tc. */
//                            scAux->lambda0 = (float *) realloc(scAux->lambda0, scAux->nBandas * sizeof(float));
//                            scAux->xi0 = (float *) realloc(scAux->xi0, scAux->nBandas * sizeof(float));
//                            scAux->m = (float *) realloc(scAux->m, scAux->nBandas * sizeof(float));
//                            scAux->Hc = (float *) realloc(scAux->Hc, scAux->nBandas * sizeof(float));
//                            scAux->Vf = (float *) realloc(scAux->Vf, scAux->nBandas * sizeof(float));
                        }
                        break;
                    case 2:
                        scAux->Tc[j] = atof(dato);
                        break;
                    case 3:
                        scAux->lambda0[j] = atof(dato);
                        break;
                    case 4:
                        scAux->xi0[j] = atof(dato);
                        break;
                    case 5:
                        scAux->m[j] = atof(dato);
                        break;
                    case 6:
                        scAux->Hc[j] = atof(dato);
                        break;
                    case 7:
                        scAux->Vf[j] = atof(dato);
                        break;
                    case 8:
                        scAux->gammaA = atof(dato);
                        break;
                    default:
                        printf("\n%d\n", i);
                        break;
                    }
                    if ((i > 1) && (scAux->nBandas > 1) && (i < 8))
                        if ((i == 2) && (j < scAux->nBandas) || (j < scAux->nBandas -1)) j++;
                        else
                        {
                            j = 0;
                            i++;
                        }
                    else i++;
                    leyoElDato = true;
                }
                if (excep != NULL) Agrega(&nuevaExc, excep, lin, col);
            }
        }
        while ((caracter != '\n') && (caracter != EOF)); /* Descartando líneas en blanco. */
        if (encontroElMaterial)
        {
            sc[ns++] = scAux;
            scAux = (superconductor *) malloc(sizeof(superconductor));
        }
    }

    free(dato);
    free(scAux->Tc);
    free(scAux->lambda0);
    free(scAux->xi0);
    free(scAux->m);
    free(scAux->Hc);
    free(scAux->Vf);
    free(scAux);
    if (ns < 1) Agrega(&nuevaExc, "Error de sintaxis: el material no existe en la biblioteca componentes.txt.\n", 0, 0);

    /** Impresión de mensajes de error y/o advertencias respecto al archivo componentes.txt. */
    printf("\nErrores y/o advertencias:");
    nuevaExc = nodoAux;
    if (nuevaExc->sgte == NULL) printf(" no se encontraron.\n");
    while ((nuevaExc = nuevaExc->sgte) != NULL)
    {
        free(nodoAux);
        printf("\n(lín: %d, col: %d) -> %s", nuevaExc->lin, nuevaExc->col, nuevaExc->msje);
        nodoAux = nuevaExc;
    }
    free(nodoAux);

    if (ns > 0) return sc[0];
    else return NULL;
}

#endif // SUPERC_H_INCLUDED
