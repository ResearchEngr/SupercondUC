#ifndef FILERW_H_INCLUDED
#define FILERW_H_INCLUDED

typedef struct nodo nodo;
typedef struct superconductor superconductor;

struct nodo
{
    int lin, col;
    char *msje;
    nodo *sgte;
};

struct superconductor
/* Esta estructura almacena información sobre los parámetros de un material superconductor multibandas. */
{
    int n_bandas;
    float *temp_c, w, gamma_a;
};

/* Prototipos de funciones. */
void agrega(nodo**, char*, int, int);
FILE* crear_archivo(char*, char*, char*, char, char*, char*, char*, short, unsigned, unsigned, unsigned);
void error(int, char*);
char *formato_datos(short, short);
char *formato_result(char*, short, char);
nodo* nuevo_nodo(char *msje, int, int);
superconductor *obtener_parametros(char*);

#endif
