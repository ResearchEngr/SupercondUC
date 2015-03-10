#ifndef PIPES_H_INCLUDED
#define PIPES_H_INCLUDED

typedef struct info_grafico info_grafico;

struct info_grafico
{
    char *titulo, *subtitulo, *etq_x, *etq_y, **leyenda;
    short n;
};

/* Prototipos de funciones. */
void crear_grafico(info_grafico, char*, char*, short, short, short, char*, unsigned);
void limites(char, short, char*, char*);
char* tipo(char*);

#endif
