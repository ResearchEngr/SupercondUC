#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "filerw.h"

static const char *sc_library = "sc_materials.txt";

/*
 * Display an error message according to the error code received as an argument
 */
void display_err_msg(int code, char *msg)
{
    printf("\n");
    switch(code)
    {
    case 0:
        printf("* Error de sintaxis: el material ``%s'' no existe en la biblioteca sc_materials.txt", msg);
        break;
    case 1:
        printf("* Error de sintaxis: identificador `%c' inválido", *msg);
        break;
    case 2:
        printf("* Error de sintaxis: se esperaba el signo \"=\" después del identificador `%c'", *msg);
        break;
    case 3:
        printf("* Error de sintaxis: identificador `%c' duplicado", *msg);
        break;
    case 4:
        printf("* Valor fuera de rango: %s", msg);
        break;
    case 5:
        printf("* Inconsistencia: %s", msg);
        break;
    case 6:
        printf("* Error de cálculo: %s", msg);
        break;
    default:
        printf("* %s", msg);
        break;
    }
    printf(".");
}

/*
 * Get part of the correct print format for output data files
 */
char *data_format(short num_of_data_cols, short dec_digits)
{
    char *data = malloc(3 + pow(2, dec_digits/10)),
          *format_tmp = malloc(num_of_data_cols*(4 + pow(2, dec_digits/10)) + 1); /* length(num_of_data_cols*(length("%.g ") = 4 + dec_digits) - a_blank_space_not_required_in_last_term + 2 (due to "\n")) */
    short i;

    sprintf(data, "%%.%dg", dec_digits);
    sprintf(format_tmp, "\n");
    for (i = 0; i < num_of_data_cols; i++)
    {
        strcat(format_tmp, data);
        if (i < num_of_data_cols - 1) strcat(format_tmp, " ");
    }

    return format_tmp;
}

/*
 * Create a FILE pointer and stores the heading information for the data file specified
 */
FILE* create_file(char *dir, char *file_name, char *ext, char char_com, char *material, char *plot_soft_name, char *user, short calc_opt, unsigned ni, unsigned n_iter, unsigned n_val)
{
    time_t seconds;
    struct tm *curr_date;
    char *path = malloc(strlen(dir) + strlen(file_name) + 5);
    FILE *file;

    /* Data file Heading */
    sprintf(path, "%s/%s.%s", dir, file_name, ext);
    file = fopen(path, "wt");
    fprintf(file, "%c Archivo de datos generado con SupercondUC.\n", char_com);
    if (strcmp(plot_soft_name, "") != 0) fprintf(file, "%c Formato compatible de archivo: %s.\n", char_com, plot_soft_name);
    fprintf(file, "%c Usuario: %s.\n", char_com, user);

    /* System date and time */
    seconds = time(NULL);
    curr_date = localtime(&seconds);
    fprintf(file, "%c Fecha (día/mes/año): %02d/%02d/%04d.\n", char_com, curr_date->tm_mday, curr_date->tm_mon + 1, curr_date->tm_year + 1900);
    fprintf(file, "%c Hora: %02dh %02dm %02ds.\n%c\n", char_com, curr_date->tm_hour, curr_date->tm_min, curr_date->tm_sec, char_com);

    fprintf(file, "%c Compuesto: %s.\n", char_com, material);
    fprintf(file, "%c Método numérico: Diferencias finitas / Newton-Raphson.\n", char_com);

    if (calc_opt > 0) /* For variations of $T$ or $\tilde\gamma$ */
        fprintf(file, "%c Cantidad de puntos válidos calculados: %d de %d.\n", char_com, n_val, n_iter);
    else /* For numerical solution of GL equations */
        fprintf(file, "%c Cantidad de puntos calculados: %d.\n", char_com, ni);

    fprintf(file, "%c\n%c Observaciones del usuario:\n", char_com, char_com);

    return file;
}

char *solution_data_format(char *param, short dec_digits, char char_com)
{
    char *result = malloc(4 + pow(2, dec_digits/10) + strlen(param));
    sprintf(result, "%c %s = %s\n", char_com, param, data_format(1, dec_digits) + 1);
    return result;
}

/*
 * Create a new node struct so it can be added to a linked list data structure
 */
struct node* new_node(char *msg, int lin, int col)
{
    node *node_new = malloc(sizeof(node));
    if (node_new == NULL)
    {
        printf("Error crítico: la memoria libre requerida por el programa es insuficiente.");
        return NULL;
    }
    node_new->msg = msg;
    node_new->lin = lin;
    node_new->col = col;
    node_new->next = NULL;
    return node_new;
}

/*
 * Add a node struct to a linked list data structure
 */
void add_node(node **plista, char *msg, int lin, int col)
{
    node *node_aux;
    if ((node_aux = new_node(msg, lin, col)) == NULL) exit(1);
    if (*plista == NULL) *plista = node_aux;
    else *plista = ((*plista)->next = node_aux);
}

/*
 * Get the parameter values from the superconductor materials library ``sc_materials.txt''
 */
superconductor *get_sc_parameters(char *material)
{
    FILE *file;
    char *data, *data0, *excep, ch;
    node *new_err_msg, *node_aux;
    superconductor *sc[10], *sc_aux = malloc(sizeof(superconductor));
    unsigned lin = 0, col;
    unsigned short i, j, ns = 0;
    bool is_comment, still_data_to_read, between_parentheses, reading_data, data_read, found_material;

    file = fopen(sc_library, "rt");
    if (file == NULL)
    {
        display_err_msg(7, "No se encuentra el archivo ``sc_materials.txt''");
        printf("\n\n");
        exit(1);
    }

    /* Data size is set to longest IUPAC name length (titin)*/
    data = malloc(31);
    data0 = data;

    new_err_msg = NULL;
    /* A sentinel node is added as the list head */
    add_node(&new_err_msg, NULL, 0, 0);
    node_aux = new_err_msg;
    sc_aux->temp_c = NULL;

    while (ch != EOF)
    {
        /* Line start */
        i = j = col = 0;
        lin++;
        is_comment = still_data_to_read = between_parentheses = reading_data = data_read = found_material = false;

        do
        {
            ch = fgetc(file);
            col++;
            if ((ch != ' ') && (ch != '\r') && (ch != '\n')) still_data_to_read = true;

            if ((!is_comment) && (still_data_to_read))
            {
                excep = NULL;
                switch (ch)
                {
                case ' ':
                    if (reading_data)
                    {
                        if (between_parentheses)
                        {
                            if ((i > 2) && (j < sc_aux->bands_number - 1)) excep = "Error de sintaxis: datos incompletos.";
                            else if ((i == 2) && (j < sc_aux->bands_number)) excep = "Error de sintaxis: por favor, revise la(s) temperatura(s) crítica(s) del compuesto. Recuerde definir la temperatura crítica del sistema.";
                        }
                        reading_data = false;
                    }
                    break;
                case '(':
                    if ((i > 1) && (sc_aux->bands_number < 2)) excep = "Error de sintaxis: no debe emplearse paréntesis en materiales de una sola banda.";
                    else if (i < 2) excep = "Error de sintaxis: uso incorrecto del paréntesis.";
                    else if (between_parentheses) excep = "Error de sintaxis: no se ha cerrado un paréntesis anterior.";
                    between_parentheses = true;
                    reading_data = false;
                    break;
                case ')':
                    if ((i > 1) && (sc_aux->bands_number < 2)) excep = "Error de sintaxis: no debe emplearse paréntesis en materiales de una sola banda.";
                    else if (i < 2) excep = "Error de sintaxis: uso incorrecto del paréntesis.";
                    else if (!between_parentheses) excep = "Error de sintaxis: paréntesis de cierre sin su paréntesis de apertura correspondiente.";
                    between_parentheses = false;
                    reading_data = false;
                    break;
                case ',':
                    if (!reading_data) excep = "Error de sintaxis.";
                    else if (!between_parentheses) excep = "Error de sintaxis: emplee el caracter ',' sólo para separar los datos que se encuentran entre paréntesis.";
                    else if ((i > 2) && (j == sc_aux->bands_number - 1)) excep = "Error de sintaxis: se definieron más datos de los requeridos.";
                    else if ((i == 2) && (j == sc_aux->bands_number)) excep = "Error de sintaxis: se definieron más datos de los requeridos.";
                    reading_data = false;
                    break;
                case '#':
                    is_comment = true;
                case '\n':
                case EOF:
                    if (((i > 3) && (data == data0)) || ((i < 4) && ((reading_data) || (data_read)))) excep = "Error de sintaxis: datos incompletos.";
                    reading_data = false;
                    break;
                default:
                    if (i < 5)
                    {
                        if (!reading_data)
                        {
                            reading_data = true;
                            data_read = false;
                        }
                        *data++ = ch;
                    }
                    else excep = "Error de sintaxis: se definieron más datos de los requeridos.";
                    break;
                }
                if ((!reading_data) && (!data_read) && (i < 5))
                {
                    *data = '\0';
                    data = data0;
                    switch (i)
                    {
                    case 0:
                        if (strcmp(material, data) == 0)
                        {
                            found_material = true;
                            if (ns > 0) excep = "Advertencia: se encontró más de una definición para este compuesto. Se seleccionarán los datos de la primera definición encontrada.";
                        }
                        break;
                    case 1:
                        if (sc_aux->bands_number != atoi(data))
                        {
                            sc_aux->bands_number = atoi(data);
                            /* Deallocating memory no longer needed */
                            sc_aux->temp_c = malloc((sc_aux->bands_number + 1)*sizeof(*sc_aux->temp_c)); /* One memory space is added to include temp_c */
                        }
                        break;
                    case 2:
                        sc_aux->temp_c[j] = atof(data);
                        break;
                    case 3:
                        sc_aux->w = atof(data);
                        break;
                    case 4:
                        sc_aux->gamma_a = atof(data);
                        break;
                    default:
                        printf("\n%d\n", i);
                        break;
                    }
                    if ((i == 2) && (sc_aux->bands_number > 1) && (j < sc_aux->bands_number)) j++;
                    else
                    {
                        j = 0;
                        i++;
                    }
                    data_read = true;
                }
                if (excep != NULL) add_node(&new_err_msg, excep, lin, col);
            }
        }
        while ((ch != '\n') && (ch != EOF)); /* Blank lines are discarded */
        if (found_material)
        {
            sc[ns++] = sc_aux;
            sc_aux = malloc(sizeof(superconductor));
        }
    }

    free(data);
    free(sc_aux->temp_c);
    free(sc_aux);

    /* Printing syntax error and/or warning messages about library ``sc_materials.txt'' */
    new_err_msg = node_aux;
    if (new_err_msg->next != NULL) printf("\nErrores y/o advertencias:");
    while ((new_err_msg = new_err_msg->next) != NULL)
    {
        free(node_aux);
        printf("\n(lín: %d, col: %d) -> %s", new_err_msg->lin, new_err_msg->col, new_err_msg->msg);
        node_aux = new_err_msg;
    }

    free(node_aux);
    if (ns > 0) return *sc;
    else return NULL;
}
