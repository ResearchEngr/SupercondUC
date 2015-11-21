#ifndef FILERW_H_INCLUDED
#define FILERW_H_INCLUDED

typedef struct node node;
typedef struct superconductor superconductor;

struct node
{
    int lin, col;
    char *msg;
    node *next;
};

struct superconductor
/* This struct stores information about the multi-band superconductor parameters. */
{
    int bands_number;
    float *temp_c, w, gamma_a;
};

/* Function prototypes */
void add_node(node**, char*, int, int);
FILE* create_file(char*, char*, char*, char, char*, char*, char*, short, unsigned, unsigned, unsigned);
void display_err_msg(int, char*);
char *data_format(short, short);
char *solution_data_format(char*, short, char);
node* new_node(char *msg, int, int);
superconductor *get_sc_parameters(char*);

#endif
