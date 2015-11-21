#ifndef PIPES_H_INCLUDED
#define PIPES_H_INCLUDED

typedef struct plot_info plot_info;

struct plot_info
{
    char *title, *subtitle, *lbl_x, *lbl_y, **key;
    short n;
};

/* Function prototype */
void create_plot(plot_info, char*, char*, short, short, short, char*, unsigned);
void plot_bounds(char, short, char*, char*);
char* typo_convert(char*);

#endif
