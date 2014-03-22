#ifndef PIPES_H_INCLUDED
#define PIPES_H_INCLUDED

#include <grace_np.h>
#define GNUPLOT "gnuplot -persist"

void my_error_function(const char *msg)
{
    fprintf(stderr, "library message: \"%s\"\n", msg);
}

char* tipografia(char *material, char formArch[2])
{
    char *t;

    //Grace
    if (strcmp(formArch, "Gr") == 0)
    {
        if (strcmp(material, "MgB2") == 0)
        {
            t = "MgB\\s2\\N";
        }
        else if (strcmp(material, "OsB2") == 0)
        {
            t = "OsB\\s2\\N";
        }
        else if (strcmp(material, "FeSe1-x") == 0)
        {
            t = "FeSe\\s1-x\\N";
        }
        else if (strcmp(material, "V3Si") == 0)
        {
            t = "V\\s3\\NSi";
        }
        else
        {
            t = material;
        }
    }
    //gnuplot
    else if (strcmp(formArch, "gp") == 0)
    {
        if (strcmp(material, "MgB2") == 0)
        {
            t = "MgB_2";
        }
        else if (strcmp(material, "OsB2") == 0)
        {
            t = "OsB_2";
        }
        else if (strcmp(material, "FeSe1-x") == 0)
        {
            t = "FeSe_{1-x}";
        }
        else if (strcmp(material, "V3Si") == 0)
        {
            t = "V_3Si";
        }
        else
        {
            t = material;
        }
    }
    return t;
}

#endif // PIPES_H_INCLUDED
