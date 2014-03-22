#include <complex.h>
#include "superc.h"
#include "pipes.h"

int main(int argc, char *argv[])
{
    double complex *eta0, *eta, a1, a2, a3, b1, b2, b3, b4, b5, b6;
    double temp, temp2, tempF, epsilon, prctjErr,
           r, rn, factor, rf, h, auxh,
           f1, f1_2, f1n, d_f1, d_f1n, d2_f1, d_f1a,
           f2, f2_2, f2n, d_f2, d_f2n, d2_f2, d_f2a,
           f1b,f2b,ab,f10,f20,a0,m,f1a,f2a,aa,aux,
           A, B, C, B01, B02, B0,
           a, a_1r, na1r_2, a1_r, an, d_a, d_aa, d_an, d2_a,
           etaR, c1, c2, c3, c4, c5, c6, c7, c8, c9,
           h10, h11, h12, h13, h20, h21, h22, h23,
           w, w1, w2, w3, w4, k1[6], k2[6], k3[6], k4[6],
           alfa, beta, mi, gamma, gammaT, kappa1, nBandas, k, aux1, aux2, aux3, menorLambda0, *xi0, *lambda, xi_ef, lambda_ef, b;
    int s, n, i, j, ni, auxni, precision, codErr[argc - 1], nErr, escala;
    char *material, *nomArch, *archSc, charComent, formArch[2], *nomSoftGraf, *usuario, metodoNum[5], *nomMetodoNum, *comentarios, *valor, formato[50], ident[argc - 1], *msjErr[argc - 1], *comando;
    time_t segundos;
    struct tm *tiempo;
    superconductor* sc;
    FILE *archivo, *igs;


    /**<  Si el usuario introduce el nombre del programa por consola y presiona Enter (sin argumentos), se activa un menú que guíe
    /**<  en la selección de los argumentos de entrada para realizar el cálculo. */
    switch(argc)
    {
        {
        /* Inicio de bloque - Interpretación de comandos. */
        case 1:
            /* Ejecutar asistente para guiar al usuario. De momento sólo mostramos la ayuda al usuario si se cumple esta condición. */
            mostrarAyuda("vortex.hlp");
            break;

        default:

            /* Definición de valores predeterminados. */
            temp = 0;
            escala = 0;
            factor = 4;
            ni = 1000;
            factor = 4; /* Factor por el que se multiplica la menor longitud de coherencia de las bandas para así hallar el valor de rf. */
            precision = 15;
            strcpy(formArch, "Gr");
            strcpy(metodoNum, "rk4");
            nomArch = "vortex.dat";
            usuario = "Anónimo";
            comentarios = "algoritmo en fase de pruebas";

            /* Verificación de la sintaxis y validación de los argumentos introducidos por el usuario. */
            //      Todavía queda un detalle por arreglar. Cuando un identificador aparece definido más de dos veces, el mensaje de identificador duplicado se muestra
            //      más de una vez. Sería correcto que sólo apareciera una vez, siguiendo la misma convención que se empleó para otros mensajes de error.

            /* Inicializa la cadena de caracteres "ident" en un valor nulo. De lo contrario, la cadena contendrá basura y generará errores. */
            strncpy(ident, "", argc);
            nErr = 0;

            for(i = 1; i < argc; i++)
            {
                if (*argv[i] == '?') mostrarAyuda("vortex.hlp");

                /* Se verifica que el primer caracter sea un identificador válido.
                 * m = material.
                 * M = método numérico.
                 * n = cantidad de puntos.
                 * h = paso.
                 * T = temperatura.
                 * e = especifica si T se escala a Tc.
                 * p = cantidad de decimales.
                 * a = nombre del archivo de salida.
                 * f = formato de archivo.
                 * x = factor.
                 * u = usuario.
                 * c = comentario.
                 */
                if (strchr("mMnhTexpafuc", *argv[i]) == NULL)
                {
                    /* Verifica si el caracter inválido ya había aparecido, para no repetir el mensaje de error correspondiente. */
                    if (strchr(ident, *argv[i]) == NULL)
                    {
                        codErr[nErr] = 1;
                        msjErr[nErr] = argv[i];
                        nErr++;
                    }
                    ident[i - 1] = *argv[i];
                    /* Avanza a la siguiente iteración para no mostrar ningún otro mensaje de error relativo a un identificador inválido. */
                    continue;
                }
                /* Se verifica que el segundo caracter sea el signo igual. */
                if (*(argv[i] + 1) != '=')
                {
                    codErr[nErr] = 2;
                    msjErr[nErr] = argv[i];
                    nErr++;
                }
                /* Se verifica que no existan identificadores duplicados. */
                if (strchr(ident, *argv[i]) != NULL)
                {
                    codErr[nErr] = 3;
                    msjErr[nErr] = argv[i];
                    nErr++;
                }

                valor = argv[i] + 2;

                switch(*argv[i])
                {
                case 'm':
                    if (strchr(ident, 'm') == NULL) material = valor;
                    break;


                case 'n':

                    ni = atoi(valor);
                    if ((ni < 10) || (ni > 1e9))
                    {
                        codErr[nErr] = 4;
                        msjErr[nErr] = "la cantidad de puntos debe ser un número entero positivo: 10 ≤ n ≤ 10^9.";
                        nErr++;
                    }
                    break;


                case 'x':

                    factor = atof(valor);
                    if ((factor < 1e-6) || (h > 100))
                    {
                        codErr[nErr] = 4;
                        msjErr[nErr] = "el múltiplo \"x\" de la longitud de coherencia debe ser un número real positivo: 10^-6 ≤ x ≤ 1000.";
                        nErr++;
                    }
                    break;


                case 'h':

                    h = atof(valor);
                    if ((h < 1e-6) || (h > 1))
                    {
                        codErr[nErr] = 4;
                        msjErr[nErr] = "el paso \"h\" debe ser un número real positivo: 10^-6 ≤ h ≤ 1.";
                        nErr++;
                    }
                    break;


                case 'T':

                    temp = atof(valor);
                    break;


                case 'e':

                    escala = atoi(valor);
                    if ((escala < 0) || (escala > 1))
                    {
                        codErr[nErr] = 4;
                        msjErr[nErr] = "la escala \"e\" sólo puede ser 0 o 1.";
                        nErr++;
                    }
                    break;


                case 'p':

                    precision = atoi(valor);

                    if ((precision < 1) || (precision > 20))
                    {
                        codErr[nErr] = 4;
                        msjErr[nErr] = "la precisión \"p\" debe ser un número entero positivo: 1 ≤ p ≤ 20.";
                        nErr++;
                    }
                    break;


                case 'a':

                    nomArch = valor;
                    break;


                case 'M':

                    // Es necesario convertir la cadena de caracteres apuntada por "valor" a minúsculas.
                    if ((strcmp(valor, "eul") != 0) && (strcmp(valor, "heun") != 0) && (strcmp(valor, "rk4") != 0) && (strcmp(valor, "abm") != 0) && (strcmp(valor, "imp") != 0))
                    {
                        codErr[nErr] = 4;
                        msjErr[nErr] = "debe indicar un formato de archivo correspondiente a algunas de las aplicaciones para representación gráfica soportadas.";
                        nErr++;
                    }
                    strcpy(metodoNum, valor);
                    break;


                case 'f':

                    if ((strcmp(valor, "Gr") != 0) && (strcmp(valor, "gp") != 0))
                    {
                        codErr[nErr] = 4;
                        msjErr[nErr] = "debe indicar un formato de archivo correspondiente a algunos de los softwares de representación gráfica soportados.";
                        nErr++;
                    }
                    strcpy(formArch, valor);
                    break;


                case 'u':

                    usuario = valor;
                    break;


                case 'c':

                    comentarios = valor;
                    break;

                }
                ident[i - 1] = *argv[i];
            }

            /* Se verifica que no existen argumentos mutuamente excluyentes ("h" y "n"). */
            if ((strchr(ident, 'n') != NULL) && (strchr(ident, 'h') != NULL))
            {
                codErr[nErr] = 5;
                msjErr[nErr] = "los argumentos \"h\" y \"n\" son mutuamente excluyentes. Emplee sólo uno de estos parámetros.";
                nErr++;
            }

            if ((strchr(ident, 'e') != NULL) && (strchr(ident, 'T') != NULL))
            {
                codErr[nErr] = 5;
                msjErr[nErr] = "los argumentos \"e\" y \"T\" son mutuamente excluyentes. Emplee sólo uno de estos parámetros.";
                nErr++;
            }

            /* Verifica que el usuario haya introducido el argumento material. */
            //      Debería crear una función que verifique si el material existe en la biblioteca componentes.txt y retorne un valor NULL en caso de que no exista
            //      o un puntero a la posición dentro del archivo donde se encuentra (dependiendo del tipo de archivo).

            if (strchr(ident, 'm') != NULL)
            {
                /* Obtención de los parámetros del material superconductor. */

                sc = obtenerParametros(material);
                /* Verifica si se está realizando un barrido de temperaturas o sólo calculando los parámetros para una temperatura puntual. */
                if (strchr(ident, 'T') == NULL) tempF = *sc->Tc;
                else tempF = temp;

                /* Validación de la temperatura. */
                if ((temp < 0) || (temp > tempF))
                {
                    codErr[nErr] = 4;
                    msjErr[nErr] = "la temperatura debe ser un número real positivo o cero y no mayor a la temperatura crítica del material: 0 ≤ T ≤ Tc.";
                    nErr++;
                }
            }
            else
            {
                codErr[nErr] = 7;
                msjErr[nErr] = "El argumento material no es opcional. Debe especificar el material a través del identificador \"m\".";
                nErr++;
            }

            /* Impresión de mensajes de error e interrupción del programa si corresponde. */
            for(i = 0; i < nErr; i++) error(codErr[i], msjErr[i]);
            if (nErr > 0)
            {
                printf("\n\nPara mayor información consulte la ayuda del programa, empleando el identificador \"?\".\n\n");
                exit(1);
            }
        }/* Fin de bloque - Interpretación de comandos. */

        /** Nombre del archivo de salida. */
        /* Reservamos espacio para el nombre del archivo si éste no fue definido; el número de caracteres de la extensión es de 4. No
         * creemos necesario usar la función 'sizeof' porque para cualquier arquitectura, el tamaño de un dato tipo 'char' debería ser 1.
         */
        if (strchr(ident, 'a') == NULL)
        {
            nomArch = (char *) malloc(strlen(material) + 6);
            if (strchr(ident, 'T') == NULL) strcat(strcpy(nomArch, material), ".dat");
            else strcat(strcpy(nomArch, material), "_v.dat");
        }



        /* Formato de archivo de salida. */
        if (strcmp(formArch, "gp") == 0)
        {
            nomSoftGraf = "gnuplot";
            charComent = '#';
            igs = popen(GNUPLOT,"w");
        }
        else if (strcmp(formArch, "Gr") == 0)
        {
            nomSoftGraf = "Grace";
            charComent = '#';

        }



        /* Definición del paso y el número de cálculos realizados. */
        if (strchr(ident, 'h') == NULL) h = (tempF - temp) / ni;
        else ni = ceil((tempF - temp) / h);
        /* Se almacenan los valores de h y ni, por si acaso sean utilizados en la solución numérica. */
        auxh = h;
        auxni = ni;

        /* Evita que ocurra un bucle infinito si la estimación se realiza en un intervalo muy cerrado del núcleo del vórtice (h ≈). */
        if (h < 1e-9)
        {
            h = 1;
            ni = 1;
        }



        /* Nombre del método seleccionado por el usuario. */
        if (strcmp(metodoNum, "eul") == 0) nomMetodoNum = "Euler hacia adelante";
        else if (strcmp(metodoNum, "heun") == 0) nomMetodoNum = "Heun";
        else if (strcmp(metodoNum, "rk4") == 0) nomMetodoNum = "Runge-Kutta de cuarto orden (RK4)";
        else if (strcmp(metodoNum, "abm") == 0) nomMetodoNum = "Adams-Bashforth-Moulton de cuarto orden (ABM)";
        else if (strcmp(metodoNum, "imp") == 0) nomMetodoNum = "Método implícito";
        /* Obtención de la fecha y hora del sistema. */
        segundos = time(NULL);
        tiempo = localtime(&segundos);


        /** Generación del archivo de datos. */
        archivo = fopen(nomArch, "wt");
        fprintf(archivo, "%c Archivo de datos generado con GNU SupercondUC v0.001.\n", charComent);
        fprintf(archivo, "%c Formato compatible de archivo: %s.\n", charComent, nomSoftGraf);
        fprintf(archivo, "%c Usuario: %s.\n%c\n", charComent, usuario, charComent);
        fprintf(archivo, "%c Material: %s.\n", charComent, material);
        /* Imprime el resto sobre la información del archivo si se está realizando un barrido de la temperatura */
        if (strchr(ident, 'T') == NULL)
        {
            fprintf(archivo, "%c\n%c Número de puntos calculados: %d.\n", charComent, charComent, ni);
            /* Formato de representación de los datos: temp s alfa beta gamma gammaT mi kappa1 n */
            sprintf(formato, "\n%%.%dg %%d %%.%dg %%.%dg %%.%dg %%.%dg %%.%dg %%.%dg %%d", precision, precision, precision, precision, precision, precision, precision);
        }
        else
        {
            fprintf(archivo, "%c Temperatura: %g K.\n%c\n", charComent, temp, charComent);
            printf("\nHola\n");
            fprintf(archivo, "%c Método numérico: %s.\n", charComent, nomMetodoNum);
            fprintf(archivo, "%c Número de iteraciones previstas: %d.\n", charComent, auxni);
            /* Formato de representación de los datos: r f1 f2 a */
//            sprintf(formato, "\n%%.%dg %%.%dg %%.%dg %%.%dg %%.%dg", precision, precision, precision, precision, precision);
            sprintf(formato, "\n%%.%dg %%.%dg %%.%dg %%.%dg", precision, precision, precision, precision);
        }
        fprintf(archivo, "%c Fecha (día/mes/año): %02i/%02i/%04i.\n", charComent, tiempo->tm_mday, tiempo->tm_mon+1, tiempo->tm_year+1900);
        fprintf(archivo, "%c Hora: %02ih %02im %02is.\n%c\n", charComent, tiempo->tm_hour, tiempo->tm_min, tiempo->tm_sec, charComent);
        fprintf(archivo, "%c Comentarios: %s.\n", charComent, comentarios);

        /**< Obtención de los parámetros del material superconductor. */
        /* Cálculo de los parámetros normalizados que son "independientes" de la temperatura. */
        w = pow(*(sc->Vf + 1) / *sc->Vf, 2);
        k = sqrt(1 - 4 * sc->gammaA * sc->gammaA);
        beta = (1 + k) / (1 - k) * pow(*(sc->Tc + 1) / *(sc->Tc + 2), 2);// beta2/beta1
        mi =  beta * w;//*sc->m / *(sc->m + 1);
        gammaT = 0;// beta12/beta1
        kappa1 = *sc->lambda0 / *sc->xi0;// lambda10/xi10
        n = 1;// Vorticidad
        do
        {
            /* Variables auxiliares para optimizar los cálculos. */
            aux1 = 1 - temp / *(sc->Tc + 1);
            aux2 = 1 - temp / *(sc->Tc + 2);
            k = sqrt(1 - 4 * sc->gammaA * sc->gammaA / (aux1 * aux2));
            aux3 = fabs((-1 + k) * aux1);

            /* Cálculo de los parámetros normalizados que dependen de la temperatura. */
            s = -1;//copysign(1, -aux1);
            alfa = (-1 - k) * aux2 / aux3;
            gamma = 2 * sc->gammaA / aux3;

            temp2 = temp;
            if (escala > 0) temp2 /= *sc->Tc;
            /* Almacenamiento en el archivo de datos. */
            if (strchr(ident, 'T') == NULL)
            {
                if (!(isfinite(alfa)) || !(isfinite(beta)) || !(isfinite(gamma)))
                    fprintf(archivo,"\n%c Error de cálculo: valor indeterminado.", charComent);
                else
                    fprintf(archivo, formato, temp2, s, alfa, beta, gamma, gammaT, mi, kappa1, n);
            }

            temp = temp + h;
        }
        while (temp <= tempF);



        if (strchr(ident, 'T') != NULL)
        {
            h = auxh;
            ni = auxni;
            printf("\nParámetros: %g %d %g %g %g %g %g %g %d\n", temp, s, alfa, beta, gamma, gammaT, mi, kappa1, n);
            /** Condiciones iniciales (en r = 0+). En la región del core. */
            b = (alfa - mi*s + sqrt((alfa - mi*s)*(alfa - mi*s) + 4*mi*gamma*gamma))/(2*gamma);
            gammaT = (mi*b*b - beta)/(mi - b*b);
            printf("\n1.- La inversa de la relación de pendientes b = b1 / b2, es igual a: %.17g.\n", 1/b);
            /* Cálculo de B0. */
            A = 24*b*pow(kappa1,4)/mi*(gammaT*mi + beta - (gammaT + mi)*b*b);
            B = -16*pow(kappa1,3)/mi*((kappa1*kappa1 + 1/4)*b*b*b*b + ((-gammaT*s - alfa)*kappa1*kappa1 + mi*s/4 - alfa/4)*b*b*b - (gamma - 1)/4*(-4*gammaT*kappa1*kappa1 + mi)*b*b + ((alfa*gammaT + beta*s)*kappa1*kappa1 + mi*(mi*s - alfa)/4)*b - gamma/4*(4*beta*kappa1*kappa1 + mi*mi));
            C = 8*pow(kappa1,4)/mi/mi*((kappa1*kappa1 + 1/4)*(mi*s + alfa)*b*b*b*b + (mi*mi/4 - kappa1*kappa1*(gamma*gamma + gammaT)*mi + (-gammaT*gamma*gamma - alfa*alfa)*kappa1*kappa1 - alfa*alfa/4)*b*b*b - (gamma - 1)/4*(-4*gammaT*kappa1*kappa1 + mi)*(mi*s + alfa)*b*b + (mi*mi*mi/4 + ((gammaT*gamma*gamma + beta)*kappa1*kappa1 - alfa*alfa/4)*mi + kappa1*kappa1*(alfa*alfa*gammaT + beta*gamma*gamma))*b - gamma/4*(4*beta*kappa1*kappa1 + mi*mi)*(mi*s + alfa));
            B01 = (-B - sqrt(B*B - 4*A*C))/A/2;
            B02 = (-B + sqrt(B*B - 4*A*C))/A/2;
            if (B01 > 0) B0 = B01;
            else B0 = B02;
            printf("\nB01 = %.17g, B02 = %.17g.\n", B01, B02);
            /* Método 1: considerando que sólo existen b1 y b2 (parámetro de orden cuyo comportamiento es lineal en el núcleo). */
            //d_f1 = kappa1*B0*b*sqrt((1 - b)/((1 + 4*kappa1*kappa1)*b*b*b + (4*kappa1*kappa1*gammaT/mi - 1)*b*b + (mi - 4*kappa1*kappa1*gammaT)*b - mi - 4*kappa1*kappa1*gammaT*beta/mi));//f10/sqrt(2)/xi_ef;
            //d_f2 = d_f1/b;
            /* Método 2: considerando un parámetro de orden cuyo comportamiento es: fi = bi*r + di*r^3.
             * Requiere resolver un sistema no lineal para las incógnitas b1 y b2, empleando métodos numéricos.
             * Se facilita el cálculo, procediendo de manera similar con el cálculo de la relación entre los parámetros de orden; para este caso: b = b1 / b2.
             * Resulta necesario entonces, hallar las raíces de un polinomio de orden cuatro.
             */
            h10 =  kappa1*kappa1*(kappa1*kappa1*(gamma*gamma/mi + 1) + B0*(3*B0 - 2*kappa1*s));
            h11 =  8*kappa1*kappa1 + 2;
            h12 =  2*m - 8*kappa1*kappa1*gammaT;
            h13 =  kappa1*kappa1*kappa1*gamma*(kappa1*(s + alfa/mi) - 2*B0);
            h20 =  kappa1*kappa1/mi*(kappa1*kappa1*(alfa*alfa/mi + gamma*gamma) + B0*(3*mi*B0 - 2*kappa1*alfa));
            h21 =  2 - 8*kappa1*kappa1*gammaT/mi;
            h22 =  8*kappa1*kappa1*beta/mi + 2*mi;
            h23 =  h13/mi;
            //b = raizRealPos(h11*h23, h10*h21 - h11*h20, h12*h23 - h13*h21, h10*h22 - h12*h20, -h13*h22);
            //printf("\n2.- La inversa de la relación de pendientes b = b1 / b2, es igual a: %.17g.\n", 1/b);
            d_f1 = 1e5;
            d_f2 = d_f1/b;
            d_a = 0;
            d_f1a = d_f1;
            d_aa = d_a;
            /* Método 3 */
            //f20/sqrt(2)/xi_ef;
            f1a = 0;
            f2a = 0;
            aa = 1 - 1/n;


            /**< En la región del bulk. */
            etaR = raizRealPos(beta*gamma, -s*beta - alfa*gammaT, gamma*(gammaT - 1), alfa + s, -gamma);
            f1b = sqrt((gamma*etaR - s)/(1 - gammaT*etaR*etaR));
            f2b = etaR*f1b;
            ab = 0;


            /* Definición del intervalo donde se va a realizar el cálculo numérico. */
            menorLambda0 = *sc->lambda0 = f1b / d_f1; // Realmente así no se calcula lambda sino xi.
            *(sc->lambda0 + 1) = f2b / d_f2;
            printf("\nb1 = %.17g, b2 = %.17g.\n", d_f1, d_f2);
            printf("\nf10 = %.17g, f20 = %.17g.\n", f1b, f2b);
            printf("\nxi1 = %.17g, xi2 = %.17g.\n", f1b / d_f1, f2b / d_f2);
            printf("temp = %lg, s = %d, alfa = %lg, beta = %lg, gamma = %lg, gammaT = %lg, mi = %lg, kappa1 = %lg, n = %d", temp, s, alfa, beta, gamma, gammaT, mi, kappa1, n);

            /* Se determina la menor de las longitudes de coherencias */
            for(i = 1; i < nBandas; i++) if (*(sc->xi0 + i) < menorLambda0) menorLambda0 = *(sc->lambda0 + i);
            printf("\nMenorLambda = %.17g.\n", menorLambda0);
            rf = factor * menorLambda0;
            /* Valor final de la variable independiente. */
            r = rf/ni;//menorXi0/ni;


            /* Definición del paso y el número de cálculos realizados. */
            if (strchr(ident, 'h') == NULL) h = (rf - r) / ni;
            else ni = ceil(fabs(rf - r) / h);

            printf("\nr = %g, rf = %g\n", r, rf);





            /** AQUÍ ESTÁ EL PROBLEMA. SI h ES MUY PEQUEÑO, LA CAPACIDAD COMPUTACIONAL NO PERMITE VARIAR LOS PASOS.
                PODRÍA RESOLVERSE ESTE PROBLEMA, UTILIZANDO UNA VARIABLE x QUE ESCALE r Y AJUSTAR LAS ECUACIONES
                PARA QUE LOS CÁLCULOS SEAN CORRECTOS. */
            /* Evita que ocurra un bucle infinito si la estimación se realiza en un intervalo muy cerrado del núcleo del vórtice (h ≈). */
//            if (h < 1e-9)
//            {
//                h = 1;
//                ni = 1;
//            }





            /** Determinación de las longitudes (coherencia y profundidad de penetración) efectivas. */
            lambda_ef = 1/(f10*sqrt(1 + m*etaR*etaR));
            xi_ef = *sc->xi0; // Falta por definir.


            /** Parámetros constantes del sistema de EDO a resolver (optimización de cálculos). */
            c2 = n*n;
            c3 = kappa1*kappa1;
            c1 = c3*gammaT;
            c4 = c3*gamma;
            aux = c3/mi;
            c5 = aux*gammaT;
            c6 = aux*alfa;
            c7 = aux*beta;
            c8 = aux*gamma;
            c9 = c3*s;


            /** Error relativo máximo. */
            prctjErr = .5; // 5% de error relativo.


            f1 = f1a;
            f2 = f2a;
            a = aa;
            B = 1; // Ficticio. Tengo que ver cómo lo calculo o lo estimo...

            /* Primera línea del archivo de datos. */
            //fprintf(archivo, formato, r, f1, f2, a, B);
            fprintf(archivo, formato, r, f1, f2, a);


            switch(metodoNum[0])
            {
            case 'e':

                /** Método de Euler hacia adelante. */

                for(i=1; i<ni; i++)
                {
                    /* Minimizando los cálculos repetitivos para optimizar el tiempo de ejecución. */
                    aux = 1/r;
                    f1_2 = f1*f1;
                    f2_2 = f2*f2;
                    a_1r = (a - 1)*aux;
                    na1r_2 = c2*a_1r*a1_r;

                    d2_f1 = -aux*d_f1 + na1r_2*f1 + c3*(s*f1 + f1*f1_2 - gamma*f2 - gammaT*f1*f2_2);
                    d2_f2 = -aux*d_f2 + na1r_2*f2 + c3/mi*(alfa*f2 + beta*f2*f2_2 - gamma*f1 - gammaT*f1_2*f2);
                    d2_a = aux*d_a + (a - 1)*(f1_2 + mi*f2_2);

                    //printf("\nf1 = %g, f2 = %g, a = %g, d2_f1 = %g, d2_f2 = %g, d2_a = %g.", f1, f2, a, d2_f1, d2_f2, d2_a);

//                d2_f1 = -aux*d_f1 + na1r_2 + f1*(c9 + c3*f1_2 - c1*f2_2) - c4*f2;
//                d2_f2 = -aux*d_f2 + na1r_2 + f2*(c6 + c7*f2_2 - c5*f1_2) - c8*f1;
//                d2_a = aux*(d_a + a_1r/aux*(f1_2 + mi*f2_2));

                    f1 += h*d_f1;
                    f2 += h*d_f2;
                    a += h*d_a;

                    if (!(isfinite(f1)) || !(isfinite(f2)) || !(isfinite(a)))
                    {
                        fprintf(archivo,"\n\n%c Ha ocurrido una excepción y los cálculos se detuvieron en este punto.", charComent);
                        fprintf(archivo,"\n%c Consulte las líneas anteriores para analizar el comportamiento de los datos.", charComent);
                        fprintf(archivo,"\n%c Iteraciones completadas: %d.", charComent, i);
                        break;
                    }

                    r += h;

                    fprintf(archivo, formato, r, f1, f2, a);

                    d_f1 += h*d2_f1;
                    d_f2 += h*d2_f2;
                    d_a += h*d2_a;
                }
                break;

            case 'h':

                /** Método de Heun. */

                for(i=1; i<ni; i++)
                {
                    /* Minimizando los cálculos repetitivos para optimizar el tiempo de ejecución. */
                    aux = 1/r;
                    f1_2 = f1*f1;
                    f2_2 = f2*f2;
                    a_1r = (a - 1)*aux;
                    na1r_2 = c2*a_1r*a1_r;

                    d2_f1 = -aux*d_f1 + na1r_2 + f1*(c9 + c3*f1_2 - c1*f2_2) - c4*f2;
                    d2_f2 = -aux*d_f2 + na1r_2 + f2*(c6 + c7*f2_2 - c5*f1_2) - c8*f1;
                    d2_a = aux*d_a + (a - 1)*(f1_2 + mi*f2_2);

                    for(j=1; j < 3; j++)
                    {
                        if (j == 1)
                        {
                            f1n = f1;
                            d_f1n = d_f1;
                            f1 += h*d_f1;
                            d_f1 += h*d2_f1;

                            f2n = f2;
                            d_f2n = d_f2;
                            f2 += h*d_f2;
                            d_f2 += h*d2_f2;

                            an = a;
                            d_an = d_a;
                            a += h*d_a;
                            d_a += h*d2_a;

                            r += h;
                        }
                    }

                    f1 = f1n + h*(d_f1n + d_f1)/2;
                    f2 = f2n + h*(d_f2n + d_f2)/2;
                    a = an + h*(d_an + d_a)/2;

                    if (!(isfinite(f1)) || !(isfinite(f2)) || !(isfinite(a)))
                    {
                        fprintf(archivo,"\n\n%c Ha ocurrido una excepción y los cálculos se detuvieron en este punto.", charComent);
                        fprintf(archivo,"\n%c Consulte las líneas anteriores para analizar el comportamiento de los datos.", charComent);
                        fprintf(archivo,"\n%c Iteraciones completadas: %d.", charComent, i);
                        break;
                    }

                    fprintf(archivo, formato, r, f1, f2, a);
                }
                break;

            case 'r':
                /** Método de Runge-Kutta de cuarto orden. */

                /* Definiendo los parámetros de Runge-Kutta de cuarto orden. */
                w1 = 1;
                w2 = 2;
                w3 = 2;
                w4 = 1;
                w = 1/(w1 + w2 + w3 + w4);

                a1 = .5;
                a2 = .5;
                a3 = 1;

                b1 = .5;
                b2 = 0;
                b3 = .5;
                b4 = 0;
                b5 = 0;
                b6 = 1;

                for (i = 1; i < ni; i++)
                {
                    f1n = f1;
                    d_f1n = d_f1;
                    f2n = f2;
                    d_f2n = d_f2;
                    an = a;
                    d_an = d_a;
                    rn = r;

                    for(j = 1; j < 5; j++)
                    {
                        /* Minimizando los cálculos repetitivos para optimizar el tiempo de ejecución. */
                        aux = 1/r;
                        f1_2 = f1*f1;
                        f2_2 = f2*f2;
                        a_1r = (a - 1)*aux;
                        na1r_2 = c2*a_1r*a1_r;

                        d2_f1 = -aux*d_f1 + na1r_2*f1 + c3*(s*f1 + f1*f1_2 - gamma*f2 - gammaT*f1*f2_2);
                        d2_f2 = -aux*d_f2 + na1r_2*f2 + c3/mi*(alfa*f2 + beta*f2*f2_2 - gamma*f1 - gammaT*f1_2*f2);
                        d2_a = aux*d_a + (a - 1)*(f1_2 + mi*f2_2);

                        switch(j)
                        {
                        case 1:

                            k1[0] = d_f1n*h;
                            k1[1] = d2_f1*h;
                            k1[2] = d_f2n*h;
                            k1[3] = d2_f2*h;
                            k1[4] = d_an*h;
                            k1[5] = d2_a*h;

                            f1 = f1n + b1*k1[0];
                            d_f1 = d_f1n + b1*k1[1];
                            f2 = f2n + b1*k1[2];
                            d_f2 = d_f2n + b1*k1[3];
                            a = an + b1*k1[4];
                            d_a = d_an + b1*k1[5];

                            r = rn + a1*h;

                            break;

                        case 2:

                            k2[0] = (d_f1n + b1*k1[1])*h;
                            k2[1] = d2_f1*h;
                            k2[2] = (d_f2n + b1*k1[3])*h;
                            k2[3] = d2_f2*h;
                            k2[4] = (d_an + b1*k1[5])*h;
                            k2[5] = d2_a*h;

                            f1 = f1n + b2*k1[0] + b3*k2[0];
                            d_f1 = d_f1n + b2*k1[1] + b3*k2[1];
                            f2 = f2n + b2*k1[2] + b3*k2[2];
                            d_f2 = d_f2n + b2*k1[3] + b3*k2[3];
                            a = an + b2*k1[4] + b3*k2[4];
                            d_a = d_an + b2*k1[5] + b3*k2[5];

                            r = rn + a2*h;

                            break;

                        case 3:

                            k3[0]=(d_f1n + b2*k1[1]+b3*k2[1])*h;
                            k3[1]=d2_f1*h;
                            k3[2]=(d_f2n+b2*k1[3]+b3*k2[3])*h;
                            k3[3]=d2_f2*h;
                            k3[4]=(d_an+b2*k1[5]+b3*k2[5])*h;
                            k3[5]=d2_a*h;

                            f1 = f1n + b4*k1[0] + b5*k2[0] + b6*k3[0];
                            d_f1 = d_f1n + b4*k1[1]+b5*k2[1]+b6*k3[1];
                            f2=f2n+b4*k1[2]+b5*k2[2]+b6*k3[2];
                            d_f2=d_f2n+b4*k1[3]+b5*k2[3]+b6*k3[3];
                            a=an+b4*k1[4]+b5*k2[4]+b6*k3[4];
                            d_a=d_an+b4*k1[5]+b5*k2[5]+b6*k3[5];

                            r=rn+a3*h;

                            break;

                        case 4:

                            k4[0]=(d_f1n+b4*k1[1]+b5*k2[1]+b6*k3[1])*h;
                            k4[1]=d2_f1*h;
                            k4[2]=(d_f2n+b4*k1[3]+b5*k2[3]+b6*k3[3])*h;
                            k4[3]=d2_f2*h;
                            k4[4]=(d_an+b4*k1[5]+b5*k2[5]+b6*k3[5])*h;
                            k4[5]=d2_a*h;

                            r=rn+h;
                            break;
                        }
                    }

                    f1 = f1n + w*(w1*k1[0] + w2*k2[0] + w3*k3[0] + w4*k4[0]);
                    f2 = f2n + w*(w1*k1[2] + w2*k2[2] + w3*k3[2] + w4*k4[2]);
                    a = an + w*(w1*k1[4] + w2*k2[4] + w3*k3[4] + w4*k4[4]);

//                    printf("%g %g %g %g", r, f1, f2, d2_f1);

                    if (!(isfinite(f1)) || !(isfinite(f2)) || !(isfinite(a)))
                    {
                        fprintf(archivo,"\n\n%c Ha ocurrido una excepción y los cálculos se detuvieron en este punto.", charComent);
                        fprintf(archivo,"\n%c Consulte las líneas anteriores para analizar el comportamiento de los datos.", charComent);
                        fprintf(archivo,"\n%c Iteraciones completadas: %d.", charComent, i);
                        break;
                    }
                    B0=1;
                    C = 1 - 1/n*cosh(.5*d_f1a/b*sqrt(mi + b*b)*r*r) + kappa1*B0*b/n/d_f1a/sqrt(mi + b*b)*sinh(.5*d_f1a/b*sqrt(mi + b*b)*r*r);
                    //fprintf(archivo, formato, r, f1, f2, a, B);
                    d_f1 = d_f1n + w*(w1*k1[1] + w2*k2[1] + w3*k3[1] + w4*k4[1]);
                    d_f2 = d_f2n + w*(w1*k1[3] + w2*k2[3] + w3*k3[3] + w4*k4[3]);
                    d_a = d_an + w*(w1*k1[5] + w2*k2[5] + w3*k3[5] + w4*k4[5]);
                    B = n/kappa1*d_a;

                    fprintf(archivo, formato, r, f1, f2, B);


                }
                break;

            case 'i':

                /** Método implícito. */


                break;

            case 'a':

                /** Método Adams-Bashforth-Moulton. */

                while (fabs(1 - f1 / f1b)*100 > prctjErr)
                {

                    f1 = f1a;
                    f2 = f2a;
                    a = aa;

                    for(i = 1; i < ni; i++)
                    {
                        /* Minimizando los cálculos repetitivos para optimizar el tiempo de ejecución. */
                        aux = 1/r;
                        f1_2 = f1*f1;
                        f2_2 = f2*f2;
                        a_1r = --a*aux;
                        na1r_2 = c2*a_1r*a1_r;

                        d2_f1 = -aux*d_f1 + na1r_2 + f1*(c9 + c3*f1_2 - c1*f2_2) - c4*f2;
                        d2_f2 = -aux*d_f2 + na1r_2 + f2*(c6 + c7*f2_2 - c5*f1_2) - c8*f1;
                        d2_a = aux*(d_a + a_1r*(f1_2 + mi*f2_2));

                        f1 += h*d_f1;
                        f2 += h*d_f2;
                        a += h*d_a;

                        if (!(isfinite(f1)) || !(isfinite(f2)) || !(isfinite(a)))
                        {
                            fprintf(archivo,"\n\n%c Ha ocurrido una excepción y los cálculos se detuvieron en este punto.", charComent);
                            fprintf(archivo,"\n%c Consulte las líneas anteriores para analizar el comportamiento de los datos.", charComent);
                            fprintf(archivo,"\n%c Iteraciones completadas: %d.", charComent, i);
                            break;
                        }

                        r += h;

                        fprintf(archivo, formato, r, f1, f2, a);

                        d_f1 += h*d2_f1;
                        d_f2 += h*d2_f2;
                        d_a += h*d2_a;
                    }
                }
                break;
            }
        }



        /** Gráfico de los resultados. */

        /* Selección del software de gráficos. */
        // gnuplot
        if (strcmp(formArch, "gp") == 0)
        {

            if ((igs = popen(GNUPLOT,"w")) == NULL)
            {
                printf("El software gnuplot parece no haber sido instalado. Si ya está instalado, pruebe removiendo y haciendo una instalación nueva del paquete.\n");
                exit(0);
            }
            else
            {
                fprintf(igs, "set termopt enhanced\n");
                fprintf(igs, "set grid\n");
                if (strchr(ident, 'T') == NULL)
                {
                    /* Impresión de textos */
                    fprintf(igs, "set title 'Parámetros de las ecuaciones de GL\n");
                    if (escala < 1) fprintf(igs, "set xlabel 'T'\n");
                    else fprintf(igs, "set xlabel 'T / T_c'\n");
                    fprintf(igs, "set ylabel 'Superconductor: %s'\n", tipografia(material, formArch));

                    fprintf(igs, "plot '%s' using 1:2 with lines title 's', '' using 1:3 with lines title '{/Symbol=a}', '' using 1:4 with lines title '{/Symbol=b}', '' using 1:5 with lines title '{/Symbol=g}', '' using 1:6 with lines title '~{/Symbol=g}{.5-}', '' using 1:7 with lines title '{/Symbol=m}', '' using 1:8 with lines title '{/Symbol=k}_1', '' using 1:9 with lines title 'n'\n", nomArch);
                    fclose(igs);
                }
                else
                {
                    /* Impresión de textos */
                    fprintf(igs, "set title 'Solución a las ecuaciones de Ginzburg-Landau'\n");
                    fprintf(igs, "set xlabel 'r'\n");
                    fprintf(igs, "set ylabel 'Superconductor: %s \\@ T = %g K'\n", tipografia(material, formArch), temp2);

//                    fprintf(igs, "plot '%s' using 1:4 with lines title 'a(r)'\n", nomArch);
                    fprintf(igs, "plot '%s' using 1:2 with lines title '|{/Symbol=y}_1(r)|', '' using 1:3 with lines title '|{/Symbol=y}_2(r)|', '' using 1:4 with lines title 'a(r)'\n", nomArch);
//                    fprintf(igs, "plot '%s' using 1:2 with lines title '|{/Symbol=y}_1(r)|', '' using 1:3 with lines title '|{/Symbol=y}_2(r)|', '' using 1:4 with lines title 'a(r)', '' using 1:5 with lines title 'B(r)'\n", nomArch);
                    fclose(igs);
                }
            }
        }



        //Grace
        else if (strcmp(formArch, "Gr") == 0)
        {
            fprintf(archivo, "\n&"); // El símbolo & marca el final de un archivo con formato para Grace.
            GraceRegisterErrorFunction(my_error_function);

            /* Start Grace with a buffer size of 2048 and open the pipe */
            if (GraceOpen(2048) == -1)
            {
                fprintf(stderr, "No se puede ejecutar Grace. \n");
                exit(EXIT_FAILURE);
            }
            if (strchr(ident, 'T') == NULL)
            {
                /* Impresión de textos */
                GracePrintf("title \"Par\\ca\\Cmetros de las ecuaciones de GL\"");
                if (escala < 1) GracePrintf("xaxis label \"T\"");
                else GracePrintf("xaxis label \"T / T\\sc\"");
                GracePrintf("subtitle \"Superconductor: %s\"", tipografia(material, formArch));

                GracePrintf("s0 legend \"\\0s\"");
                GracePrintf("s1 legend \"\\xa\"");
                GracePrintf("s2 legend \"\\xb\"");
                GracePrintf("s3 legend \"Acoplamiento (\\xg)\"");
                GracePrintf("s4 legend \"\\xg\\h{-0.45}\\v{0.45}~\"");
                GracePrintf("s5 legend \"Masa efectiva (\\xm)\"");
                GracePrintf("s6 legend \"Par\\ca\\Cmetro de GL (\\xk\\s1\\N)\"");
                GracePrintf("s7 legend \"\Vorticidad (n)\"");
            }
            else
            {
                /* Impresión de textos */
                GracePrintf("title \"Soluci\\cs\\Cn a las ecuaciones de Ginzburg-Landau\"");
                GracePrintf("xaxis label \"r\"");
                GracePrintf("subtitle \"Superconductor: %s @ T = %g K\"", tipografia(material, formArch), temp2);

                GracePrintf("s0 legend \"\\x|y\\s1\\N\\f{}(r)\\x|\"");
                GracePrintf("s1 legend \"\\x|y\\s2\\N\\f{}(r)\\x|\"");
                GracePrintf("s2 legend \"a(r)\"");
                //GracePrintf("s3 legend \"B(r)\"");
            }

            GracePrintf("read nxy \"%s\"", nomArch);
            GracePrintf("legend on");
            GracePrintf("autoscale");
        }



    }
    fclose(archivo);
    return 0;
}

