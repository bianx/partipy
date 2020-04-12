#include <stdio.h>
#include <stdlib.h>
#include <barnes-hut.h>

static const char *me = "barnes_hut/example/print";
enum { SIZE = 999 };

static void
usg(void)
{
    fprintf(stderr, "%s points query\n", me);
    exit(1);
}

int
main(int argc, char **argv)
{
    char line[SIZE];
    double xc;
    double yc;
    double w;
    double *x;
    double *y;
    double mass;
    long cap;
    long n;
    struct BarnesHut *barnes_hut;

    (void) argc;

    xc = 0;
    yc = 0;
    w = 2;
    mass = 1;

    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }

    cap = 1;
    if ((x = malloc(cap * sizeof(*x))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((y = malloc(cap * sizeof(*y))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((barnes_hut = barnes_hut_ini(xc, yc, w)) == NULL) {
        fprintf(stderr, "%s:%d: barnes_hut_ini failed\n", __FILE__,
                __LINE__);
        exit(1);
    }

    n = 0;
    while (fgets(line, SIZE, stdin) != NULL) {
        if (n == cap) {
            cap *= 2;
            x = realloc(x, cap * sizeof(*x));
            if (x == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(1);
            }
            y = realloc(y, cap * sizeof(*y));
            if (y == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(1);
            }
        }
        if (sscanf(line, "%lf %lf\n", &x[n], &y[n]) != 2) {
            fprintf(stderr, "%s: fail to parse '%s'\n", me, line);
            exit(1);
        }
        n++;
    }

    long i;

    for (i = 0; i < n; i++)
        barnes_hut_insert(barnes_hut, x[i], y[i], mass, i);
    barnes_hut_print(barnes_hut, stdout);

    free(x);
    free(y);
    barnes_hut_fin(barnes_hut);
}
