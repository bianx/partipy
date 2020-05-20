#include <stdio.h>
#include <stdlib.h>
#include <tree.h>

static const char *me = "tree/example/print";
enum { SIZE = 999 };

static void
usg(void)
{
    fprintf(stderr, "%s < points\n", me);
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
    long i;
    struct Tree *tree;
    struct TreeParam param;

    (void) argc;

    xc = 0;
    yc = 0;
    w = 2;
    mass = 1;
    param.cap = 0;

    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        case 'c':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -c\n", me);
                exit(2);
            }
            param.cap = atoi(argv[0]);
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }

    if (param.cap == 0) {
        fprintf(stderr, "%s: -c is not set\n", me);
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
    if ((tree = tree_ini(param, xc, yc, w)) == NULL) {
        fprintf(stderr, "%s:%d: tree_ini failed\n", __FILE__, __LINE__);
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

    for (i = 0; i < n; i++)
        tree_insert(tree, x[i], y[i], mass, i);
    tree_print(tree, stdout);

    free(x);
    free(y);
    tree_fin(tree);
}
