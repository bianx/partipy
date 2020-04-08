#include <stdio.h>
#include <stdlib.h>
#include <quadtree.h>

static const char *me = "quadtree/example/main";
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
    const long *list;
    double B;
    double L;
    double R;
    double T;
    double *x;
    double *y;
    FILE *f;
    int *color;
    long cap;
    long i;
    long m;
    long n;
    struct Quadtree *quadtree;

    (void) argc;

    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }
    if (argv[0] == NULL || argv[1] == NULL) {
        fprintf(stderr, "%s: need two arguments\n", me);
        exit(1);
    }

    if ((f = fopen(argv[0], "r")) == NULL) {
        fprintf(stderr, "%s: fail to open '%s'\n", me, argv[0]);
        exit(1);
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
    if ((quadtree = quadtree_ini()) == NULL) {
        fprintf(stderr, "%s:%d: quadtree_ini failed\n", __FILE__,
                __LINE__);
        exit(1);
    }

    n = 0;
    while (fgets(line, SIZE, f) != NULL) {
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
    if ((color = malloc(n * sizeof(*color))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    for (i = 0; i < n; i++)
        color[i] = 0;
    for (i = 0; i < n; i++)
        if ((quadtree_insert(quadtree, x[i], y[i], i)) != 0) {
            fprintf(stderr, "%s: quadtree_insert failed\n", me);
            exit(1);
        }

    if ((f = fopen(argv[1], "r")) == NULL) {
        fprintf(stderr, "%s: fail to open '%s' in '%s'\n", me, argv[1],
                argv[0]);
        exit(1);
    }
    while (fgets(line, SIZE, f) != NULL) {
        if (sscanf(line, "%lf %lf %lf %lf\n", &L, &R, &B, &T) != 4) {
            fprintf(stderr, "%s: fail to parse '%s' in '%s'\n", me, line,
                    argv[1]);
            exit(1);
        }
        if (quadtree_rectangle(quadtree, L, R, B, T, &m, &list) != 0) {
            fprintf(stderr, "%s: quadtree_rectangle failed\n", me);
            exit(1);
        }
        fprintf(stderr, "%s: %ld\n", me, m);
        for (i = 0; i < m; i++)
            color[list[i]] = 1;
    }
    for (i = 0; i < n; i++)
        printf("%g %g %d\n", x[i], y[i], color[i]);

    fclose(f);
    free(x);
    free(y);
    free(color);
    quadtree_fin(quadtree);
}
