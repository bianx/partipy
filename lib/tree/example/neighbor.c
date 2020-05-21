#include <stdio.h>
#include <stdlib.h>
#include <tree.h>

static const char *me = "tree/example/neighbor";
enum { SIZE = 999 };

static void
usg(void)
{
    fprintf(stderr, "%s -c int -p float float < points\n", me);
    exit(1);
}

int
main(int argc, char **argv)
{
    char line[SIZE];
    double mass;
    double w;
    double wc;
    double *x;
    double xc;
    double xh;
    double xl;
    double xp;
    double *y;
    double yc;
    double yh;
    double yl;
    double yp;
    int Pflag;
    long cap;
    long i;
    long n;
    struct TreeParam param;
    struct Tree *tree;
    struct Node *node;
    struct Node **neighbor;

    (void) argc;

    xc = 0;
    yc = 0;
    w = 2;
    mass = 1;
    param.cap = 0;
    Pflag = 0;

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
        case 'p':
            argv++;
            if (argv[0] == NULL || argv[1] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -p\n", me);
                exit(2);
            }
            xp = atof(argv[0]);
            argv++;
            yp = atof(argv[0]);
            Pflag = 1;
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }

    if (param.cap == 0) {
        fprintf(stderr, "%s: -c is not set\n", me);
        exit(2);
    }
    if (Pflag == 0) {
        fprintf(stderr, "%s: -p is not set\n", me);
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

    if ((tree = tree_ini(param, xc, yc, w)) == NULL) {
        fprintf(stderr, "%s:%d: tree_ini failed\n", __FILE__, __LINE__);
        exit(1);
    }

    for (i = 0; i < n; i++)
        tree_insert(tree, x[i], y[i], mass, i);

    if ((node = tree_node(tree, xp, yp)) == NULL) {
        fprintf(stderr, "%s:%d: tree_box failed\n", __FILE__, __LINE__);
        exit(1);
    }

    if (tree_neighbor(tree, node, &neighbor) != 0) {
        fprintf(stderr, "%s:%d: tree_neighbor failed\n", __FILE__, __LINE__);
        exit(1);
    }
    for (; *neighbor != NULL; neighbor++) {
        if (node_print(*neighbor, stdout) != 0) {
            fprintf(stderr, "%s:%d: node_print failed\n", __FILE__, __LINE__);
            exit(1);
        }
    }
    free(x);
    free(y);
    tree_fin(tree);
}
