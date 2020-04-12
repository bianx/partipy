#include <stdio.h>
#include <stdlib.h>
#include <barnes-hut.h>

static const char *me = "barnes_hut/example/force";
enum { SIZE = 999 };

static void
usg(void)
{
    fprintf(stderr, "%s -p x y -t theta < points\n", me);
    exit(1);
}

static int
function(double x, double y, double *fx, double *fy, void *data)
{
    (void) x;
    (void) y;
    (void) data;

    *fx = 1;
    *fy = 0;
    return 0;
}

int
main(int argc, char **argv)
{
    char line[SIZE];
    double fx;
    double fx0;
    double fy;
    double fy0;
    double *m;
    double mass;
    double theta;
    double *u;
    double *v;
    double w;
    double *x;
    double xc;
    double xp;
    double *y;
    double yc;
    double yp;
    int Pflag;
    int Tflag;
    long cap;
    long cnt;
    long i;
    long j;
    long n;
    struct BarnesHut *barnes_hut;

    (void) argc;

    xc = 0;
    yc = 0;
    w = 2;
    Pflag = Tflag = 0;
    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        case 't':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -t\n", me);
                exit(2);
            }
            theta = atof(argv[0]);
            Tflag = 1;
            break;
        case 'p':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -p\n", me);
                exit(2);
            }
            xp = atof(argv[0]);
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -p\n", me);
                exit(2);
            }
            yp = atof(argv[0]);
            Pflag = 1;
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }
    if (Pflag == 0) {
        fprintf(stderr, "%s: -p is not set\n", me);
        exit(2);
    }
    if (Tflag == 0) {
        fprintf(stderr, "%s: -t is not set\n", me);
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

    if ((u = malloc(n * sizeof(*u))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((v = malloc(n * sizeof(*v))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((m = malloc(n * sizeof(*m))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }

    mass = 1;
    for (i = 0; i < n; i++)
        barnes_hut_insert(barnes_hut, x[i], y[i], mass, i);
    barnes_hut_interaction(barnes_hut, theta, -1, xp, yp, &cnt, u, v, m);

    fx = 0;
    fy = 0;
    for (j = 0; j < cnt; j++) {
        function(xp - u[j], yp - u[j], &fx0, &fy0, NULL);
        fx += m[j] * fx0;
        fy += m[j] * fy0;
    }

    printf("%.16e %.16e\n", fx, fy);

    free(x);
    free(y);
    free(u);
    free(v);
    free(m);
    barnes_hut_fin(barnes_hut);
}
