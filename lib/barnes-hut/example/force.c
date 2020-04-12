#include <math.h>
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
    double f;
    double r;
    double r2;

    (void) data;

    r2 = x * x + y * y;
    if (r2 == 0) {
        fprintf(stderr, "%s: r2 == 0\n", me);
        return 1;
    }

    f = 1 / r2;
    r = sqrt(r2);

    *fx = f * x / r;
    *fy = f * y / r;
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
    double *mass;
    double theta;
    double *u;
    double *v;
    double *x;
    double xp;
    double *y;
    double yp;
    int Pflag;
    int Tflag;
    long cap;
    long cnt;
    long i;
    long n;
    struct BarnesHut *barnes_hut;
    struct BarnesHutInfo *info;
    struct BarnesHutInfo *in;

    (void) argc;

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
    if ((mass = malloc(cap * sizeof(*mass))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    n = 0;
    while (fgets(line, SIZE, stdin) != NULL) {
        if (n == cap) {
            cap *= 2;
            x = realloc(x, cap * sizeof *x);
            if (x == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(1);
            }
            y = realloc(y, cap * sizeof *y);
            if (y == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(1);
            }
            mass = realloc(mass, cap * sizeof *mass);
            if (mass == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(1);
            }
        }
        if (sscanf(line, "%lf %lf %lf\n", &x[n], &y[n], &mass[n]) != 3) {
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
    if ((info = malloc(n * sizeof(*info))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((barnes_hut = barnes_hut_build(n, x, y, mass)) == NULL) {
        fprintf(stderr, "%s:%d: barnes_hut_ini failed\n", __FILE__,
                __LINE__);
        exit(1);
    }
    barnes_hut_interaction(barnes_hut, theta, -1, xp, yp, &cnt, u, v, m);
    barnes_hut_info(barnes_hut, theta, -1, xp, yp, &cnt, info);

    fx = 0;
    fy = 0;
    fprintf(stderr, "%s: cnt = %ld / %ld\n", me, cnt, n);
    for (i = 0; i < cnt; i++) {
        in = &info[i];
        function(xp - in->mx / in->m, yp - in->my / in->m, &fx0, &fy0,
                 NULL);
        fx += in->m * fx0;
        fy += in->m * fy0;
    }
    printf("%.16e %.16e\n", fx, fy);

    fx = 0;
    fy = 0;
    for (i = 0; i < n; i++) {
        function(xp - x[i], yp - y[i], &fx0, &fy0, NULL);
        fx += mass[i] * fx0;
        fy += mass[i] * fy0;
    }
    printf("%.16e %.16e\n", fx, fy);

    free(x);
    free(y);
    free(mass);
    free(info);
    free(u);
    free(v);
    free(m);
    barnes_hut_fin(barnes_hut);
}
