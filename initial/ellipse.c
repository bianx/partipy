#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static char me[] = "initial/ellipse";
static const double pi = 3.141592653589793;
enum { SIZE = 999 };
enum { N_DBL = 100, M_DBL = 100 };

static void
usg(void)
{
    fprintf(stderr, "%s < points > initial\n", me);
    exit(1);
}

static double dblint(double (*f)(double, double, void *), void *, double r,
                     double s, double a, double b);
static double cross(double, double, void *);
static double gauss(double, void *);

struct Core {
  double (*fun)(double, void *);
  void *param;
};
struct CoreParam {
    double delta;
};
struct CrossParam {
  double xi;
  double yi;
  double xj;
  double yj;
  struct Core *core;
};

static struct Core Core[] = {
    { gauss, NULL },
};

static const char *CoreName[] = {
    "gauss",
};

int
main(int argc, char **argv)
{
    (void) argc;
    int i;
    int n;
    double a;
    double b;
    double *x;
    double *y;
    double *A;
    double *B;
    double delta;
    int ncap;
    char line[SIZE];
    int Aflag;
    int Bflag;
    int Dflag;
    struct Core *core;
    struct CoreParam core_param;
    struct CrossParam cross_param;

    core = NULL;
    Aflag = Bflag = Dflag = 0;
    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        case 'd':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -d needs an argument\n", me);
                exit(2);
            }
            delta = atof(argv[0]);
            Dflag = 1;
            break;
        case 'a':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -a needs an argument\n", me);
                exit(2);
            }
            a = atof(argv[0]);
            Aflag = 1;
            break;
        case 'b':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -b needs an argument\n", me);
                exit(2);
            }
            b = atof(argv[0]);
            Bflag = 1;
            break;
        case 'c':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -c needs an argument\n", me);
                exit(2);
            }
            for (i = 0;; i++) {
                if (i == sizeof(Core) / sizeof(*Core)) {
                    fprintf(stderr, "%s: unknown core '%s'\n", me,
                            argv[0]);
                    exit(2);
                }
                if (strncmp(argv[0], CoreName[i], SIZE) == 0) {
                    core = &Core[i];
                    break;
                }
            }
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }
    if (Dflag == 0) {
        fprintf(stderr, "%s: -d is not given\n", me);
        exit(2);
    }
    if (Aflag == 0) {
        fprintf(stderr, "%s: -a is not given\n", me);
        exit(2);
    }
    if (Bflag == 0) {
        fprintf(stderr, "%s: -b is not given\n", me);
        exit(2);
    }
    if (core == NULL) {
        fprintf(stderr, "%s: -c is not given\n", me);
        exit(2);
    }

    ncap = 1;
    if ((x = malloc(ncap * sizeof(*x))) == NULL) {
        fprintf(stderr, "%s: malloc failed\n", me);
        exit(2);
    }
    if ((y = malloc(ncap * sizeof(*y))) == NULL) {
        fprintf(stderr, "%s: malloc failed\n", me);
        exit(2);
    }

    i = 0;
    while (fgets(line, SIZE, stdin) != NULL) {
        if (i == ncap) {
            ncap *= 2;
            if ((x = realloc(x, ncap * sizeof(*x))) == NULL) {
                fprintf(stderr, "%s: realloc failed\n", me);
                exit(2);
            }
            if ((y = realloc(y, ncap * sizeof(*y))) == NULL) {
                fprintf(stderr, "%s: realloc failed\n", me);
                exit(2);
            }
        }
        if (sscanf(line, "%lf %lf", &x[i], &y[i]) != 2) {
            fprintf(stderr, "%s: fail to parse '%s'\n", line, me);
            exit(2);
        }
        i++;
    }
    n = i;
    if ((A = malloc(n * n * sizeof(*A))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }
    if ((B = malloc(n * sizeof(*B))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }

    core_param.delta = delta;
    core->param = &core_param;
    cross_param.core = core;
    cross_param.xi = 0;
    cross_param.yi = 0;
    cross_param.xj = 0;
    cross_param.yj = 0;    
    
    fprintf(stderr, "%g\n", dblint(cross, &cross_param, 0, 2, 0, 2 * pi));

    free(x);
    free(y);
    free(A);
    free(B);
}

static double
cross(double r, double phi, void *p0)
{
  double x;
  double y;
  double xi;
  double yi;
  double xj;
  double yj;
  double ri;
  double rj;
  struct CrossParam *p;
  struct Core *core;

  p = p0;
  core = p->core;
  xi = p->xi;
  yi = p->yi;
  xj = p->xj;
  yj = p->yj;  

  x = r * cos(phi);
  y = r * sin(phi);
  xi = x - p->xi;
  yi = y - p->yi;
  xj = x - p->xj;
  yj = y - p->yj;

  ri = sqrt(xi*xi + yi*yi);
  rj = sqrt(xj*xj + yj*yj);

  return core->fun(ri, core->param) * core->fun(rj, core->param);
}

static double
dblint(double (*f)(double, double, void *), void *param, double a,
       double b, double s, double r)
{
    double dx;
    double dy;
    double sum;
    double x;
    double y;
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;

    n = N_DBL;
    m = N_DBL;
    dx = (b - a) / n;
    dy = (r - s) / m;
    sum = 0;

    for (i = 0; i < n + 1; i++) {
        x = a + dx * i;
        k = (i == 0 || i == n) ? 1 : (i % 2 == 0) ? 2 : 4;
        for (j = 0; j < m + 1; j++) {
            l = (j == 0 || j == m) ? 1 : (j % 2 == 0) ? 2 : 4;
            y = s + dy * j;
            sum += f(x, y, param) * k * l;
        }
    }
    return dx * dy * sum / 9;

}

static double
gauss(double r, void *p0)
{
    double d2;
    struct CoreParam *p;

    p = p0;
    d2 = p->delta * p->delta;
    return exp(-r * r / d2) / (pi * d2);
}
