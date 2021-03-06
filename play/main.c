#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

static char me[] = "play/main";
static void
usg(void)
{
    fprintf(stderr, "%s -n int > out\n", me);
    exit(1);
}

enum { M = 20 };
static const double pi = 3.141592653589793;
static const double L = 1.0;
static const double t1 = 4.0;
static const double epsrel = 1e-6;
static const double epsabs = 0;
double Ham[M + 1];
double MomX[M + 1];
double MomY[M + 1];
double MomAz[M + 1];
int iH = 0;


struct Param {
    int n;
    double eps;
    double dt;
};

int
func(double t, const double *z, double *f, void *params0)
{
    (void) (t);
    const double *x;
    const double *y;
    double *fx;
    double *fy;
    double dx;
    double dy;
    double den;
    double eps;
    int i;
    int j;
    int n;
    struct Param *params;

    params = params0;
    n = params->n;
    eps = params->eps;
    x = z;
    y = &z[n];
    fx = f;
    fy = &f[n];

    for (i = 0; i < n; i++)
        fx[i] = fy[i] = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            if (i != j) {
                dx = x[i] - x[j];
                dy = y[i] - y[j];
                den = cosh(2 * pi * dy) - cos(2 * pi * dx) + eps * eps;
                fx[i] -= sinh(2 * pi * dy) / den;
                fy[i] += sin(2 * pi * dx) / den;
            }
    }

    for (i = 0; i < n; i++) {
        fx[i] /= 2 * n;
        fy[i] /= 2 * n;
    }
    return GSL_SUCCESS;
}

int
post(const double *z, void *params0)
{
    const double *x;
    const double *y;
    double dx;
    double dy;
    double den;
    double eps;
    int i;
    int j;
    int n;
    struct Param *params;

    double Ht;
    double mx;
    double my;

    params = params0;
    n = params->n;
    eps = params->eps;
    x = z;
    y = &z[n];

    Ht = 0;
    for (i = 0; i < n; i++) {
        mx = 0;
        my = 0;
        for (j = 0; j < n; j++)
            if (i != j) {
                dx = x[i] - x[j];
                dy = y[i] - y[j];
                den = cosh(2 * pi * dy) - cos(2 * pi * dx) + eps * eps;
                mx -= sinh(2 * pi * dy) / den;
                my += sin(2 * pi * dx) / den;
                if (j > i) {
                    Ht += log(den);
                }
            }
        MomX[iH] += mx;
        MomY[iH] += my;
        MomAz[iH] += x[i] * my - y[i] * mx;

    }
    Ham[iH++] = -Ht / 4 / pi / n / n;

    return 0;
}

int
main(int argc, char **argv)
{
    (void) argc;
    int i;
    int j;
    double t;
    int n;
    double *z;
    double *x;
    double *y;
    double h;
    double ti;
    double dt, dt_out;
    double x0;
    double eps;
    gsl_odeiv2_driver *driver;
    double dt_start;
    gsl_odeiv2_system sys;
    struct Param param;
    int Nflag;
    int Eflag;
    int Dtflag;

    Nflag = Eflag = Dtflag = 0;
    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        case 'n':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -n needs an argument\n", me);
                exit(2);
            }
            n = atoi(argv[0]);
            Nflag = 1;
            break;
        case 'e':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -e needs an argument\n", me);
                exit(2);
            }
            eps = atof(argv[0]);
            Eflag = 1;
            break;
        case 'd':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -d needs an argument\n", me);
                exit(2);
            }
            dt = atof(argv[0]);
            Dtflag = 1;
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }
    if (Nflag == 0) {
        fprintf(stderr, "%s: -n is not given\n", me);
        exit(2);
    }
    if (Eflag == 0) {
        fprintf(stderr, "%s: -e is not given\n", me);
        exit(2);
    }
    if (Dtflag == 0) {
        fprintf(stderr, "%s: -d is not given\n", me);
        exit(2);
    }

    z = malloc(2 * n * sizeof(*z));
    if (z == NULL) {
        fprintf(stderr, "%s: malloc failed\n", me);
        exit(2);
    }
    param.n = n;
    param.eps = eps;
    param.dt = dt;
    //dt_start = 0.0001;
    dt_start = dt;
    sys.function = func;
    sys.jacobian = NULL;
    sys.dimension = 2 * n;
    sys.params = &param;
    driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4,
                                           dt_start, epsrel, epsabs);
    if (driver == NULL) {
        fprintf(stderr, "%s: driver allocation failed\n", me);
        exit(2);
    }
    x = z;
    y = &z[n];
    h = L / n;
    for (i = 0; i < n; i++) {
        x0 = i * h;
        x[i] = x0 + 0.01 * sin(2 * pi * x0);
        y[i] = -0.01 * sin(2 * pi * x0);
    }
    dt_out = t1 / M;
    t = 0;
    for (i = 0; i <= M; i++) {
        ti = dt_out * i;
        fprintf(stderr, "%s: %g\n", me, ti);

        if (post(z, &param) != 0) {
            fprintf(stderr, "%s: post failed\n", me);
            exit(2);
        }

        if (gsl_odeiv2_driver_apply(driver, &t, ti, z) != GSL_SUCCESS) {
            fprintf(stderr, "%s: driver failed\n", me);
            exit(2);
        }

        if (i > 0)
            printf("\n");
        //printf ("#t=%.5g\n", ti);
        for (j = 0; j < n; j++)
            printf("%.16g %.16g\n", x[j], y[j]);
    }

    printf("#This is Ham, MomX, MomY, MomAz\n");
    for (i = 0; i <= M; i++) {
        ti = dt_out * i;
        printf("%.5g %.16g %.16g %.16g %.16g\n", ti, Ham[i], MomX[i],
               MomY[i], MomAz[i]);
    }
    free(z);
    gsl_odeiv2_driver_free(driver);
}
