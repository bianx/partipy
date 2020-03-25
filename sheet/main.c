#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

static char me[] = "sheet/main";
static void
usg(void)
{
    fprintf(stderr, "%s -n N -e eps -s [rk2 rk4 rk8pd rkck rkf45] > out\n", me);
    exit(1);
}

struct Param {
    int n;
    double eps;
};
enum { SIZE = 999 };
enum { M = 10 };
static const double pi = 3.141592653589793;
static const double L = 1.0;
static const double t1 = 4.0;
static const double epsrel = 1e-6;
static const double epsabs = 0;

static const gsl_odeiv2_step_type **Type[] = {
    &gsl_odeiv2_step_rk2,
    &gsl_odeiv2_step_rk4,
    &gsl_odeiv2_step_rk8pd,
    &gsl_odeiv2_step_rkck,
    &gsl_odeiv2_step_rkf45,
};

static const char *Name[] = {
    "rk2",
    "rk4",
    "rk8pd",
    "rkck",
    "rkf45",
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
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (i != j) {
                dx = x[i] - x[j];
                dy = y[i] - y[j];
                den = cosh(2 * pi * dy) - cos(2 * pi * dx) + eps * eps;
                fx[i] -= sinh(2 * pi * dy) / den;
                fy[i] += sin(2 * pi * dx) / den;
            }
    for (i = 0; i < n; i++) {
        fx[i] /= 2 * n;
        fy[i] /= 2 * n;
    }
    return GSL_SUCCESS;
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
    double dt;
    double x0;
    double eps;
    gsl_odeiv2_driver *driver;
    double dt_start;
    gsl_odeiv2_system sys;
    struct Param param;
    int Nflag;
    int Eflag;
    const char *scheme = "rk4";
    Nflag = Eflag = 0;
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
        case 's':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -s needs an argument\n", me);
                exit(2);
            }
            scheme = argv[0];
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
    z = malloc(2 * n * sizeof(*z));
    if (z == NULL) {
        fprintf(stderr, "%s: malloc failed\n", me);
        exit(2);
    }
    param.n = n;
    param.eps = eps;
    dt_start = 0.001;
    sys.function = func;
    sys.jacobian = NULL;
    sys.dimension = 2 * n;
    sys.params = &param;

    for (i = 0; ; i++) {
      if (i == sizeof(Type)/sizeof(*Type)) {
	fprintf(stderr, "%s: unknown scheme '%s'\n", me, scheme);
	exit(2);
      }
      if (strncmp(scheme, Name[i], SIZE) == 0) {
	driver = gsl_odeiv2_driver_alloc_y_new(&sys, *Type[i],
					       dt_start, epsrel, epsabs);
	if (driver == NULL) {
	  fprintf(stderr, "%s: driver allocation failed\n", me);
	  exit(2);
	}
	break;
      }
    }
    x = z;
    y = &z[n];
    h = L / (n - 1);
    for (i = 0; i < n; i++) {
        x0 = i * h;
        x[i] = x0 + 0.01 * sin(2 * pi * x0);
        y[i] = -0.01 * sin(2 * pi * x0);
    }
    dt = t1 / (M - 1);
    t = 0;
    for (i = 0; i < M; i++) {
        ti = dt * i;
        if (gsl_odeiv2_driver_apply(driver, &t, ti, z) != GSL_SUCCESS) {
            fprintf(stderr, "%s: driver failed\n", me);
            exit(2);
        }
        if (i > 0)
            printf("\n");
        for (j = 0; j < n; j++)
            printf("%.16e %.16e\n", x[j], y[j]);
    }
    free(z);
    gsl_odeiv2_driver_free(driver);
}
