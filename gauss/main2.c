#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
//#include "float/real.h"
#include "double/real.h"

#define pi (3.141592653589793)
static char me[] = "gauss/main2";

static void
usg(void)
{
    fprintf(stderr,
            "%s -t time -m M -d delta -s [rk2 rk4 rk8pd rkck rkf45] -o [punto|skel|off|gnuplot] [> punto] < initial\n",
            me);
    exit(1);
}

static int function(const real *, real *, void *);

static int dpsi(real, real, real *, real *, void *);
static const real dpsi_coef = 1 / (2 * pi);

static int punto_write(int n, const real *, const real *, int step);
static int skel_write(int n, const real *, const real *, int step);
static int off_write(int n, const real *, const real *, int step);
static int gnuplot_write(int n, const real *, const real *, int step);

struct Ode;
struct OdeParam {
    int n;
    int (*function)(const real *, real *, void *);
    void *param;
    real dt;
};
static int ode_ini(char **, struct OdeParam *, struct Ode **);
static int ode_step(struct Ode *, real * y);
static int ode_fin(struct Ode *);

struct PsiParam {
    real delta;
};

struct Param {
    int n;
    void *psi_param;
    real *ksi;
};

enum { SIZE = 999 };

static const char *WriteName[] = {
    "gnuplot",
    "off",
    "punto",
    "skel",
};

static int (*const WriteFun[])(int, const real *, const real *, int) = {
    gnuplot_write,
    off_write,
    punto_write,
    skel_write,
};

int
main(int argc, char **argv)
{
    (void) argc;
    char line[SIZE];
    real dt;
    real delta;
    real t1;
    real *x;
    real *y;
    real *z;
    real *buf;
    real *ksi;
    int Dflag;
    int i;
    int j;
    int m;
    int Mflag;
    int n;
    int ncap;
    int Tflag;
    int (*write)(int, const real *, const real *, int);
    struct Param param;
    struct PsiParam psi_param;
    struct Ode *ode;
    struct OdeParam ode_param;

    Dflag = Mflag = Tflag = 0;
    write = NULL;
    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        case 'm':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -m needs an argument\n", me);
                exit(2);
            }
            m = atoi(argv[0]);
            Mflag = 1;
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
        case 't':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -t needs an argument\n", me);
                exit(2);
            }
            t1 = atof(argv[0]);
            Tflag = 1;
            break;
        case 'o':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -o needs an argument\n", me);
                exit(2);
            }
            for (i = 0;; i++) {
                if (i == sizeof(WriteFun) / sizeof(*WriteFun)) {
                    fprintf(stderr, "%s: unknown output type '%s'\n", me,
                            argv[0]);
                    exit(2);
                }
                if (strncmp(argv[0], WriteName[i], SIZE) == 0) {
                    write = WriteFun[i];
                    break;
                }
            }
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }
    if (Mflag == 0) {
        fprintf(stderr, "%s: -m is not given\n", me);
        exit(2);
    }
    if (Tflag == 0) {
        fprintf(stderr, "%s: -t is not given\n", me);
        exit(2);
    }
    if (Dflag == 0) {
        fprintf(stderr, "%s: -d is not given\n", me);
        exit(2);
    }
    if (write == NULL) {
        fprintf(stderr, "%s: -o is not given\n", me);
        exit(2);
    }

    ncap = 1;
    if ((buf = malloc(ncap * sizeof(*buf))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }
    n = 0;
    while (fgets(line, SIZE, stdin) != NULL) {
        while (ncap <= 3 * n + 2) {
            ncap *= 2;
            if ((buf = realloc(buf, ncap * sizeof(*buf))) == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(2);
            }
        }
        if (sscanf
            (line, "%" FM " %" FM " %" FM, &buf[3 * n], &buf[3 * n + 1],
             &buf[3 * n + 2]) != 3) {
            fprintf(stderr, "%s: fail to parse '%s'\n", me, line);
            exit(2);
        }
        n++;
    }
    if ((z = malloc(2 * n * sizeof(*z))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }
    if ((ksi = malloc(n * sizeof(*ksi))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }
    x = z;
    y = &z[n];
    for (i = j = 0; i < n; i++) {
        x[i] = buf[j++];
        y[i] = buf[j++];
        ksi[i] = buf[j++];
    }

    psi_param.delta = delta;
    param.n = n;
    param.psi_param = &psi_param;
    param.ksi = ksi;
    dt = t1 / (m - 1);
    ode_param.n = 2 * n;
    ode_param.function = function;
    ode_param.param = &param;
    ode_param.dt = dt;
    if (ode_ini(argv, &ode_param, &ode) != 0) {
        fprintf(stderr, "%s: ode_ini failed\n", me);
        exit(2);
    }
    write(n, x, y, 0);
    fprintf(stderr, "delta: %g\n", delta);
    for (i = 1; i < m; i++) {
        if (ode_step(ode, z) != 0) {
            fprintf(stderr, "%s: ode_step failed\n", me);
            exit(2);
        }
	write(n, x, y, i);
    }
    free(z);
    free(buf);
    free(ksi);
    ode_fin(ode);
}

static int
function(const real * z, real * f, void *params0)
{
    const real *ksi;
    const real *x;
    const real *y;
    real dx;
    real dy;
    real *fx;
    real *fy;
    real gx;
    real gy;
    int i;
    int j;
    int n;
    struct Param *params;
    void *psi_param;

    params = params0;
    n = params->n;
    psi_param = params->psi_param;
    ksi = params->ksi;
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
                dpsi(dx, dy, &gx, &gy, psi_param);
                fx[i] -= ksi[j] * gy;
                fy[i] += ksi[j] * gx;
            }
    for (i = 0; i < n; i++) {
        fx[i] *= dpsi_coef;
        fy[i] *= dpsi_coef;
    }
    return 0;
}

static int
punto_write(int n, const real * x, const real * y, int step)
{
    int j;

    if (step > 0)
        printf("\n");
    for (j = 0; j < n; j++)
        printf("%.16e %.16e\n", x[j], y[j]);
    return 0;
}

static int
skel_write(int n, const real * x, const real * y, int step)
{
    char path[SIZE];
    real z;
    FILE *f;
    int i;
    int npolylines;

    npolylines = 1;
    z = 0;
    snprintf(path, SIZE, "%06d.skel", step);
    fprintf(stderr, "%s: write '%s'\n", me, path);
    if ((f = fopen(path, "w")) == NULL) {
        fprintf(stderr, "%s: fail to open '%s'\n", me, path);
        exit(2);
    }
    if (fputs("SKEL\n", f) == EOF) {
        fprintf(stderr, "%s: fail to write '%s'\n", me, path);
        exit(2);
    }
    fprintf(f, "%d %d\n", n, npolylines);
    for (i = 0; i < n; i++)
        fprintf(f, "%.16g %.16g %.16g\n", x[i], y[i], z);
    fprintf(f, "%d", n);
    for (i = 0; i < n; i++)
        fprintf(f, " %d", i);
    fprintf(f, "\n");

    if (fclose(f) != 0) {
        fprintf(stderr, "%s: fail to close '%s'\n", me, path);
        exit(2);
    }
    return 0;
}

static int
off_write(int n, const real * x, const real * y, int step)
{
#define NTRI (50)
    char path[SIZE];
    real u[NTRI];
    real v[NTRI];
    real z;
    real h;
    FILE *f;
    int i;
    int j;
    int k;
    int m;
    static const real r = 0.0005;

    m = NTRI;
    h = 2 * pi / m;
    for (i = 0; i < m; i++) {
        u[i] = r * cosr(i * h);
        v[i] = r * sinr(i * h);
    }
    z = 0;
    snprintf(path, SIZE, "%06d.off", step);
    fprintf(stderr, "%s: write '%s'\n", me, path);
    if ((f = fopen(path, "w")) == NULL) {
        fprintf(stderr, "%s: fail to open '%s'\n", me, path);
        exit(2);
    }
    if (fputs("OFF\n", f) == EOF) {
        fprintf(stderr, "%s: fail to write '%s'\n", me, path);
        exit(2);
    }
    fprintf(f, "%d %d 0\n", (1 + m) * n, m * n);
    for (i = 0; i < n; i++) {
        fprintf(f, "%.16g %.16g %.16g\n", x[i], y[i], z);
        for (j = 0; j < m; j++)
            fprintf(f, "%.16g %.16g %.16g\n", x[i] + u[j], y[i] + v[j], z);
    }
    for (i = 0; i < n; i++) {
        k = (m + 1) * i;
        for (j = 0; j < m - 1; j++)
            fprintf(f, "3 %d %d %d\n", k, k + j + 1, k + j + 2);
        fprintf(f, "3 %d %d %d\n", k, k + m, k + 1);
    }
    if (fclose(f) != 0) {
        fprintf(stderr, "%s: fail to close '%s'\n", me, path);
        exit(2);
    }
    return 0;
}

static int
gnuplot_write(int n, const real * x, const real * y, int step)
{
    char path[SIZE];
    int i;
    FILE *f;

    snprintf(path, SIZE, "%06d.dat", step);
    fprintf(stderr, "%s: write '%s'\n", me, path);
    if ((f = fopen(path, "w")) == NULL) {
        fprintf(stderr, "%s: fail to open '%s'\n", me, path);
        exit(2);
    }
    for (i = 0; i < n; i++)
        fprintf(f, "%.16g %.16g\n", x[i], y[i]);
    if (fclose(f) != 0) {
        fprintf(stderr, "%s: fail to close '%s'\n", me, path);
        exit(2);
    }
    return 0;
}

static int
dpsi(real x, real y, real * u, real * v, void *p0)
{
    real delta;
    real r2;
    real coef;
    struct PsiParam *p;

    p = p0;
    delta = p->delta;
    r2 = x * x + y * y;
    if (r2 > delta * delta) {
        coef = 1 / r2;
        *u = coef * x;
        *v = coef * y;
    } else if (r2 > 10 * DBL_MIN) {
        coef = 1 / sqrtr(r2) / delta;
        *u = coef * x;
        *v = coef * y;
    } else {
        *u = 0;
        *v = 0;
    }
    return 0;
}

struct Ode {
    int n;
    int (*function)(const real *, real *, void *);
    void *param;
    real dt;
    real *k;

    real *y0;
    real *ytmp;
    int (*step)(struct Ode *, real *);
};


static int step_euler(struct Ode *, real *);
static int step_rk4(struct Ode *, real *);

static int
ode_ini(char **argv, struct OdeParam *p, struct Ode **pq)
{
    real *k;
    real *y0;
    real *ytmp;
    struct Ode *q;

    if ((q = malloc(sizeof(*q))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((k = malloc(p->n * sizeof(*k))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((y0 = malloc(p->n * sizeof(*y0))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((ytmp = malloc(p->n * sizeof(*ytmp))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    q->dt = p->dt;
    q->k = k;
    q->function = p->function;
    q->n = p->n;
    q->param = p->param;
    q->step = step_euler;
    q->y0 = y0;
    q->ytmp = ytmp;

    *pq = q;
    return 0;
}

static int
ode_step(struct Ode *q, real * y)
{
    return q->step(q, y);
}

static int
ode_fin(struct Ode *q)
{
    free(q->k);
    free(q->y0);
    free(q->ytmp);
    free(q);
    return 0;
}

static int
step_euler(struct Ode *q, real * y)
{
    real *k;
    real dt;
    int i;

    k = q->k;
    dt = q->dt;
    if (q->function(y, k, q->param) != 0) {
        fprintf(stderr, "%s:%d: function failed\n", __FILE__, __LINE__);
        return 1;
    }
    for (i = 0; i < q->n; i++)
        y[i] += dt * k[i];
    return 0;
}

static int
step_rk4(struct Ode *q, real * y)
{
    real h;
    real *k;
    real *y0;
    real *ytmp;
    int i;
    int n;

#define EVAL(y, k) \
  if (q->function((y), (k), q->param) != 0) {			     \
    fprintf(stderr, "%s:%d: function failed\n", __FILE__, __LINE__); \
    return 1; \
  } \

    n = q->n;
    k = q->k;
    y0 = q->y0;
    ytmp = q->ytmp;
    h = q->dt;
    for (i = 0; i < n; i++)
        y0[i] = y[i];
    EVAL(y, k);
    for (i = 0; i < n; i++) {
        y[i] += h * k[i] / 6;
        ytmp[i] = y0[i] + h * k[i] / 2;
    }
    EVAL(ytmp, k);
    for (i = 0; i < n; i++) {
        y[i] += h * k[i] / 3;
        ytmp[i] = y0[i] + h * k[i] / 2;
    }
    EVAL(ytmp, k);
    for (i = 0; i < n; i++) {
        y[i] += h * k[i] / 3;
        ytmp[i] = y0[i] + h * k[i] / 2;
    }
    EVAL(ytmp, k);
    for (i = 0; i < n; i++) {
        y[i] += h * k[i] / 6;
    }
    return 0;
}
