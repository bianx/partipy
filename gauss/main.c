#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "float/real.h"
#include "double/real.h"

#define pi (3.141592653589793)
static char me[] = "gauss/main";
enum { nmax = 99999 };

static void
usg(void)
{
    fprintf(stderr,
	    "%s -t time -m M -e every -d delta -c [chorin gauss hald j0 krasny] -s [euler rk4] -o [punto|skel|off|gnuplot] [> punto] < initial\n",
	    me);
    exit(1);
}

static int function(int n, const real *, real *, void *);

static real gauss_psi(real, real, void *);
static int gauss_dpsi(real, real, real *, real *, void *);
static real gauss_coef(void *);
static int chorin_dpsi(real, real, real *, real *, void *);
static real chorin_psi(real, real, void *);
static real chorin_coef(void *);
static int hald_dpsi(real, real, real *, real *, void *);
static real hald_coef(void *);
static int krasny_dpsi(real, real, real *, real *, void *);
static real krasny_coef(void *);
static int j0_dpsi(real, real, real *, real *, void *);
static real j0_coef(void *);
static int grid(void * p0, int n, const real * x, const real * y, const real * ksi, int step);

struct Core {
    real (*psi)(real, real, void *);
    int (*dpsi)(real, real, real *, real *, void *);
    real (*coef) (void *);
    void *param;
};

static struct Core Core[] = {
			     { chorin_psi, chorin_dpsi, chorin_coef, NULL },
			     { gauss_psi, gauss_dpsi, gauss_coef, NULL },
			     { NULL, hald_dpsi, hald_coef, NULL },
			     { NULL, j0_dpsi, j0_coef, NULL },
			     { NULL, krasny_dpsi, krasny_coef, NULL },
};

static const char *CoreName[] = {
    "chorin",
    "gauss",
    "hald",
    "j0",
    "krasny",
};

static int punto_write(int n, const real *, const real *, const real *,
		       int step);
static int skel_write(int n, const real *, const real *, const real *,
		      int step);
static int off_write(int n, const real *, const real *, const real *,
		     int step);
static int gnuplot_write(int n, const real *, const real *, const real *,
			 int step);

struct Ode;
struct OdeParam {
    int (*function)(int n, const real *, real *, void *);
    void *param;
    real dt;
    const char *scheme;
};
struct RemeshParam {
  int nx;
  int ny;
  real xlo;
  real xhi;
  real ylo;
  real yhi;
  real eps;
  struct Core *core;
};
static int ode_ini(char **, struct OdeParam *, struct Ode **);
static int ode_step(struct Ode *, int n, real * z);
static int ode_fin(struct Ode *);

static int remesh_m4(void *, int *, real *, real *, real *);
static int remesh_psi(void *, int *, real *, real *, real *);

struct PsiParam {
    real delta;
};

struct Param {
    int n;
    struct Core *core;
    real *ksi;
};

enum { SIZE = 999 };

static const char *WriteName[] = {
    "gnuplot",
    "off",
    "punto",
    "skel",
};

static int (*const WriteFun[])(int, const real *, const real *,
			       const real *, int) = {
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
    const char *scheme;
    int Dflag;
    int Eflag;
    int every;
    int i;
    int j;
    int m;
    int Mflag;
    int n;
    int ncap;
    int nremesh;
    int Sflag;
    int Tflag;
    int (*write)(int, const real *, const real *, const real *, int);
    int (*remesh)(void *, int *, real *, real *, real *);
    real *buf;
    real delta;
    real dt;
    real *ksi;
    real t1;
    real *x;
    real *y;
    real *z;
    struct Core *core;
    struct Ode *ode;
    struct OdeParam ode_param;
    struct Param param;
    struct PsiParam psi_param;
    struct RemeshParam remesh_param;

    Eflag = Dflag = Mflag = Tflag = Sflag = 0;
    write = NULL;
    core = NULL;
    scheme = NULL;
    nremesh = 0;
    remesh = remesh_psi;

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
	case 's':
	    argv++;
	    if (argv[0] == NULL) {
		fprintf(stderr, "%s: -s needs an argument\n", me);
		exit(2);
	    }
	    scheme = argv[0];
	    Sflag = 1;
	    break;
	case 'e':
	    argv++;
	    if (argv[0] == NULL) {
		fprintf(stderr, "%s: -e needs an argument\n", me);
		exit(2);
	    }
	    every = atoi(argv[0]);
	    Eflag = 1;
	    break;
	case 'r':
	    argv++;
	    if (argv[0] == NULL) {
		fprintf(stderr, "%s: -r needs an argument\n", me);
		exit(2);
	    }
	    nremesh = atoi(argv[0]);
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
    if (Mflag == 0) {
	fprintf(stderr, "%s: -m is not given\n", me);
	exit(2);
    }
    if (Sflag == 0) {
	fprintf(stderr, "%s: -s is not given\n", me);
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
    if (Eflag == 0) {
	fprintf(stderr, "%s: -e is not given\n", me);
	exit(2);
    }
    if (write == NULL) {
	fprintf(stderr, "%s: -o is not given\n", me);
	exit(2);
    }
    if (core == NULL) {
	fprintf(stderr, "%s: -c is not given\n", me);
	exit(2);
    }

    ncap = 1;
    if ((buf = malloc(nmax * sizeof(*buf))) == NULL) {
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
    if ((z = malloc(2 * nmax * sizeof(*z))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	exit(2);
    }
    if ((x = malloc(nmax * sizeof(*x))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	exit(2);
    }
    if ((y = malloc(nmax * sizeof(*y))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	exit(2);
    }
    if ((ksi = malloc(nmax * sizeof(*ksi))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	exit(2);
    }
    for (i = j = 0; i < n; i++) {
	x[i] = buf[j++];
	y[i] = buf[j++];
	ksi[i] = buf[j++];
    }

    core->param = &psi_param;
    psi_param.delta = delta;
    param.n = n;
    param.core = core;
    param.ksi = ksi;
    dt = t1 / (m - 1);
    ode_param.function = function;
    ode_param.param = &param;
    ode_param.dt = dt;
    ode_param.scheme = scheme;

    remesh_param.nx = 3.2 / delta * 2;
    remesh_param.ny = 3.2 / delta * 2;
    remesh_param.xlo = -1.6;
    remesh_param.xhi = 1.6;
    remesh_param.ylo = -1.6;
    remesh_param.yhi = 1.6;
    remesh_param.eps = 1e-4;
    remesh_param.core = core;

    if (ode_ini(argv, &ode_param, &ode) != 0) {
	fprintf(stderr, "%s: ode_ini failed\n", me);
	exit(2);
    }
    for (i = 0; ; i++) {
      if (i > 0 && nremesh > 0 && i % nremesh == 0)
	remesh(&remesh_param, &n, x, y, ksi);
      if (every > 0 && i % every == 0) {
	if (write(n, x, y, ksi, i) != 0) {
	  fprintf(stderr, "%s: write failed\n", me);
	  exit(2);
	}
	grid(&remesh_param, n, x, y, ksi, i);
      }
      if (i == m) break;

      for (j = 0; j < n; j++) {
	z[j] = x[j];
	z[j + n] = y[j];
      }
      if (ode_step(ode, 2 * n, z) != 0) {
	fprintf(stderr, "%s: ode_step failed\n", me);
	exit(2);
      }
      for (j = 0; j < n; j++) {
	x[j] = z[j];
	y[j] = z[j + n];
      }
    }
    free(z);
    free(x);
    free(y);
    free(buf);
    free(ksi);
    ode_fin(ode);
}

static int
function(int n, const real * z, real * f, void *params0)
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
    real coef;
    int i;
    int j;
    struct Param *params;
    struct Core *core;

    params = params0;
    core = params->core;
    ksi = params->ksi;
    coef = core->coef(core->param);
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
		core->dpsi(dx, dy, &gx, &gy, core->param);
		fx[i] -= ksi[j] * gy;
		fy[i] += ksi[j] * gx;
	    }
    for (i = 0; i < n; i++) {
	fx[i] *= coef;
	fy[i] *= coef;
    }
    return 0;
}

static int
punto_write(int n, const real * x, const real * y, const real * ksi,
	    int step)
{
    int j;

    if (step > 0)
	printf("\n");
    for (j = 0; j < n; j++)
	printf("%.16e %.16e %.16e\n", x[j], y[j], ksi[j]);
    return 0;
}

static int
skel_write(int n, const real * x, const real * y, const real * ksi,
	   int step)
{
    char path[SIZE];
    real z;
    FILE *f;
    int i;
    int npolylines;

    (void) ksi;

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
off_write(int n, const real * x, const real * y, const real * ksi,
	  int step)
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
    static const real r = 0.05;

    (void) ksi;

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
gnuplot_write(int n, const real * x, const real * y, const real * ksi,
	      int step)
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
      fprintf(f, "%.16g %.16g %.16g\n", x[i], y[i], ksi[i]);
    if (fclose(f) != 0) {
	fprintf(stderr, "%s: fail to close '%s'\n", me, path);
	exit(2);
    }
    return 0;
}

static real
hald_coef(void *p)
{
    (void) p;
    return 1 / (2 * pi);
}

static int
hald_dpsi(real x, real y, real * u, real * v, void *p0)
{
    real coef;
    real delta;
    real r2;
    real r;
    struct PsiParam *p;

    p = p0;
    delta = p->delta;
    r2 = x * x + y * y;
    if (r2 > 10 * DBL_MIN) {
	r = sqrtr(r2) / delta;
	coef = 1 - 8.0 / 3.0 * j1(2 * r) / (2 * r) + 2.0 / 3 * j1(r) / r;
	coef /= r2;
	*u = coef * x;
	*v = coef * y;
    } else {
	*u = 0;
	*v = 0;
    }
    return 0;
}

static real
gauss_coef(void *p)
{
    (void) p;
    return 1 / (2 * pi);
}

static real
gauss_psi(real x, real y, void *p0)
{
    real d2;
    real r2;
    struct PsiParam *p;

    p = p0;
    d2 = p->delta * p->delta;
    r2 = x * x + y * y;
    return exp(-r2 / d2) / (pi * d2);
}


static int
gauss_dpsi(real x, real y, real * u, real * v, void *p0)
{
    real coef;
    real d2;
    real delta;
    real r2;
    struct PsiParam *p;

    p = p0;
    delta = p->delta;
    d2 = delta * delta;
    r2 = x * x + y * y;
    if (r2 > 10 * DBL_MIN) {
	coef = (1 - exp(-r2 / d2)) / r2;
	*u = coef * x;
	*v = coef * y;
    } else {
	*u = 0;
	*v = 0;
    }
    return 0;
}

static real
j0_coef(void *p)
{
    (void) p;
    return 1 / (2 * pi);
}

static int
j0_dpsi(real x, real y, real * u, real * v, void *p0)
{
    real coef;
    real delta;
    real r2;
    real r;
    struct PsiParam *p;

    p = p0;
    delta = p->delta;
    r2 = x * x + y * y;
    if (r2 > 10 * DBL_MIN) {
	r = sqrtr(r2);
	coef = (1 - j0(2 * r / delta)) / r2;
	*u = coef * x;
	*v = coef * y;
    } else {
	*u = 0;
	*v = 0;
    }
    return 0;
}

static real
chorin_psi(real x, real y, void *p0)
{
    real delta;
    real d3;
    real r;
    struct PsiParam *p;

    p = p0;
    delta = p->delta;
    d3 = delta * delta * delta;
    r = sqrtr(x * x + y * y);
    return r < delta ? 3 * r / (2 * pi * d3)  : 0;
}

static real
chorin_coef(void *p)
{
    (void) p;
    return 1 / (2 * pi);
}

static int
chorin_dpsi(real x, real y, real * u, real * v, void *p0)
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


static real
krasny_coef(void *p)
{
    (void) p;
    return 1.0;
}

static int
krasny_dpsi(real x, real y, real * u, real * v, void *p0)
{
    real delta;
    real den;
    struct PsiParam *p;

    p = p0;
    delta = p->delta;
    den = cosh(2 * pi * y) - cos(2 * pi * x) + delta * delta;
    if (fabs(den) < 10 * DBL_MIN) {
	fprintf(stderr, "den is too small\n");
	return 1;
    }
    *u = sin(2 * pi * x) / den;
    *v = sinh(2 * pi * y) / den;
    return 0;
}

struct Ode {
    int (*function)(int n, const real *, real *, void *);
    void *param;
    real dt;
    real *k;

    real *y0;
    real *ytmp;
    int (*step)(struct Ode *, int n, real *);
};

static int step_euler(struct Ode *, int n, real *);
static int step_rk4(struct Ode *, int n, real *);

static const char *OdeName[] = {
    "euler",
    "rk4",
};

static int (*const OdeStep[])(struct Ode *, int n, real *) = {
    step_euler,
    step_rk4,
};

static int
ode_ini(char **argv, struct OdeParam *p, struct Ode **pq)
{
    real *k;
    real *y0;
    real *ytmp;
    int i;
    struct Ode *q;

    (void) argv;
    if ((q = malloc(sizeof(*q))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	return 1;
    }
    if ((k = malloc(nmax * sizeof(*k))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	return 1;
    }
    if ((y0 = malloc(nmax * sizeof(*y0))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	return 1;
    }
    if ((ytmp = malloc(nmax * sizeof(*ytmp))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	return 1;
    }
    for (i = 0;; i++) {
	if (i == sizeof(OdeName) / sizeof(*OdeName)) {
	    fprintf(stderr, "%s: unknown scheme type '%s'\n", me,
		    p->scheme);
	    goto err;
	}
	if (strncmp(p->scheme, OdeName[i], SIZE) == 0) {
	    q->step = OdeStep[i];
	    break;
	}
    }
    q->dt = p->dt;
    q->k = k;
    q->function = p->function;
    q->param = p->param;
    q->y0 = y0;
    q->ytmp = ytmp;

    *pq = q;
    return 0;
  err:
    return 1;
}

static int
ode_step(struct Ode *q, int n, real * y)
{
  return q->step(q, n, y);
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
step_euler(struct Ode *q, int n, real * y)
{
    real *k;
    real dt;
    int i;

    k = q->k;
    dt = q->dt;
    if (q->function(n/2, y, k, q->param) != 0) {
	fprintf(stderr, "%s:%d: function failed\n", __FILE__, __LINE__);
	return 1;
    }
    for (i = 0; i < n; i++)
	y[i] += dt * k[i];
    return 0;
}

static int
step_rk4(struct Ode *q, int n, real * y)
{
    real h;
    real *k;
    real *y0;
    real *ytmp;
    int i;

#define EVAL(y, k) \
    if (q->function((n)/2, (y), (k), q->param) != 0) {			     \
    fprintf(stderr, "%s:%d: function failed\n", __FILE__, __LINE__); \
    return 1; \
  } \

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

static real
m4(real cutoff, real x)
{
  real a;
  real b;
  x = fabs(x)/cutoff;
  a = 2 - x;
  b = 1 - x;
  if (x < 1)
    return (a*a*a/6 - 4*b*b*b/6) / cutoff;
  else if (x < 2)
    return a*a*a/6/ cutoff;
  else
    return 0;
}

static int
remesh_m4(void * p0, int * pn, real * x, real * y, real * ksi) {
  int i;
  int j;
  int k;
  int l;
  int n;
  int m;
  int nx;
  int ny;
  real coef;
  real dx;
  real dy;
  real ksi0[99999];
  real u;
  real v;
  real xhi;
  real xlo;
  real yhi;
  real ylo;
  real eps;
  struct RemeshParam *p;

  n = *pn;
  p = p0;
  nx = p->nx;
  ny = p->ny;
  xlo = p->xlo;
  xhi = p->xhi;
  ylo = p->ylo;
  yhi = p->yhi;
  eps = p->eps;

  dx = (xhi - xlo) / nx;
  dy = (yhi - ylo) / ny;
  m = nx * ny;
  for (i = 0; i < m; i++)
    ksi0[i] = 0;

  for (k = 0; k < n; k++) {
    l = 0;
    for (i = 0; i < nx; i++) {
      u = xlo + (i + 0.5) * dx;
      for (j = 0; j < ny; j++) {
	v = ylo + (j + 0.5) * dy;
	ksi0[l] += ksi[k] * m4(dx, x[k] - u) * m4(dy, y[k] - v);
	l++;
      }
    }
  }
  coef = dx * dy;
  for (i = 0; i < m; i++)
    ksi[i] = ksi0[i] * coef;

  l = 0;
  for (i = 0; i < nx; i++) {
      u = xlo + (i + 0.5) * dx;
      for (j = 0; j < ny; j++) {
	v = ylo + (j + 0.5) * dy;
	x[l] = u;
	y[l] = v;
	l++;
      }
  }

  j = 0;
  for (i = 0; i < m; i++)
    if (fabs(ksi[i]) > eps) {
      x[j] = x[i];
      y[j] = y[i];
      ksi[j] = ksi[i];
      j++;
    }
  *pn = j;
  return 0;
}


static int
remesh_psi(void * p0, int * pn, real * x, real * y, real * ksi) {
  int i;
  int j;
  int k;
  int l;
  int n;
  int m;
  int nx;
  int ny;
  real coef;
  real dx;
  real dy;
  real ksi0[999999];
  real u;
  real v;
  real xhi;
  real xlo;
  real yhi;
  real ylo;
  real eps;
  struct Core *core;
  struct RemeshParam *p;

  n = *pn;
  p = p0;
  nx = p->nx;
  ny = p->ny;
  xlo = p->xlo;
  xhi = p->xhi;
  ylo = p->ylo;
  yhi = p->yhi;
  core = p->core;
  eps = p->eps;

  dx = (xhi - xlo) / nx;
  dy = (yhi - ylo) / ny;
  m = nx * ny;
  for (i = 0; i < m; i++)
    ksi0[i] = 0;

  for (k = 0; k < n; k++) {
    l = 0;
    for (i = 0; i < nx; i++) {
      u = xlo + (i + 0.5) * dx;
      for (j = 0; j < ny; j++) {
	v = ylo + (j + 0.5) * dy;
	ksi0[l] += ksi[k] * core->psi(x[k] - u, y[k] - v, core->param);
	l++;
      }
    }
  }
  coef = dx * dy;
  for (i = 0; i < m; i++)
    ksi[i] = ksi0[i] * coef;

  l = 0;
  for (i = 0; i < nx; i++) {
      u = xlo + (i + 0.5) * dx;
      for (j = 0; j < ny; j++) {
	v = ylo + (j + 0.5) * dy;
	x[l] = u;
	y[l] = v;
	l++;
      }
  }

  j = 0;
  for (i = 0; i < m; i++)
    if (fabs(ksi[i]) > eps) {
      x[j] = x[i];
      y[j] = y[i];
      ksi[j] = ksi[i];
      j++;
    }
  *pn = j;
  return 0;
}

static int
grid(void * p0, int n, const real * x, const real * y, const real * ksi, int step) {
  FILE *f;
  int i;
  int j;
  int k;
  int l;
  int m;
  int nx;
  int ny;
  real coef;
  real dx;
  real dy;
  real ksi0[99999];
  real u;
  real v;
  real xhi;
  real xlo;
  real yhi;
  real ylo;
  int status;
  char path[SIZE];
  struct Core *core;
  struct RemeshParam *p;

  p = p0;
  nx = p->nx;
  ny = p->ny;
  xlo = p->xlo;
  xhi = p->xhi;
  ylo = p->ylo;
  yhi = p->yhi;
  core = p->core;

  if (core->psi == NULL)
    return 0;

  dx = (xhi - xlo) / nx;
  dy = (yhi - ylo) / ny;
  m = nx * ny;
  for (i = 0; i < m; i++)
    ksi0[i] = 0;

  for (k = 0; k < n; k++) {
    l = 0;
    for (j = 0; j < ny; j++) {
	v = ylo + (j + 0.5) * dy;
	for (i = 0; i < nx; i++) {
	  u = xlo + (i + 0.5) * dx;
	  ksi0[l] += ksi[k] * gauss_psi(x[k] - u, y[k] - v, core->param);
	  l++;
	}
    }
  }
  coef = dx * dy;
  for (i = 0; i < m; i++)
    ksi0[i] *= coef;
  snprintf(path, SIZE, "a.%06d.vtk", step);
  fprintf(stderr, "%s: write '%s'\n", me, path);
  if ((f = fopen(path, "w")) == NULL) {
    fprintf(stderr, "%s: fail to open '%s'\n", me, path);
    exit(2);
  }
  status = fprintf(f, "# vtk DataFile Version 2.0\n"
		   "generated by %s\n"
		   "ASCII\n"
		   "DATASET STRUCTURED_POINTS\n"
		   "DIMENSIONS %d %d 1\n"
		   "ORIGIN %.16e %.16e 0\n"
		   "SPACING %.16e %.16e 0\n",
		   me, nx, ny, xlo + dx/2, ylo + dy/2,
		   dx, dy);
  if (status < 0) {
    fprintf(stderr, "%s: fail to write '%s'\n", me, path);
    exit(2);
  }
  fprintf(f, "POINT_DATA %d\n"
	  "SCALARS omega double\n" "LOOKUP_TABLE DEFAULT\n", m);
  for (i = 0; i < m; i++)
    fprintf(f, "%.16e\n", ksi0[i]);
  fclose(f);
  return 0;
}
