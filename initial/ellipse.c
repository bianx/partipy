#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>

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
static double lin(double r, double phi, void *p0);
static double gauss(double);
static double chorin(double);
static double hald(double);
static double j2(double);

static const double q = 2.56085;
static const double Ksi = 20.0;
static const double a = 0.8;
static const double b = 2 * 0.8;
static const double delta = 0.4;

struct CrossParam {
  double xi;
  double yi;
  double xj;
  double yj;
  double (*psi)(double);  
};

struct LinParam {
  double x;
  double y;
  double (*vor)(double, double);
  double (*psi)(double);
};

static double f(double z, double q)
{
  return exp(-(q/z)*exp(1/(z - 1)));
}

static double vorI(double x, double y)
{
  double r;
  x /= a;
  y /= b;
  r = sqrt(x*x + y*y);
  if (r < 1)
    return Ksi * (1 - f(r, q));
  else
    return 0;
}

static double vorConst(double x, double y)
{
  return Ksi;
}

int
main(int argc, char **argv)
{
    char line[SIZE];
    double *A;
    double *B;
    double *C;
    double gamma;
    double *invC;
    double *Ksi;
    double *x;
    double *y;
    gsl_matrix_view vA;
    gsl_matrix_view vC;
    gsl_matrix_view vInvC;
    gsl_permutation *p;
    gsl_vector_view vB;
    gsl_vector_view vKsi;
    int i;
    int j;
    int k;
    int n;
    int ncap;
    int s;
    int Verbose;
    struct CrossParam cross_param;
    struct LinParam lin_param;
    (void) argc;
    
    gamma = 0;
    lin_param.vor = vorConst;
    lin_param.psi = cross_param.psi = hald;
    Verbose = getenv("LOG") != NULL;
    while (*++argv != NULL && argv[0][0] == '-')
	switch (argv[0][1]) {
	case 'h':
	    usg();
	    break;
	default:
	    fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
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
    if ((C = malloc(n * n * sizeof(*C))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	exit(2);
    }
    if ((invC = malloc(n * n * sizeof(*invC))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	exit(2);
    }
    if ((Ksi = malloc(n * sizeof(*Ksi))) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	exit(2);
    }
    for (i = 0; i < n; i++) {
      if (Verbose)
	fprintf(stderr, "%s: %04d of %04d\n", me, i, n);
      for (j = i; j < n; j++) {
	cross_param.xi = x[i];
	cross_param.yi = y[i];
	cross_param.xj = x[j];
	cross_param.yj = y[j];
	A[j + i * n] = A[i + j * n] = a * b * dblint(cross, &cross_param, 0, 1, 0, 2 * pi);
      }
    }
    for (i = 0; i < n; i++) {
      lin_param.x = x[i];
      lin_param.y = y[i];
      B[i] = a * b * dblint(lin, &lin_param, 0, 1, 0, 2 * pi);
    }
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++) {
	C[j + i * n] = 0;
	for (k = 0; k < n; k++)
	  C[j + i * n] += A[i + k * n] * A[j + k * n];
      }
    for (i = 0; i < n; i++)
      C[i + i * n] += gamma;
    vA = gsl_matrix_view_array(A, n, n);
    vB = gsl_vector_view_array(B, n);
    vC = gsl_matrix_view_array(C, n, n);
    vInvC = gsl_matrix_view_array(invC, n, n);
    vKsi = gsl_vector_view_array(Ksi, n);
    if ((p = gsl_permutation_alloc(n)) == NULL) {
      fprintf(stderr, "%s:%d: gsl_permutation_alloc failed\n", __FILE__, __LINE__);
      exit(1);
    }
    gsl_linalg_LU_decomp(&vC.matrix, p, &s);
    gsl_linalg_LU_invert(&vC.matrix, p, &vInvC.matrix);
    for (i = 0; i < n; i++) {
      Ksi[i] = 0;
      for (j = 0; j < n; j++)
	for (k = 0; k < n; k++)
	  Ksi[i] += invC[j + i * n] * A[k + j * n] * B[k];
    }
    for (i = 0; i < n; i++)
       printf ("%.16e %.16e %.16e\n", x[i], y[i], Ksi[i]);
    gsl_permutation_free(p);
    free(x);
    free(y);
    free(C);
    free(A);
    free(B);
    free(invC);
    free(Ksi);
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

  p = p0;
  x = a * r * cos(phi);
  y = b * r * sin(phi);
  xi = x - p->xi;
  yi = y - p->yi;
  xj = x - p->xj;
  yj = y - p->yj;
  ri = sqrt(xi*xi + yi*yi);
  rj = sqrt(xj*xj + yj*yj);

  return r * p->psi(ri) * p->psi(rj);
}

static double
lin(double r, double phi, void *p0)
{
  double x;
  double y;
  double xi;
  double yi;
  double ri;
  struct LinParam *p;

  p = p0;
  x = a * r * cos(phi);
  y = b * r * sin(phi);
  xi = x - p->x;
  yi = y - p->y;
  ri = sqrt(xi*xi + yi*yi);
  return r * p->vor(x, y) * p->psi(ri);
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
gauss(double r)
{
    double d2;
    d2 = delta * delta;
    return exp(-r * r / d2) / (pi * d2);
}

static double
chorin(double r)
{
  double d3;
  d3 = delta * delta * delta;
  return r < delta ? 3 * r / (2 * pi * d3) : 0;
}

static double
hald(double r)
{
    double ans;
    r /= delta;
    if (r < 0.001)
      ans = 15/8.0 - 21*r*r/32.0;
    else
      ans = 1/(r*r) * (4 * j2(2 * r) - j2(r));
    return ans / (3 * pi * delta * delta);
}

static double
j2(double x)
{
  return jn(2, x);
}
