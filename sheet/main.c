#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

enum {N = 400};
enum {M = 100};

static char me[] = "sheet/main";
static const double pi = 3.141592653589793;
static const double eps = 0.5;
static const double L = 1.0;
static const double t1 = 4.0;

int
func(double t, const double *z, double *f, void *params)
{
  (void)(params);
  (void)(t);
  const double *x;
  const double *y;
  double *fx;
  double *fy;
  double dx;
  double dy;
  double den;
  int i;
  int j;

  x = z;
  y = &z[N];
  fx = f;
  fy = &f[N];
  for (i = 0; i < N; i++)
    fx[i] = fy[i] = 0;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      if (i != j) {
	dx = x[i] - x[j];
	dy = y[i] - y[j];
	den = cosh(2*pi*dy) - cos(2*pi*dx) + eps*eps;
	fx[i] -= sinh(2*pi*dy)/den;
	fy[i] += sin(2*pi*dx)/den;
      }
  for (i = 0; i < N; i++) {
    fx[i] /= 2*N;
    fy[i] /= 2*N;
  }
  return GSL_SUCCESS;
}

int
main()
{
  int i;
  int j;
  double t;
  double z[2*N];
  double *x;
  double *y;
  double h;
  double ti;
  double dt;
  double x0;
  gsl_odeiv2_driver * driver;
  gsl_odeiv2_system sys = {func, NULL, 2*N, NULL};
  driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4,
				  1e-6, 1e-6, 0.0);
  x = z;
  y = &z[N];
  h = L/(N - 1);
  for (i = 0; i < N; i++) {
    x0 = i*h;
    x[i] = x0 + 0.01*sin(2*pi*x0);
    y[i] = -0.01*sin(2*pi*x0);
  }
  dt = t1 / (M - 1);
  t = 0;
  for (i = 0; i < M; i++) {
    fprintf(stderr, "%s: %d\n", me, i);
    ti = dt * i;
    if (gsl_odeiv2_driver_apply(driver, &t, ti, z) != GSL_SUCCESS) {
      fprintf(stderr, "%s: driver failed\n", me);
      exit(2);
      }
    if (i > 0)
      printf("\n");
    for (j = 0; j < N; j++)
      printf("%.16g %.16g\n", x[j], y[j]);
  }
  gsl_odeiv2_driver_free(driver);
  return 0;
}
