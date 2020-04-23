#include "view.h"
#include "navier-stokes/centered.h"

#define MAXLEVEL (9)

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

static const char *me = "vertex";
static double contours[] = {0.25, 0.50, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

static double
f(double z, double q)
{
    double eps;
    eps = 1e-12;
    return fabs(z) < eps ? 0 : exp(-(q/z)*exp(1/(z - 1)));
}

static double
vorI(double x, double y)
{
    double a;
    double b;
    double Ksi;
    double q;
    double r;

    a = 0.8;
    b = 1.6;
    q = 2.56085;
    Ksi = 20;
    r = sqrt(sq(x/a) + sq(y/b));
    return r < 1 ? Ksi * (1 - f(r, q)) : 0;
}

static double
vorII(double x, double y)
{
    double a;
    double b;
    double Ksi;
    double r2;

    a = 0.8;
    b = 1.6;
    Ksi = 20;
    r2 = sq(x/a) + sq(y/b);
    return r2 < 1 ? 0 : Ksi * (1 - r2*r2);
}

int main()
{
    origin (-2, -2);
    size(4);
    init_grid (1 << MAXLEVEL);
    run();
}

event init (t = 0)
{
    scalar psi[], omega[];

    psi[left]   = dirichlet(0);
    psi[right]  = dirichlet(0);
    psi[top]    = dirichlet(0);
    psi[bottom] = dirichlet(0);

    foreach() {
        omega[] = vorII(x, y);
        psi[] = 0.;
    }
    boundary ({psi,omega});
    poisson (psi, omega);
    coord f = {-1.,1.};
    foreach()
        foreach_dimension()
        u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
    boundary ((scalar *){u});
}

event logfile (t <= 24; t += 1)
{
    scalar omega[];
    vorticity (u, omega);
    stats s = statsf (omega);
    fprintf (ferr, "%g %d %g %g %g %d\n", t, i, dt, s.sum, s.max, mgp.i);
}

event output (t += 0.5)
{
    enum {SIZE = 99};
    char name[SIZE];
    scalar omega[];
    FILE *f;
    int i;
    static int nf = 0;
    vorticity (u, omega);
    snprintf (name, SIZE, "omega-%03d", nf);
    if ((f = fopen (name, "w"))  == NULL) {
        fprintf(stderr, "%s: fail to open '%s'\n", me, name);
        exit(2);
    }
    output_field ({omega}, f, linear = true);
    fclose(f);

    view(width = 1200, height = 1200);
    clear();
    for (i = 0; i < sizeof contours/sizeof *contours; i++) 
        isoline ("omega", contours[i], lw = 2);
    snprintf (name, SIZE, "%03d.png", nf);
    save(name);
    nf++;
}

/*
#if TREE
event adapt (i++) {
    adapt_wavelet ((scalar *){u}, (double[]){5e-5,5e-5}, MAXLEVEL);
}
#endif
*/
