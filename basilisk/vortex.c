#include "navier-stokes/centered.h"

#define MAXLEVEL (8)

uf.n[left]   = 0.;
uf.n[right]  = 0.;
uf.n[top]    = 0.;
uf.n[bottom] = 0.;

static double
f(double z, double q) {
    double eps;
    eps = 1e-6;
    return fabs(z) < eps ? 0 : exp(-(q/z)*exp(1/(z - 1)));
}

static double
vorI(double x, double y) {
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
        omega[] = vorI(x, y);
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

event logfile (t <= 24; t += 1) {
    scalar omega[];
    vorticity (u, omega);
    stats s = statsf (omega);
    fprintf (ferr, "%g %d %g %g %g %d\n", t, i, dt, s.sum, s.max, mgp.i);
}

event output (t += 0.5) {
    static int nf = 0;
    scalar omega[];
    vorticity (u, omega);

    char name[80];
    sprintf (name, "omega-%03d", nf);
    FILE * fp = fopen (name, "w");
    output_field ({omega}, fp, linear = true);
    fclose (fp);

    nf++;
}

/*
#if TREE
event adapt (i++) {
    adapt_wavelet ((scalar *){u}, (double[]){5e-5,5e-5}, MAXLEVEL);
}
#endif */
