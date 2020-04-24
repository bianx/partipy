#include "view.h"
#include "navier-stokes/centered.h"

#define MAXLEVEL (7)

u.t[bottom] = dirichlet(0);

static const char *me = "vertex";
static double contours[] = {0.25, 0.50, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
#define RAD (pow(pow((x), 2)+(pow((y), 2)), 0.5))
#define ST (-(x)/RAD)

int main()
{
    origin (-10, -2.5);
    size(20);
    init_grid (1 << MAXLEVEL);
    mask (y > 2.5 ? top : none);
    run();
}

event init (t = 0)
{
    scalar psi[];
    coord f = {-1.,1.};
    double k;

    k = 3.8317;
    foreach() {
        psi[] = ((RAD > 1)*ST/RAD +
                 (RAD < 1)*(-2*j1(k*RAD)*ST/(k*j0(k)) + RAD*ST));
    }
    boundary ({psi});
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
