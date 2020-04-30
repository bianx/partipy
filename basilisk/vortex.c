#include "view.h"
#include "navier-stokes/centered.h"

#define MAXLEVEL (9)
static const char *me = "vertex";
static double contours[] = {0.25, 0.50, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
static double (*vor)(double, double);
static double vorI(double, double);
static double vorII(double, double);
static double f(double, double);

static void
usg(void)
{
    fprintf(stderr, "%s [-1|-2]\n", me);
    exit(2);
}

uf.n[left]   = 0;
uf.n[right]  = 0;
uf.n[top]    = 0;
uf.n[bottom] = 0;

u.n[left] = neumann(0);
u.n[right] = neumann(0);
u.n[top] = neumann(0);
u.n[bottom] = neumann(0);

int main(int argc, char **argv)
{
    int Vflag;
    (void)argc;
    Vflag = 0;
    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        case '1':
            vor = vorI;
            Vflag = 1;
            break;
        case '2':
            vor = vorII;
            Vflag = 1;
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }
    if (!Vflag) {
        fprintf(stderr, "%s: -1 or -2 must be given\n", me);
        exit(2);
    }
    origin(-4, -4);
    size(8);
    init_grid(1 << MAXLEVEL);
    run();
}

event init (t = 0)
{
    scalar psi[];
    scalar omega[];

    psi[left]   = dirichlet(0);
    psi[right]  = dirichlet(0);
    psi[top]    = dirichlet(0);
    psi[bottom] = dirichlet(0);
    foreach() {
        omega[] = vor(x, y);
        psi[] = 0;
    }
    boundary({psi,omega});
    poisson(psi, omega);
    coord f = {-1.,1.};
    foreach()
        foreach_dimension()
        u.x[] = f.x*(psi[0,1] - psi[0,-1])/(2.*Delta);
    boundary((scalar *){u});
}

event logfile (t <= 24; t += 1)
{
    scalar omega[];
    stats s;
    vorticity(u, omega);
    s = statsf(omega);
    fprintf(stderr, "%g %d %g %g %g %d\n", t, i, dt, s.sum, s.max, mgp.i);
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

    view(width = 1200, height = 1200, fov = 10);
    clear();
    for (i = 0; i < sizeof contours/sizeof *contours; i++)
        isoline ("omega", contours[i], lw = 2);
    snprintf (name, SIZE, "%03d.png", nf);
    save(name);
    nf++;
}

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
    return r2 <= 1 ? Ksi * (1 - r2 * r2) : 0;
}

/*
#if TREE
event adapt (i++) {
    adapt_wavelet ((scalar *){u}, (double[]){5e-5,5e-5}, MAXLEVEL);
}
#endif
*/
