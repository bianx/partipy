#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <barnes-hut.h>

//#include "float/real.h"
#include "double/real.h"

#define pi (3.141592653589793)
static char me[] = "gauss/main";
enum { nmax = 99999 };

static void
usg(void)
{
    fprintf(stderr,
            "%s -t time -m M -e every -d delta -c [chorin gauss hald j0 krasny] -s [euler rk4] -o [punto|skel|off|gnuplot] -a [n2|bh] -g [punto|vtk] -x nx -y ny [> punto] < initial\n",
            me);
    exit(2);
}

struct Core;
static int algorithm_n2(int n, const real *, real *, void *);
static int algorithm_bh(int n, const real *, real *, void *);
static int particle_n2(struct Core *, int n, const real * x,
                       const real * y, const real * ksi, real * ksi0);
static int particle_bh(struct Core *, int n, const real * x,
                       const real * y, const real * ksi, real * ksi0);
static real gauss_psi(real, real, void *);
static int gauss_dpsi(real, real, real *, real *, void *);
static real gauss_coef(void *);
static int chorin_dpsi(real, real, real *, real *, void *);
static real chorin_psi(real, real, void *);
static real chorin_coef(void *);
static int hald_dpsi(real, real, real *, real *, void *);
static real hald_coef(void *);
static real hald_psi(real, real, void *);
static int krasny_dpsi(real, real, real *, real *, void *);
static real krasny_coef(void *);
static int j0_dpsi(real, real, real *, real *, void *);
static real j0_coef(void *);

static int punto_grid(void *, int, const real *, const real *,
                      const real *, int step);
static int vtk_grid(void *, int, const real *, const real *,
                    const real *, int step);
static int null_grid(void *, int, const real *, const real *,
                     const real *, int step);

static real Theta;              /* Barnes-Hat parameter */

struct Core {
    real(*psi) (real, real, void *);
    int (*dpsi)(real, real, real *, real *, void *);
     real(*coef) (void *);
    void *param;
};

static struct Core Core[] = {
    { chorin_psi, chorin_dpsi, chorin_coef, NULL },
    { gauss_psi, gauss_dpsi, gauss_coef, NULL },
    { hald_psi, hald_dpsi, hald_coef, NULL },
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
                       const real *, int step);
static int vtk_write(int n, const real *, const real *, const real *,
                     const real *, int step);
static int skel_write(int n, const real *, const real *, const real *,
                      const real *, int step);
static int off_write(int n, const real *, const real *, const real *,
                     const real *, int step);
static int gnuplot_write(int n, const real *, const real *, const real *,
                         const real *, int step);

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

static int remesh_m4_n2(void *, int *, real *, real *, real *);
static int remesh_m4_bh(void *, int *, real *, real *, real *);
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
    "vtk",
};

static int (*const WriteFun[])(int, const real *, const real *,
                               const real *, const real *, int) = {
    gnuplot_write,
    off_write,
    punto_write,
    skel_write,
    vtk_write,
};

static const char *GridName[] = {
    "punto",
    "vtk",
    "null",
};

static int (*const Grid[])(void *, int, const real *, const real *,
                           const real *, int) = {
    punto_grid,
    vtk_grid,
    null_grid,
};

int
main(int argc, char **argv)
{
    char line[SIZE];
    const char *scheme;
    int Aflag;
    int Dflag;
    int Eflag;
    int every;
    int (*grid)(void *, int, const real *, const real *, const real *,
                int);
    int i;
    int j;
    int m;
    int Mflag;
    int n;
    int ncap;
    int nremesh;
    int nx;
    int ny;
    int (*remesh)(void *, int *, real *, real *, real *);
    int Sflag;
    int Tflag;
    int (*write)(int, const real *, const real *, const real *,
                 const real *, int);
    int (*algorithm)(int, const real *, real *, void *);
    int (*particle)(struct Core *, int, const real *, const real *,
                    const real *, real *);
    real *buf;
    real delta;
    real dt;
    real *ksi;
    real *ksi0;
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

    (void) argc;

    Aflag = Eflag = Dflag = Mflag = Tflag = Sflag = 0;
    write = NULL;
    core = NULL;
    scheme = NULL;
    nremesh = 0;
    grid = null_grid;
    nx = ny = 0;
    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        case 'a':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -m needs an argument\n", me);
                exit(2);
            }
            if (strncmp(argv[0], "n2", SIZE) == 0) {
                algorithm = algorithm_n2;
                particle = particle_n2;
                remesh = remesh_m4_n2;
            } else if (strncmp(argv[0], "bh", SIZE) == 0) {
                argv++;
                if (argv[0] == NULL) {
                    fprintf(stderr,
                            "%s: barnes-hut algorithm needs an argument 'theta'\n",
                            me);
                    exit(2);
                }
                Theta = atof(argv[0]);
                algorithm = algorithm_bh;
                particle = particle_bh;
                remesh = remesh_m4_bh;
            } else {
                fprintf(stderr, "%s: unknown algorithm '%s'\n", me,
                        argv[0]);
                exit(2);
            }
            Aflag = 1;
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
        case 'g':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -g needs an argument\n", me);
                exit(2);
            }
            for (i = 0;; i++) {
                if (i == sizeof(GridName) / sizeof(*GridName)) {
                    fprintf(stderr, "%s: unknown grid name '%s'\n", me,
                            argv[0]);
                    exit(2);
                }
                if (strncmp(argv[0], GridName[i], SIZE) == 0) {
                    grid = Grid[i];
                    break;
                }
            }
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
        case 'x':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -x needs an argument\n", me);
                exit(2);
            }
            nx = atoi(argv[0]);
            break;
        case 'y':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: -y needs an argument\n", me);
                exit(2);
            }
            ny = atoi(argv[0]);
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }
    if (Aflag == 0) {
        fprintf(stderr, "%s: -a is not given\n", me);
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

    if (nremesh != 0 || grid != null_grid) {
        if (nx == 0 || ny == 0) {
            fprintf(stderr, "%s: unset grid sizes x = %d, y = %d\n", me,
                    nx, ny);
            exit(2);
        }
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
    if ((ksi0 = malloc(nmax * sizeof(*ksi))) == NULL) {
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
    ode_param.function = algorithm;
    ode_param.param = &param;
    ode_param.dt = dt;
    ode_param.scheme = scheme;

    remesh_param.nx = nx;
    remesh_param.ny = ny;
    remesh_param.xlo = -1.6;
    remesh_param.xhi = 1.6;
    remesh_param.ylo = -1.6;
    remesh_param.yhi = 1.6;
    remesh_param.eps = 1e-3;
    remesh_param.core = core;

    if (ode_ini(argv, &ode_param, &ode) != 0) {
        fprintf(stderr, "%s: ode_ini failed\n", me);
        exit(2);
    }
    for (i = 0;; i++) {
        if (i > 0 && nremesh > 0 && i % nremesh == 0)
            remesh(&remesh_param, &n, x, y, ksi);
        if (every > 0 && i % every == 0) {
            particle(core, n, x, y, ksi, ksi0);
            if (write(n, x, y, ksi, ksi0, i) != 0) {
                fprintf(stderr, "%s: write failed\n", me);
                exit(2);
            }
            grid(&remesh_param, n, x, y, ksi, i);
        }
        if (i == m)
            break;

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
    free(ksi0);
    ode_fin(ode);
}

static int
algorithm_n2(int n, const real * z, real * f, void *params0)
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
algorithm_bh(int n, const real * z, real * f, void *params0)
{
    const real *ksi;
    const real *x;
    const real *y;
    real *ksi0;
    real *x0;
    real *y0;
    int i;
    int j;
    long cnt;
    real coef;
    real dx;
    real dy;
    real *fx;
    real *fy;
    real gx;
    real gy;
    struct BarnesHut *barnes_hut;
    struct Core *core;
    struct Param *params;

    params = params0;
    core = params->core;
    ksi = params->ksi;
    coef = core->coef(core->param);
    x = z;
    y = &z[n];
    fx = f;
    fy = &f[n];


    if ((x0 = malloc(n * sizeof *x0)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((y0 = malloc(n * sizeof *y0)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((ksi0 = malloc(n * sizeof *ksi0)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((barnes_hut = barnes_hut_build(n, x, y, ksi)) == NULL) {
        fprintf(stderr, "%s: barnes_hut_build failed\n", me);
        return 1;
    }
    for (i = 0; i < n; i++)
        fx[i] = fy[i] = 0;
    for (i = 0; i < n; i++) {
        if (barnes_hut_interaction
            (barnes_hut, Theta, i, x[i], y[i], &cnt, x0, y0, ksi0) != 0) {
            fprintf(stderr, "%s: barnes_hut_interaction failed\n", me);
            return 1;
        }
        //fprintf(stderr, "cnt: %ld %ld\n", cnt, n);
        for (j = 0; j < cnt; j++) {
            dx = x[i] - x0[j];
            dy = y[i] - y0[j];
            core->dpsi(dx, dy, &gx, &gy, core->param);
            fx[i] -= ksi0[j] * gy;
            fy[i] += ksi0[j] * gx;
        }
    }
    for (i = 0; i < n; i++) {
        fx[i] *= coef;
        fy[i] *= coef;
    }

    free(x0);
    free(y0);
    free(ksi0);
    barnes_hut_fin(barnes_hut);
    return 0;
}

static int
vtk_write(int n, const real * x, const real * y, const real * ksi,
          const real * ksi0, int step)
{
    char name[SIZE];
    FILE *f;
    int i;
    int status;

    if (snprintf(name, SIZE, "p.%06d.vtk", step) < 0)
        return 1;
    if ((f = fopen(name, "w")) == NULL) {
        fprintf(stderr, "%s: fail to open file '%s'\n", me, name);
        return 1;
    }
    fprintf(f,
            "# vtk DataFile Version 2.0\n"
            "generated by %s\n"
            "ASCII\n"
            "DATASET UNSTRUCTURED_GRID\n" "POINTS %d float\n", me, n);
    for (i = 0; i < n; i++) {
        status = fprintf(f, "%.16e %.16e %.16e\n", x[i], y[i], 0.0);
        if (status < 0) {
            fprintf(stderr, "%s:%d: fail to write '%s'\n", __FILE__,
                    __LINE__, name);
            return 1;
        }
    }
    fprintf(f, "CELL_TYPES %d\n", n);
    for (i = 0; i < n; i++)
        fprintf(f, "1\n");
    fprintf(f, "POINT_DATA %d\n", n);
    fprintf(f, "SCALARS ksi float 1\n" "LOOKUP_TABLE default\n");
    for (i = 0; i < n; i++)
        fprintf(f, "%.16g\n", ksi[i]);
    fprintf(f, "SCALARS ksi0 float 1\n" "LOOKUP_TABLE default\n");
    for (i = 0; i < n; i++)
        fprintf(f, "%.16g\n", ksi0[i]);

    if (fclose(f) != 0) {
        fprintf(stderr, "%s:%d: fail to close '%s'\n", __FILE__, __LINE__,
                name);
        return 1;
    }
    return 0;
}

static int
punto_write(int n, const real * x, const real * y, const real * ksi,
            const real * ksi0, int step)
{
    int j;

    if (step > 0)
        printf("\n");
    for (j = 0; j < n; j++)
        printf("%.16e %.16e %.16e %.16e\n", x[j], y[j], ksi[j], ksi0[j]);
    return 0;
}

static int
skel_write(int n, const real * x, const real * y, const real * ksi,
           const real * ksi0, int step)
{
    char path[SIZE];
    real z;
    FILE *f;
    int i;
    int npolylines;

    (void) ksi0;
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
          const real * ksi0, int step)
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
    (void) ksi0;

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
              const real * ksi0, int step)
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
        fprintf(f, "%.16g %.16g %.16g %.16g\n", x[i], y[i], ksi[i],
                ksi0[i]);
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

static real
j2(real x)
{
    return jn(2, x);
}

static real
hald_psi(real x, real y, void *p0)
{
    real d;
    real r;
    real ans;
    struct PsiParam *p;

    p = p0;
    d = p->delta;
    r = sqrtr(x * x + y * y) / d;
    if (r < 0.001)
        ans = 15 / 8.0 - 21 * r * r / 32.0;
    else
        ans = 1 / (r * r) * (4 * j2(2 * r) - j2(r));
    return ans / (3 * pi * d * d);
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
    return r < delta ? 3 * r / (2 * pi * d3) : 0;
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
    if (q->function(n / 2, y, k, q->param) != 0) {
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
lam1(real cutoff, real x)
{
    real a;

    x = fabs(x) / cutoff;
    a = 1 - x;

    if (x < 1)
        return a / cutoff;
    else
        return 0;
}

static real
lam2(real cutoff, real x)
{
    real a;
    real b;

    x = fabs(x) / cutoff;
    a = 1 - x;
    b = 2 - x;

    if (x < 0.5)
        return (1 - x * x) / cutoff;
    else if (x < 1.5)
        return a * b / 2 / cutoff;
    else
        return 0;
}

static real
lam3(real cutoff, real x)
{
    real a;
    real b;
    real c;

    x = fabs(x) / cutoff;
    a = 1 - x;
    b = 2 - x;
    c = 3 - x;

    if (x < 1)
        return (1 - x * x) * b / 2 / cutoff;
    else if (x < 2)
        return a * b * c / 6 / cutoff;
    else
        return 0;
}

static real
m4(real cutoff, real x)
{
    real a;
    real b;

    x = fabs(x) / cutoff;
    a = 2 - x;
    b = 1 - x;
    if (x < 1)
        return (a * a * a / 6 - 4 * b * b * b / 6) / cutoff;
    else if (x < 2)
        return a * a * a / 6 / cutoff;
    else
        return 0;
}

static real
m4p(real cutoff, real x)
{
    real a;

    //This is Eq.(19) (M prime kernel) in Koumoutsakos J. Comput. Phys. 1996
    x = fabs(x) / cutoff;
    a = 2 - x;

    if (x < 1)
        return (1 - 5 * x * x / 2 + 3 * x * x * 3 / 2) / cutoff;
    else if (x < 2)
        return (2 - a) * (2 - a) * (1 - a) / 2 / cutoff;
    else
        return 0;
}

static int
remesh_m4_n2(void *p0, int *pn, real * x, real * y, real * ksi)
{
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    int nx;
    int ny;
    real coef;
    real dx;
    real dy;
    real eps;
    real *ksi0;
    real u;
    real v;
    real xhi;
    real xlo;
    real yhi;
    real ylo;
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
    if ((ksi0 = malloc(m * sizeof(*ksi0))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }
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

    real ksi_m;

    ksi_m = ksi[0];
    for (i = 0; i < m; i++)
        if (ksi[i] > ksi_m)
            ksi_m = ksi[i];

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
        if (fabs(ksi[i]) > ksi_m * eps) {
            x[j] = x[i];
            y[j] = y[i];
            ksi[j] = ksi[i];
            j++;
        }
    *pn = j;
    free(ksi0);
    return 0;
}

static int
remesh_m4_bh(void *p0, int *pn, real * x, real * y, real * ksi)
{
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;
    int nx;
    int ny;
    long cnt;
    real coef;
    real dx;
    real dy;
    real eps;
    real *ksi0;
    real *ksi00;
    real ksi_m;
    real u;
    real v;
    real *x0;
    real xhi;
    real xlo;
    real *y0;
    real yhi;
    real ylo;
    struct BarnesHut *barnes_hut;
    struct RemeshParam *p;

    n = *pn;
    if ((x0 = malloc(n * sizeof *x0)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((y0 = malloc(n * sizeof *y0)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((ksi00 = malloc(n * sizeof *ksi00)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((barnes_hut = barnes_hut_build(n, x, y, ksi)) == NULL) {
        fprintf(stderr, "%s: barnes_hut_build failed\n", me);
        return 1;
    }

    
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
    if ((ksi0 = malloc(m * sizeof(*ksi0))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }

    l = 0;
    for (i = 0; i < nx; i++) {
        u = xlo + (i + 0.5) * dx;
        for (j = 0; j < ny; j++) {
            v = ylo + (j + 0.5) * dy;
            ksi0[l] = 0;
            if (barnes_hut_interaction
                (barnes_hut, Theta, -1, u, v, &cnt, x0, y0, ksi00) != 0) {
                fprintf(stderr, "%s: barnes_hut_interaction failed\n", me);
                return 1;
            }
            for (k = 0; k < cnt; k++)
                ksi0[l] += ksi00[k] * m4(dx, x0[k] - u) * m4(dy, y0[k] - v);
            l++;
        }
    }
    coef = dx * dy;
    for (i = 0; i < m; i++)
        ksi[i] = ksi0[i] * coef;

    ksi_m = ksi[0];
    for (i = 0; i < m; i++)
        if (ksi[i] > ksi_m)
            ksi_m = ksi[i];
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
        if (fabs(ksi[i]) > ksi_m * eps) {
            x[j] = x[i];
            y[j] = y[i];
            ksi[j] = ksi[i];
            j++;
        }
    *pn = j;
    free(x0);
    free(y0);
    free(ksi0);
    free(ksi00);
    return 0;
}


static int
remesh_psi(void *p0, int *pn, real * x, real * y, real * ksi)
{
    int i;
    int j;
    int k;
    int l;
    int n;
    int m;
    int nx;
    int ny;
    real coef;
    real *ksi0;
    real dx;
    real dy;
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
    if ((ksi0 = malloc(m * sizeof(*ksi0))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }
    for (i = 0; i < m; i++)
        ksi0[i] = 0;

    for (k = 0; k < n; k++) {
        l = 0;
        for (i = 0; i < nx; i++) {
            u = xlo + (i + 0.5) * dx;
            for (j = 0; j < ny; j++) {
                v = ylo + (j + 0.5) * dy;
                ksi0[l] +=
                    ksi[k] * core->psi(x[k] - u, y[k] - v, core->param);
                l++;
            }
        }
    }
    coef = dx * dy;
    for (i = 0; i < m; i++)
        ksi[i] = ksi0[i] * coef;

    real ksi_m;

    ksi_m = ksi[0];
    for (i = 0; i < m; i++)
        if (ksi[i] > ksi_m)
            ksi_m = ksi[i];

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
        if (fabs(ksi[i]) > ksi_m * eps) {
            x[j] = x[i];
            y[j] = y[i];
            ksi[j] = ksi[i];
            j++;
        }
    free(ksi0);
    *pn = j;
    return 0;
}

static int
particle_n2(struct Core *core, int n, const real * x, const real * y,
            const real * ksi, real * ksi0)
{
    int i;
    int j;

    if (core->psi == NULL)
        return 0;

    for (i = 0; i < n; i++) {
        ksi0[i] = 0;
        for (j = 0; j < n; j++)
            ksi0[i] +=
                ksi[j] * core->psi(x[i] - x[j], y[i] - y[j], core->param);
    }
    return 0;
}

static int
particle_bh(struct Core *core, int n, const real * x, const real * y,
            const real * ksi, real * ksi0)
{
    int i;
    int j;
    real *x0;
    real *y0;
    real *ksi00;
    real dx;
    real dy;
    real f;
    struct BarnesHut *barnes_hut;
    long cnt;

    if (core->psi == NULL)
        return 0;

        if ((x0 = malloc(n * sizeof *x0)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((y0 = malloc(n * sizeof *y0)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((ksi00 = malloc(n * sizeof *ksi00)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if ((barnes_hut = barnes_hut_build(n, x, y, ksi)) == NULL) {
        fprintf(stderr, "%s: barnes_hut_build failed\n", me);
        return 1;
    }
    for (i = 0; i < n; i++) {
        ksi0[i] = 0;
        if (barnes_hut_interaction
            (barnes_hut, Theta, i, x[i], y[i], &cnt, x0, y0, ksi00) != 0) {
            fprintf(stderr, "%s: barnes_hut_interaction failed\n", me);
            return 1;
        }
        for (j = 0; j < cnt; j++) {
            dx = x[i] - x0[j];
            dy = y[i] - y0[j];
            f = core->psi(dx, dy, core->param);
            ksi0[i] += ksi00[j] * f;
        }
    }

    free(x0);
    free(y0);
    free(ksi00);
    barnes_hut_fin(barnes_hut);
    return 0;
}

static int
punto_grid(void *p0, int n, const real * x, const real * y,
           const real * ksi, int step)
{
    FILE *f;
    int i;
    int j;
    int k;
    int l;
    int m;
    int nx;
    int ny;
    real dx;
    real dy;
    real *ksi0;
    real u;
    real v;
    real xhi;
    real xlo;
    real yhi;
    real ylo;
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
    if ((ksi0 = malloc(m * sizeof(*ksi0))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }
    for (i = 0; i < m; i++)
        ksi0[i] = 0;
    for (k = 0; k < n; k++) {
        l = 0;
        for (j = 0; j < ny; j++) {
            v = ylo + (j + 0.5) * dy;
            for (i = 0; i < nx; i++) {
                u = xlo + (i + 0.5) * dx;
                ksi0[l] +=
                    ksi[k] * core->psi(x[k] - u, y[k] - v, core->param);
                l++;
            }
        }
    }
    snprintf(path, SIZE, "%06d.grid", step);
    fprintf(stderr, "%s: write '%s'\n", me, path);
    if ((f = fopen(path, "w")) == NULL) {
        fprintf(stderr, "%s: fail to open '%s'\n", me, path);
        exit(2);
    }
    l = 0;
    for (j = 0; j < ny; j++) {
        v = ylo + (j + 0.5) * dy;
        for (i = 0; i < nx; i++) {
            u = xlo + (i + 0.5) * dx;
            fprintf(f, "%.16e %.16e %.16e\n", u, v, ksi0[l++]);
        }
    }

    free(ksi0);
    fclose(f);
    return 0;
}

static int
null_grid(void *p0, int n, const real * x, const real * y,
          const real * ksi, int step)
{
    (void) p0;
    (void) n;
    (void) x;
    (void) y;
    (void) ksi;
    (void) step;
    return 0;
}

static int
vtk_grid(void *p0, int n, const real * x, const real * y, const real * ksi,
         int step)
{
    char path[SIZE];
    FILE *f;
    int i;
    int j;
    int k;
    int l;
    int m;
    int nx;
    int ny;
    int status;
    real dx;
    real dy;
    real *ksi0;
    real u;
    real v;
    real xhi;
    real xlo;
    real yhi;
    real ylo;
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
    if ((ksi0 = malloc(m * sizeof(*ksi0))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(2);
    }
    for (i = 0; i < m; i++)
        ksi0[i] = 0;
    for (k = 0; k < n; k++) {
        l = 0;
        for (j = 0; j < ny; j++) {
            v = ylo + (j + 0.5) * dy;
            for (i = 0; i < nx; i++) {
                u = xlo + (i + 0.5) * dx;
                ksi0[l] +=
                    ksi[k] * core->psi(x[k] - u, y[k] - v, core->param);
                l++;
            }
        }
    }
    snprintf(path, SIZE, "g.%06d.vtk", step);
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
                     me, nx, ny, xlo + dx / 2, ylo + dy / 2, dx, dy);
    if (status < 0) {
        fprintf(stderr, "%s: fail to write '%s'\n", me, path);
        exit(2);
    }
    fprintf(f, "POINT_DATA %d\n"
            "SCALARS omega double\n" "LOOKUP_TABLE DEFAULT\n", m);
    for (i = 0; i < m; i++)
        fprintf(f, "%.16e\n", ksi0[i]);
    free(ksi0);
    fclose(f);
    return 0;
}
