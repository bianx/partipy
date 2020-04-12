#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <barnes-hut.h>

static const char *me = "barnes_hut/example/info";
enum { SIZE = 999 };

enum {
    FORCE,
    LEAF,
    TWIG,
};

static const char *OutputName[] = {
    "force",
    "leaf",
    "twig",
};

static const int Output[] = {
    FORCE,
    LEAF,
    TWIG,
};

static void
usg(void)
{
    fprintf(stderr, "%s -p x y -t theta < points\n", me);
    exit(1);
}

int
main(int argc, char **argv)
{
    char line[SIZE];
    double *mass;
    double theta;
    double *x;
    double xh;
    double xl;
    double xp;
    double *y;
    double yh;
    double yl;
    double yp;
    int Pflag;
    int Tflag;
    int Oflag;
    long cap;
    long cnt;
    long i;
    long n;
    int output;
    struct BarnesHut *barnes_hut;
    struct BarnesHutInfo *info;
    const struct BarnesHutInfo *in;

    (void) argc;

    theta = 1;

    Pflag = Tflag = Oflag = 0;
    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        case 't':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -t\n", me);
                exit(2);
            }
            theta = atof(argv[0]);
            Tflag = 1;
            break;
        case 'p':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -p\n", me);
                exit(2);
            }
            xp = atof(argv[0]);
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -p\n", me);
                exit(2);
            }
            yp = atof(argv[0]);
            Pflag = 1;
            break;
        case 'o':
            argv++;
            if (argv[0] == NULL) {
                fprintf(stderr, "%s: not enough argumetns for -o\n", me);
                exit(2);
            }
            for (i = 0;; i++) {
                if (i == (int) (sizeof OutputName / sizeof *OutputName)) {
                    fprintf(stderr, "%s: unknown output type '%s'\n", me,
                            argv[0]);
                    exit(2);
                }
                if (strncmp(OutputName[i], argv[0], SIZE) == 0) {
                    output = Output[i];
                    break;
                }
            }
            Oflag = 1;
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }
    if (Pflag == 0) {
        fprintf(stderr, "%s: -p is not set\n", me);
        exit(2);
    }
    if (Tflag == 0) {
        fprintf(stderr, "%s: -t is not set\n", me);
        exit(2);
    }
    if (Oflag == 0) {
        fprintf(stderr, "%s: -o is not set\n", me);
        exit(2);
    }

    cap = 1;
    if ((x = malloc(cap * sizeof(*x))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((y = malloc(cap * sizeof *y)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((mass = malloc(cap * sizeof *mass)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }    
    n = 0;
    while (fgets(line, SIZE, stdin) != NULL) {
        if (n == cap) {
            cap *= 2;
            x = realloc(x, cap * sizeof *x);
            if (x == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(1);
            }
            y = realloc(y, cap * sizeof *y);
            if (y == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(1);
            }
            mass = realloc(mass, cap * sizeof *mass);
            if (mass == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                exit(1);
            }            
        }
        if (sscanf(line, "%lf %lf %lf\n", &x[n], &y[n], &mass[n]) != 3) {
            fprintf(stderr, "%s: fail to parse '%s'\n", me, line);
            exit(1);
        }
        n++;
    }
    if ((info = malloc(n * sizeof(*info))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        exit(1);
    }
    if ((barnes_hut = barnes_hut_build(n, x, y, mass)) == NULL) {
        fprintf(stderr, "%s:%d: barnes_hut_ini failed\n", __FILE__,
                __LINE__);
        exit(1);
    }
    barnes_hut_info(barnes_hut, theta, -1, xp, yp, &cnt, info);
    for (i = 0; i < cnt; i++) {
        in = &info[i];
        switch (output) {
        case FORCE:
            printf("%.16g %.16g\n", xp, yp);
            printf("%.16g %.16g\n", in->mx / in->m, in->my / in->m);
            printf("\n");
            break;
        case LEAF:
            if (!in->Coarse) {
                xl = in->x - in->w / 2;
                xh = in->x + in->w / 2;
                yl = in->y - in->w / 2;
                yh = in->y + in->w / 2;
                printf("%.16g %.16g\n", xl, yl);
                printf("%.16g %.16g\n", xh, yl);
                printf("%.16g %.16g\n", xh, yh);
                printf("%.16g %.16g\n", xl, yh);
                printf("%.16g %.16g\n", xl, yl);
                printf("\n");
            }
            break;
        case TWIG:
            if (in->Coarse) {
                xl = in->x - in->w / 2;
                xh = in->x + in->w / 2;
                yl = in->y - in->w / 2;
                yh = in->y + in->w / 2;
                printf("%.16g %.16g\n", xl, yl);
                printf("%.16g %.16g\n", xh, yl);
                printf("%.16g %.16g\n", xh, yh);
                printf("%.16g %.16g\n", xl, yh);
                printf("%.16g %.16g\n", xl, yl);
                printf("\n");
            }
            break;
        }
    }

    free(x);
    free(y);
    free(info);
    barnes_hut_fin(barnes_hut);
}
