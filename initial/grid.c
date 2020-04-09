#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vof.h>

static char me[] = "initial/grid";

static void
usg(void)
{
    fprintf(stderr, "%s > initial\n", me);
    exit(1);
}


int
main(int argc, char **argv)
{
    double h;
    double ksi;
    double x;
    double y;
    int n;
    int i;
    int j;

    (void) argc;
    while (*++argv != NULL && argv[0][0] == '-')
        switch (argv[0][1]) {
        case 'h':
            usg();
            break;
        default:
            fprintf(stderr, "%s: unknown option '%s'\n", me, argv[0]);
            exit(2);
        }

    n = 5;
    h = 1.0 / n;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++) {
            ksi = vof_2d(i + 0.5, j + 0.5, n) * h * h;
            x = h * (i + 0.5);
            y = h * (j + 0.5);
            if (ksi > 10 * DBL_MIN) {
                printf("%.16e %.16e %.16e\n", x, y, ksi);
                printf("%.16e %.16e %.16e\n", x, -y, ksi);
                printf("%.16e %.16e %.16e\n", -x, y, ksi);
                printf("%.16e %.16e %.16e\n", -x, -y, ksi);
            }
        }
}
