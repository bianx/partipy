#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vof.h>

static char me[] = "initial/grid";

static void
usg(void)
{
    fprintf(stderr,
            "%s > initial\n",
            me);
    exit(1);
}


int
main(int argc, char **argv)
{
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
}
