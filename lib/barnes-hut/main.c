#include <stdlib.h>
#include <stdio.h>
#include "barnes-hut.h"

enum { NE, NW, SW, SE };
struct Node {
    double x;
    double y;
    double w;
    double m;
    double mx;
    double my;
    long id;
    long cnt;
    struct Node *elm[4];
};
