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

struct BarnesHut {
    struct Node **node;
    long cnt;
    long cap;
};

struct BarnesHut *
barnes_hut_ini(long n, const double *x, const double *y)
{
  struct BarnesHut *q;

  if ((q = malloc(sizeof(*q))) == NULL) {
    fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
    return NULL;
  }

  q->cap = 1;
  q->cnt = 0;
  if ((q->node = malloc(q->cap * sizeof *q->node)) == NULL) {
	fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
	return NULL;
  }



  return q;
}
