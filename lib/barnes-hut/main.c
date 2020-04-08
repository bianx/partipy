#include <stdlib.h>
#include <stdio.h>
#include "barnes-hut.h"

enum { NE, NW, SW, SE };
struct Node {
    double x;
    double y;
    long data;
    struct Node *elm[4];
};

struct Quadtree {
    long *out;
    long out_cnt;
    long out_capacity;
    struct Node **node;
    long node_cnt;
    long node_capacity;
};
struct Region {
    double L;
    double R;
    double B;
    double T;
};
static int compare(struct Node *, struct Node *);
static int in_region(double, double, const struct Region *rectangle);
static int rectangle_overlaps_region(double, double, double, double,
                                     const struct Region *);
static int insert(struct Node *, struct Node *root);
static int search(struct Quadtree *q, struct Node *P, double L, double R,
                  double B, double T, const struct Region *);
static struct Node *node_ini(double x, double y, long data);

struct Quadtree *
quadtree_ini(void)
{
    struct Quadtree *q;
    long *out;
    long out_capacity;
    long node_capacity;
    struct Node **node;

    if ((q = malloc(sizeof(*q))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return NULL;
    }
    out_capacity = 1;
    node_capacity = 1;
    if ((out = malloc(out_capacity * sizeof(*out))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return NULL;
    }
    if ((node = malloc(out_capacity * sizeof(*node))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return NULL;
    }
    q->node = node;
    q->node_capacity = node_capacity;
    q->node_cnt = 0;
    q->out = out;
    q->out_capacity = out_capacity;
    return q;
}

int
quadtree_fin(struct Quadtree *q)
{
    long i;

    for (i = 0; i < q->node_cnt; i++)
        free(q->node[i]);
    free(q->out);
    free(q->node);
    free(q);
    return 0;
}

int
quadtree_build(struct Quadtree *q, long n, const double *x,
               const double *y)
{
    long i;

    for (i = 0; i < n; i++)
        if (quadtree_insert(q, x[i], y[i], i) != 0)
            return 1;
    return 0;
}

int
quadtree_insert(struct Quadtree *q, double x, double y, long d)
{
    struct Node *no;
    long node_cnt;
    long node_capacity;
    struct Node **node;

    node_cnt = q->node_cnt;
    node_capacity = q->node_capacity;
    node = q->node;
    if (node_cnt == node_capacity) {
        node_capacity *= 2;
        node = realloc(node, node_capacity * sizeof(*node));
        if (node == NULL) {
            fprintf(stderr, "%s:%d: realloc failed\n", __FILE__, __LINE__);
            return 1;
        }
    }
    if ((no = node_ini(x, y, d)) == NULL)
        return 1;
    node[node_cnt] = no;
    if (node_cnt > 0)
        insert(node[node_cnt], node[0]);
    node_cnt++;
    q->node_cnt = node_cnt;
    q->node_capacity = node_capacity;
    q->node = node;
    return 0;
}

int
quadtree_rectangle(struct Quadtree *q, double xlo, double xhi, double ylo,
                   double yhi, long *pn, const long **pout)
{
    struct Region region;

    region.L = xlo;
    region.R = xhi;
    region.B = ylo;
    region.T = yhi;

    q->out_cnt = 0;
    if (search(q, q->node[0], xlo, xhi, ylo, yhi, &region) != 0) {
        fprintf(stderr, "%s:%d: search failed\n", __FILE__, __LINE__);
        return 1;
    }

    *pn = q->out_cnt;
    *pout = q->out;
    return 0;
}


static int
insert(struct Node *K, struct Node *R)
{
    int direction;

    direction = compare(R, K);
    while (R->elm[direction] != NULL) {
        R = R->elm[direction];
        direction = compare(R, K);
    }
    R->elm[direction] = K;
    return 0;
}


static int
search(struct Quadtree *q, struct Node *P, double L, double R, double B,
       double T, const struct Region *region)
{
    double x;
    double y;
    struct Node **elm;

    elm = P->elm;
    x = P->x;
    y = P->y;
    if (in_region(x, y, region)) {  /* found */
        if (q->out_cnt == q->out_capacity) {
            q->out_capacity *= 2;
            q->out = realloc(q->out, q->out_capacity * sizeof(*q->out));
            if (q->out == NULL) {
                fprintf(stderr, "%s:%d: realloc failed\n", __FILE__,
                        __LINE__);
                return 1;
            }
        }
        q->out[q->out_cnt] = P->data;
        q->out_cnt++;
    }
    if (elm[NE] != NULL && rectangle_overlaps_region(x, R, y, T, region))
        if (search(q, elm[NE], x, R, y, T, region) != 0)
            return 1;
    if (elm[NW] != NULL && rectangle_overlaps_region(L, x, y, T, region))
        if (search(q, elm[NW], L, x, y, T, region) != 0)
            return 1;
    if (elm[SW] != NULL && rectangle_overlaps_region(L, x, B, y, region))
        if (search(q, elm[SW], L, x, B, y, region) != 0)
            return 1;
    if (elm[SE] != NULL && rectangle_overlaps_region(x, R, B, y, region))
        if (search(q, elm[SE], x, R, B, y, region) != 0)
            return 1;
    return 0;
}

static struct Node *
node_ini(double x, double y, long data)
{
    struct Node *node;

    if ((node = malloc(sizeof(*node))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return NULL;
    }
    node->x = x;
    node->y = y;
    node->data = data;
    node->elm[NE] = node->elm[NW] = node->elm[SW] = node->elm[SE] = NULL;
    return node;
}

static int
compare(struct Node *a, struct Node *b)
{
    if (b->x > a->x)
        return b->y > a->y ? NE : SE;
    else
        return b->y > a->y ? NW : SW;
}

static int
in_region(double x, double y, const struct Region *rectangle)
{
    double BP;
    double LP;
    double RP;
    double TP;

    BP = rectangle->B;
    LP = rectangle->L;
    RP = rectangle->R;
    TP = rectangle->T;
    return LP <= x && x <= RP && BP <= y && y <= TP;
}

static int
rectangle_overlaps_region(double L, double R, double B, double T,
                          const struct Region *region)
{
    double BP;
    double LP;
    double RP;
    double TP;

    BP = region->B;
    LP = region->L;
    RP = region->R;
    TP = region->T;
    return L <= RP && R >= LP && B <= TP && T >= BP;
}
