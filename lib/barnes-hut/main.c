#include <stdlib.h>
#include <stdio.h>
#include "barnes-hut.h"

enum { NE, NW, SW, SE };
struct Node {
    double m;
    double x;
    double y;
    long id;
    double xl;
    double xh;
    double yl;
    double yh;
    struct Node *elm[4];
};

struct BarnesHut {
    double xh;
    double xl;
    double yh;
    double yl;
    long cap;
    long cnt;
    struct Node **node;
};

static int quadrant(struct Node *, double, double);
static int box(struct Node *, int, double *, double *, double *, double *);
static struct Node *node_ini(struct BarnesHut *);
static int insert(struct BarnesHut *, struct Node *, double, double,
                  double, long);
static int walk(struct BarnesHut *, struct Node *, int (*)(struct BarnesHut *, struct Node *, void *), void *);

struct BarnesHut *
barnes_hut_ini(double xl, double xh, double yl, double yh)
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
    q->xl = xl;
    q->xh = xh;
    q->yl = yl;
    q->yh = yh;
    return q;
}

int
barnes_hut_fin(struct BarnesHut *q)
{
    long i;

    for (i = 0; i < q->cnt; i++)
        free(q->node[i]);
    free(q->node);
    free(q);
    return 0;
}

int
barnes_hut_insert(struct BarnesHut *q, double x, double y, double m,
                  long id)
{
    struct Node *no;

    if (q->cnt == 0) {
        if ((no = node_ini(q)) == NULL) {
            fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
            return 1;
        }
        no->m = m;
        no->x = x;
        no->y = y;
        no->id = id;
        no->xl = q->xl;
        no->xh = q->xh;
        no->yl = q->yl;
        no->yh = q->yh;
        no->elm[NE] = NULL;
        no->elm[NW] = NULL;
        no->elm[SW] = NULL;
        no->elm[SE] = NULL;
    } else {
        insert(q, q->node[0], x, y, m, id);
    }
    return 0;
}

static int
print(struct BarnesHut *q, struct Node *n, void *f0)
{
    FILE *f;
    (void)q;

    f = f0;

    if (n->elm[NE] == NULL && n->elm[NW] == NULL
        && n->elm[SW] == NULL && n->elm[SE] == NULL) {
        if (fprintf(f, "%.16g %.16g\n", n->xl, n->yl) < 0)
            return 1;
        if (fprintf(f, "%.16g %.16g\n", n->xh, n->yl) < 0)
            return 1;
        if (fprintf(f, "%.16g %.16g\n", n->xh, n->yh) < 0)
            return 1;
        if (fprintf(f, "%.16g %.16g\n", n->xl, n->yh) < 0)
            return 1;
        if (fprintf(f, "%.16g %.16g\n", n->xl, n->yl) < 0)
            return 1;
        if (fprintf(f, "\n") < 0)
            return 1;
    }
    return 0;
}

int
barnes_hut_print(struct BarnesHut *q, FILE *f)
{
    if (q->cnt > 0)
        return walk(q, q->node[0], print, f);
    else
        return 0;
}

static int
quadrant(struct Node *n, double x, double y)
{
    double u;
    double v;

    u = (n->xl + n->xh) / 2;
    v = (n->yl + n->yh) / 2;
    if (x > u)
        return y > v ? NE : SE;
    else
        return y > v ? NW : SW;
}

static int
box(struct Node *n, int i, double *xl, double *xh, double *yl, double *yh)
{
    double x;
    double y;

    x = (n->xl + n->xh) / 2;
    y = (n->yl + n->yh) / 2;
    switch (i) {
    case NE:
        *xl = x;
        *xh = n->xh;
        *yl = y;
        *yh = n->yh;
        break;
    case SE:
        *xl = x;
        *xh = n->xh;
        *yl = n->yl;
        *yh = y;
        break;
    case NW:
        *xl = n->xl;
        *xh = x;
        *yl = y;
        *yh = n->yh;
        break;
    case SW:
        *xl = n->xl;
        *xh = x;
        *yl = n->yl;
        *yh = y;
        break;
    }
    return 0;
}

static struct Node *
node_ini(struct BarnesHut *q)
{
    struct Node *no;

    if ((no = malloc(sizeof *no)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return NULL;
    }
    if (q->cnt == q->cap) {
        q->cap *= 2;
        if ((q->node = realloc(q->node, q->cap * sizeof *q->node)) == NULL) {
            fprintf(stderr, "%s:%d: realloc failed\n", __FILE__, __LINE__);
            return NULL;
        }
    }
    no->elm[NE] = NULL;
    no->elm[NW] = NULL;
    no->elm[SW] = NULL;
    no->elm[SE] = NULL;
    q->node[q->cnt++] = no;
    return no;
}

static int
insert(struct BarnesHut *q, struct Node *root, double x, double y,
       double m, long id)
{
    int i;
    struct Node *no;

    if (root->elm[NE] == NULL && root->elm[NW] == NULL
        && root->elm[SW] == NULL && root->elm[SE] == NULL) {
        i = quadrant(root, root->x, root->y);
        if ((no = node_ini(q)) == NULL) {
            fprintf(stderr, "%s:%d: node_ini failed\n", __FILE__,
                    __LINE__);
            return 1;
        }
        no->m = root->m;
        no->x = root->x;
        no->y = root->y;
        no->id = root->id;
        box(root, i, &no->xl, &no->xh, &no->yl, &no->yh);
        root->elm[i] = no;
        insert(q, root, x, y, m, id);
    } else {
        i = quadrant(root, x, y);
        if (root->elm[i] == NULL) {
            if ((no = node_ini(q)) == NULL) {
                fprintf(stderr, "%s:%d: node_ini failed\n", __FILE__,
                        __LINE__);
                return 1;
            }
            no->m = m;
            no->x = x;
            no->y = y;
            no->id = id;
            box(root, i, &no->xl, &no->xh, &no->yl, &no->yh);
            root->elm[i] = no;
        } else
            insert(q, root->elm[i], x, y, m, id);
        root->x = (root->x + x * m) / (root->m + m);
        root->y = (root->y + y * m) / (root->m + m);
    }
    return 0;
}

static int
walk(struct BarnesHut *q, struct Node *node, int (*fun)(struct BarnesHut *, struct Node *, void *data), void *data)
{
    if (node != NULL) {
        if (fun(q, node, data) != 0)
            return 1;
        if (walk(q, node->elm[NE], fun, data) != 0)
            return 1;
        if (walk(q, node->elm[SE], fun, data) != 0)
            return 1;
        if (walk(q, node->elm[NW], fun, data) != 0)
            return 1;
        if (walk(q, node->elm[SW], fun, data) != 0)
            return 1;        
    }
    return 0;
}
