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
box(int i, struct Node *n, double *xl, double *xh, double *yl, double *yh)
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

static int
insert(struct BarnesHut *q, struct Node *root, double x, double y,
       double m, long id)
{
    int i;
    struct Node *no;

    if ((no = malloc(sizeof *no)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return 1;
    }
    if (q->cnt == q->cap) {
        q->cap *= 2;
        if ((q->node = realloc(q->node, q->cap * sizeof *q->node)) == NULL) {
            fprintf(stderr, "%s:%d: realloc failed\n", __FILE__, __LINE__);
            return 1;
        }
    }
    no->elm[NE] = NULL;
    no->elm[NW] = NULL;
    no->elm[SW] = NULL;
    no->elm[SE] = NULL;
    no->m = m;
    no->id = id;
    q->node[q->cnt++] = no;
    if (root->elm[NE] == NULL && root->elm[NW] == NULL
        && root->elm[SW] == NULL && root->elm[SE] == NULL) {
        i = quadrant(root, root->x, root->y);
        root->elm[i] = no;
        box(i, root, &no->xl, &no->xh, &no->yl, &no->yh);
        /* */
    } else {
    }
    return 0;
}

int
barnes_hut_insert(struct BarnesHut *q, double x, double y, double m,
                  long id)
{
    double xh;
    double xl;
    double yh;
    double yl;
    long cap;
    long cnt;
    struct Node **node;
    struct Node *no;

    xl = q->xl;
    xh = q->xh;
    yl = q->yl;
    yh = q->yh;
    cap = q->cap;
    cnt = q->cnt;
    node = q->node;
    if (cnt == 0) {
        if ((no = malloc(sizeof *no)) == NULL) {
            fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
            return 1;
        }
        no->m = m;
        no->x = x;
        no->y = y;
        no->id = id;
        no->xl = xl;
        no->xh = xh;
        no->yl = yl;
        no->yh = yh;
        no->elm[NE] = NULL;
        no->elm[NW] = NULL;
        no->elm[SW] = NULL;
        no->elm[SE] = NULL;
        node[cnt++] = no;
    } else
        insert(q, node[0], x, y, m, id);
    q->cap = cap;
    q->cnt = cnt;
    q->node = node;
    return 0;
}
