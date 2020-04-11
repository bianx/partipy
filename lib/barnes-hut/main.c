#include <stdlib.h>
#include <stdio.h>
#include "barnes-hut.h"

enum { NE, NW, SW, SE };
struct Node {
    double m;
    double x;
    double y;
    long cnt;
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
insert(struct BarnesHut *q, struct Node *root, double x, double y,
       double m, long id)
{
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
        no->cnt = 1;
        no->id = 1;
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
