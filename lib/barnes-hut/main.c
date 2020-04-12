#include <stdlib.h>
#include <stdio.h>
#include "barnes-hut.h"

enum { NE, NW, SW, SE };
static const int Dir[] = { NE, NW, SW, SE };

struct Node {
    double m;
    double mx;
    double my;
    double w;
    double x;
    double y;
    long id;
    struct Node *elm[4];
};

struct BarnesHut {
    double x;
    double y;
    double w;
    long cap;
    long cnt;
    struct Node **node;
};

struct InfoParam {
    struct BarnesHutInfo *info;
    long cnt;
};

struct Param {
    double *x;
    double *y;
    double *m;
    long cnt;
};

static int trace(struct Node *, int, void *);
static int info(struct Node *, int, void *);
static int quadrant(struct Node *, double, double);
static int box(struct Node *, int, double *, double *, double *);
static struct Node *node_ini(struct BarnesHut *);
static int insert(struct BarnesHut *, struct Node *, double, double,
                  double, long);
static int walk(struct BarnesHut *, struct Node *,
                int (*)(struct BarnesHut *, struct Node *, void *),
                void *);

struct BarnesHut *
barnes_hut_ini(double x, double y, double w)
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
    q->x = x;
    q->y = y;
    q->w = w;
    return q;
}

static double
max(double a, double b)
{
    return a > b ? a : b;
}

static double
min(double a, double b)
{
    return a < b ? a : b;
}

struct BarnesHut *
barnes_hut_build(long n, const double *x, const double *y, const double *m)
{
    double w;
    double xc;
    double xh;
    double xl;
    double yc;
    double yh;
    double yl;
    long i;
    struct BarnesHut *q;

    if (n < 1) {
        fprintf(stderr, "%s:%d: n < 1\n", __FILE__, __LINE__);
        return NULL;
    }

    xl = xh = x[0];
    yl = yh = y[0];
    for (i = 1; i < n; i++) {
        xl = min(xl, x[i]);
        xh = max(xh, x[i]);
        yl = min(yl, y[i]);
        yh = max(yh, y[i]);        
    }
    w = max(xh - xl, yh - yl);
    xc = (xl + xh)/2;
    yc = (yl + yh)/2;

    if ((q = barnes_hut_ini(xc, yc, w)) == NULL) {
        fprintf(stderr, "%s:%d: barnes_hut_ini failed\n", __FILE__, __LINE__);
        return NULL;
    }

    for (i = 0; i < n; i++)
        if (barnes_hut_insert(q, x[i], y[i], m[i], i) != 0) {
            fprintf(stderr, "%s:%d: barnes_hut_insert failed\n", __FILE__, __LINE__);
            return NULL;            
        }
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
        no->mx = m * x;
        no->my = m * y;
        no->id = id;
        no->x = q->x;
        no->y = q->y;
        no->w = q->w;
        no->elm[NE] = NULL;
        no->elm[NW] = NULL;
        no->elm[SW] = NULL;
        no->elm[SE] = NULL;
    } else {
        insert(q, q->node[0], x, y, m, id);
    }
    return 0;
}

static double
sq(double x)
{
    return x * x;
}

static int
force(struct Node *n, double theta, int (*fun)(struct Node *, int, void *),
      long id, double x, double y, void *data)
{
    double r2;
    double w2;
    int i;
    int j;
    int Coarse;

    r2 = sq(x - n->mx / n->m) + sq(y - n->my / n->m);
    w2 = n->w * n->w;
    Coarse = w2 < theta * theta * r2;
    if (n->elm[NE] == NULL && n->elm[NW] == NULL && n->elm[SW] == NULL
        && n->elm[SE] == NULL) {
        if (id != n->id) {
            if (fun(n, Coarse, data)) {
                fprintf(stderr, "%s:%d: function() failed\n", __FILE__,
                        __LINE__);
                return 1;
            }
        }
    } else {
        if (Coarse) {
            if (fun(n, Coarse, data)) {
                fprintf(stderr, "%s:%d: function() failed\n", __FILE__,
                        __LINE__);
                return 1;
            }
        } else {
            for (i = 0; i < (int) (sizeof Dir / sizeof *Dir); i++) {
                j = Dir[i];
                if (n->elm[j] != NULL) {
                    if (force(n->elm[j], theta, fun, id, x, y, data)) {
                        fprintf(stderr, "%s:%d: function() failed\n",
                                __FILE__, __LINE__);
                        return 1;
                    }
                }
            }
        }
    }
    return 0;
}


int
barnes_hut_interaction(struct BarnesHut *q, double theta,
                       long id, double x, double y, long *cnt, double *px,
                       double *py, double *m)
{
    struct Param param;

    param.cnt = 0;
    param.x = px;
    param.y = py;
    param.m = m;
    if (q->cnt > 0) {
        if (force(q->node[0], theta, trace, id, x, y, &param) != 0) {
            fprintf(stderr, "%s:%d: force() failed\n", __FILE__, __LINE__);
            return 1;
        }
    }
    *cnt = param.cnt;
    return 0;
}

int
barnes_hut_info(struct BarnesHut *q, double theta,
                long id, double x, double y, long *cnt,
                struct BarnesHutInfo *info0)
{
    struct InfoParam param;

    param.cnt = 0;
    param.info = info0;
    if (q->cnt > 0) {
        if (force(q->node[0], theta, info, id, x, y, &param) != 0) {
            fprintf(stderr, "%s:%d: force() failed\n", __FILE__, __LINE__);
            return 1;
        }
    }
    *cnt = param.cnt;
    return 0;
}

static int
info(struct Node *n, int coarse, void *p0)
{
    struct BarnesHutInfo info;
    struct InfoParam *p;

    p = p0;
    info.id = n->id;
    info.m = n->m;
    info.mx = n->mx;
    info.my = n->my;
    info.w = n->w;
    info.x = n->x;
    info.y = n->y;
    info.Leaf = n->elm[NE] == NULL && n->elm[NW] == NULL
        && n->elm[SW] == NULL && n->elm[SE] == NULL;
    info.Coarse = coarse;
    p->info[p->cnt++] = info;
    return 0;
}

static int
trace(struct Node *n, int coarse, void *p0)
{
    struct Param *p;
    double m;

    (void) coarse;
    p = p0;
    m = n->m;
    p->x[p->cnt] = m != 0 ? n->mx / m : 0;
    p->y[p->cnt] = m != 0 ? n->my / m : 0;
    p->m[p->cnt] = m;
    p->cnt++;
    return 0;
}

static int
print(struct BarnesHut *q, struct Node *n, void *f0)
{
    FILE *f;
    double xl;
    double xh;
    double yl;
    double yh;

    (void) q;
    f = f0;
    xl = n->x - n->w / 2;
    xh = n->x + n->w / 2;
    yl = n->y - n->w / 2;
    yh = n->y + n->w / 2;
    if (n->elm[NE] == NULL && n->elm[NW] == NULL
        && n->elm[SW] == NULL && n->elm[SE] == NULL) {
        if (fprintf(f, "%.16g %.16g\n", xl, yl) < 0)
            return 1;
        if (fprintf(f, "%.16g %.16g\n", xh, yl) < 0)
            return 1;
        if (fprintf(f, "%.16g %.16g\n", xh, yh) < 0)
            return 1;
        if (fprintf(f, "%.16g %.16g\n", xl, yh) < 0)
            return 1;
        if (fprintf(f, "%.16g %.16g\n", xl, yl) < 0)
            return 1;
        if (fprintf(f, "\n") < 0)
            return 1;
    }
    return 0;
}

int
barnes_hut_print(struct BarnesHut *q, FILE * f)
{
    if (q->cnt > 0)
        return walk(q, q->node[0], print, f);
    else
        return 0;
}

static int
quadrant(struct Node *n, double x, double y)
{
    return (x > n->x) ? (y > n->y) ? NE : SE : (y > n->y) ? NW : SW;
}

static int
box(struct Node *n, int i, double *px, double *py, double *pw)
{
    double l;

    *px = n->x;
    *py = n->y;
    *pw = n->w / 2;
    l = n->w / 4;
    switch (i) {
    case NE:
        *px += l;
        *py += l;
        break;
    case SE:
        *px += l;
        *py -= l;
        break;
    case NW:
        *px -= l;
        *py += l;
        break;
    case SW:
        *px -= l;
        *py -= l;
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
    double xc;
    double yc;

    if (root->elm[NE] == NULL && root->elm[NW] == NULL
        && root->elm[SW] == NULL && root->elm[SE] == NULL) {
        xc = root->m == 0 ? root->mx : root->mx / root->m;
        yc = root->m == 0 ? root->my : root->my / root->m;
        i = quadrant(root, xc, yc);
        if ((no = node_ini(q)) == NULL) {
            fprintf(stderr, "%s:%d: node_ini failed\n", __FILE__,
                    __LINE__);
            return 1;
        }
        no->m = root->m;
        no->mx = root->mx;
        no->my = root->my;
        no->id = root->id;
        box(root, i, &no->x, &no->y, &no->w);
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
            no->mx = m * x;
            no->my = m * y;
            no->id = id;
            box(root, i, &no->x, &no->y, &no->w);
            root->elm[i] = no;
        } else
            insert(q, root->elm[i], x, y, m, id);
        root->mx += m * x;
        root->my += m * y;
        root->m += m;
    }
    return 0;
}

static int
walk(struct BarnesHut *q, struct Node *node,
     int (*fun)(struct BarnesHut *, struct Node *, void *data), void *data)
{
    int i;
    int j;

    if (node != NULL) {
        if (fun(q, node, data) != 0)
            return 1;
        for (i = 0; i < (int) (sizeof Dir / sizeof *Dir); i++) {
            j = Dir[i];
            if (walk(q, node->elm[j], fun, data) != 0)
                return 1;
        }
    }
    return 0;
}
