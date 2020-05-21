#include <stdlib.h>
#include <stdio.h>
#include "tree.h"

enum { NE, NW, SW, SE };
static const int Dir[] = { NE, NW, SW, SE };

struct Particle {
    double x;
    double y;
    double m;
    long id;
};

struct Node {
    double w;
    double x;
    double y;
    long cnt;
    struct Node *elm[4];
    struct Particle *particle;
    struct Node *mother;
};

struct Tree {
    double x;
    double y;
    double w;
    long cap;
    long cnt;
    struct Node **node;
    struct TreeParam param;
};

static int quadrant(struct Node *, double, double);
static int box(struct Node *, int, double *, double *, double *);
static struct Node *node_ini(struct Tree *);
static int node_fin(struct Node *);
static int insert(struct Tree *, struct Node *, struct Particle *);
static int leaf(struct Node *);
static int insert(struct Tree *, struct Node *, struct Particle *);
static int insert0(struct Tree *, struct Node *, struct Particle *);
static int walk(struct Tree *, struct Node *,
                int (*)(struct Tree *, struct Node *, void *), void *);
static int print(struct Tree *, struct Node *, void *);
static struct Node *find_node(struct Node *, double, double);

struct Tree *
tree_ini(const struct TreeParam param, double x, double y, double w)
{
    struct Tree *q;

    if (param.cap <= 0) {
        fprintf(stderr, "%s:%d: param.cap <= 0\n", __FILE__, __LINE__);
        return NULL;
    }
    if ((q = malloc(sizeof(*q))) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return NULL;
    }
    q->param = param;
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

struct Tree *
tree_build(struct TreeParam param, long n, const double *x,
           const double *y, const double *m)
{
    double w;
    double xc;
    double xh;
    double xl;
    double yc;
    double yh;
    double yl;
    long i;
    struct Tree *q;

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
    xc = (xl + xh) / 2;
    yc = (yl + yh) / 2;

    if ((q = tree_ini(param, xc, yc, w)) == NULL) {
        fprintf(stderr, "%s:%d: tree_ini failed\n", __FILE__, __LINE__);
        return NULL;
    }

    for (i = 0; i < n; i++)
        if (tree_insert(q, x[i], y[i], m[i], i) != 0) {
            fprintf(stderr, "%s:%d: tree_insert failed\n", __FILE__,
                    __LINE__);
            return NULL;
        }
    return q;
}

int
tree_fin(struct Tree *q)
{
    long i;

    for (i = 0; i < q->cnt; i++)
        node_fin(q->node[i]);
    free(q->node);
    free(q);
    return 0;
}

int
tree_insert(struct Tree *q, double x, double y, double m, long id)
{
    struct Node *no;
    struct Particle particle;

    if (q->cnt == 0) {
        if ((no = node_ini(q)) == NULL) {
            fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
            return 1;
        }
        no->mother = NULL;
        no->x = q->x;
        no->y = q->y;
        no->w = q->w;
    }
    particle.x = x;
    particle.y = y;
    particle.m = m;
    particle.id = id;
    insert(q, q->node[0], &particle);
    return 0;
}

int
tree_print(struct Tree *q, FILE * f)
{
    if (q->cnt > 0)
        return walk(q, q->node[0], print, f);
    else
        return 0;
}

struct Node*
tree_node(struct Tree *q, double x, double y)
{
    return q->cnt == 0 ? NULL : find_node(q->node[0], x, y);
}

int
node_print(struct Node *q, FILE *f)
{
    double xl;
    double xh;
    double yl;
    double yh;
    if (q == NULL)
        return 1;
    xl = q->x - q->w / 2;
    xh = q->x + q->w / 2;
    yl = q->y - q->w / 2;
    yh = q->y + q->w / 2;
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
    return 0;
}

struct Node *
find_node(struct Node *q, double x, double y)
{
    int i;

    i = quadrant(q, x, y);
    if (q->elm[i] == NULL) {
        return q;
    } else
        return find_node(q->elm[i], x, y);
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
node_ini(struct Tree *q)
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
    if ((no->particle =
         malloc(q->param.cap * sizeof *no->particle)) == NULL) {
        fprintf(stderr, "%s:%d: malloc failed\n", __FILE__, __LINE__);
        return NULL;
    }
    no->cnt = 0;
    no->elm[NE] = NULL;
    no->elm[NW] = NULL;
    no->elm[SW] = NULL;
    no->elm[SE] = NULL;
    q->node[q->cnt++] = no;
    return no;
}

static int
node_fin(struct Node *q)
{
    free(q->particle);
    free(q);
    return 0;
}

static int
insert0(struct Tree *q, struct Node *root, struct Particle *particle)
{
    struct Node *no;
    int i;

    i = quadrant(root, particle->x, particle->y);
    if (root->elm[i] == NULL) {
        if ((no = node_ini(q)) == NULL) {
            fprintf(stderr, "%s:%d: node_ini failed\n", __FILE__,
                    __LINE__);
            return 1;
        }
        no->mother = root;
        box(root, i, &no->x, &no->y, &no->w);
        root->elm[i] = no;
    }
    return insert(q, root->elm[i], particle);
}

static int
insert(struct Tree *q, struct Node *root, struct Particle *particle)
{
    long k;

    if (leaf(root)) {
        if (root->cnt < q->param.cap) {
            root->particle[root->cnt] = *particle;
            root->cnt++;
        } else {
            for (k = 0; k < root->cnt; k++) {
                if (insert0(q, root, &root->particle[k]) != 0) {
                    fprintf(stderr, "%s:%d: insert0 failed\n", __FILE__,
                            __LINE__);
                    return 1;
                }
            }
            root->cnt = 0;
            if (insert0(q, root, particle) != 0) {
                fprintf(stderr, "%s:%d: insert0 failed\n", __FILE__,
                        __LINE__);
                return 1;
            }
        }
    } else if (insert0(q, root, particle) != 0) {
        fprintf(stderr, "%s:%d: insert0 failed\n", __FILE__, __LINE__);
        return 1;
    }
    return 0;
}

static int
leaf(struct Node *q)
{
    return q->elm[NE] == NULL && q->elm[NW] == NULL
        && q->elm[SW] == NULL && q->elm[SE] == NULL;
}

static int
walk(struct Tree *q, struct Node *node,
     int (*fun)(struct Tree *, struct Node *, void *data), void *data)
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

static int
print(struct Tree *q, struct Node *n, void *f0)
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
