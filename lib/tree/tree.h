struct Tree;
struct TreeParam {
    long cap;
};
struct Tree *tree_ini(struct TreeParam, double, double, double);
struct Tree *tree_build(struct TreeParam, long, const double *, const double *, const double *);

int tree_fin(struct Tree *);
int tree_insert(struct Tree *, double, double, double, long);

int tree_print(struct Tree *, FILE *);
int tree_interaction(struct Tree *q, double theta,
                           long id, double x, double y, long *cnt,
                           double *px, double *py, double *m);
struct TreeInfo {
    double m;
    double mx;
    double my;
    double w;
    double x;
    double y;
    int Leaf;
    int Coarse;
    long id;
};
int tree_info(struct Tree *q, double theta,
                    long id, double x, double y, long *cnt,
                    struct TreeInfo *info);
