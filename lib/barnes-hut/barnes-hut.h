struct BarnesHut;

struct BarnesHut *barnes_hut_ini(double, double, double);
struct BarnesHut *barnes_hut_build(long, const double *, const double *, const double *);

int barnes_hut_fin(struct BarnesHut *);
int barnes_hut_insert(struct BarnesHut *, double, double, double, long);

int barnes_hut_print(struct BarnesHut *, FILE *);
int barnes_hut_interaction(struct BarnesHut *q, double theta,
                           long id, double x, double y, long *cnt,
                           double *px, double *py, double *m);
struct BarnesHutInfo {
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
int barnes_hut_info(struct BarnesHut *q, double theta,
                    long id, double x, double y, long *cnt,
                    struct BarnesHutInfo *info);
