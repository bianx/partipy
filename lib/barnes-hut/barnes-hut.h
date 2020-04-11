struct BarnesHut;
struct BarnesHut *barnes_hut_ini(double, double, double, double);
int barnes_hut_fin(struct BarnesHut *);
int barnes_hut_insert(struct BarnesHut *, double, double, double, long);
int barnes_hut_build(struct BarnesHut *, long, const double *,
                     const double *);
int barnes_hut_rectangle(struct BarnesHut *, double xlo, double xhi,
                         double ylo, double yhi, long *, const long **);
