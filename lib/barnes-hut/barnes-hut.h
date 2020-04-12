struct BarnesHut;
struct BarnesHut *barnes_hut_ini(double, double, double);
int barnes_hut_fin(struct BarnesHut *);
int barnes_hut_insert(struct BarnesHut *, double, double, double, long);
int barnes_hut_print(struct BarnesHut *, FILE *);
