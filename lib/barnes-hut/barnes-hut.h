struct Quadtree;
struct Quadtree *quadtree_ini(void);
int quadtree_fin(struct Quadtree *);
int quadtree_insert(struct Quadtree *, double, double, long);
int quadtree_build(struct Quadtree *, long, const double *,
                   const double *);
int quadtree_rectangle(struct Quadtree *, double xlo, double xhi,
                       double ylo, double yhi, long *, const long **);
