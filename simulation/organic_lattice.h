#ifndef ORGANIC_LATTICE
#define ORGANIC_LATTICE

unsigned *construct_triangular(unsigned l1, unsigned l2, unsigned *ns, unsigned *nn);

double *construct_contraction_mat(double *a, double *b, unsigned n);

double contr_gr_mat(double *gr, double *mat, unsigned n);
double contr_gr_triangle(double *gr, unsigned n);
double contract_greens(double *gr, double *mat, unsigned nn, unsigned tau);

#endif
