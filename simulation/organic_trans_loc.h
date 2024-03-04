#ifndef ORGANIC_TRANS_LOC
#define ORGANIC_TRANS_LOC

void hist_trans_loc(double *x, double *t, double *kappa, double mu, double *contr_mat, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, double width, double *results, double *greens, FILE *res_out);

unsigned run_trans_loc(double *greens, unsigned l1, unsigned l2, unsigned nt, double *t, double *kappa, double mu, double *contr_mat, unsigned meas, const char *res_name, const char *corr_name);

#endif
