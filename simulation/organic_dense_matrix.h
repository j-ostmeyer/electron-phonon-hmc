#ifndef DENSE_MATRIX
#define DENSE_MATRIX

void mat_mul(double *m, double *x, double *y, unsigned n);
void mat_mul_tr(double *m, double *x, double *y, unsigned n);
void mat_mat_mul(double *m1, double *m2, double *y, unsigned n);

void construct_id(double *m, unsigned n);
void transpose(double *m, unsigned n);
void transpose_rect(double *x, double *y, unsigned n1, unsigned n2);

double trace(double *m, unsigned n);

void apply_mu_K_to_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn);
void apply_mu_K_to_M(double *x, double *M, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn);

void set_taylor_pair_order();

void apply_K_plus_1_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn);
void apply_2nd_order_mu_K_to_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, double complex a);
void apply_pair_expansion_mu_K_to_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int revert, int stop);
void apply_exp_mu_K_to_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int revert, int stop);
void apply_exp_mu_K_to_M(double *x, double *M, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int revert, int stop);

void construct_hamilton(double *x, double *ham, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn);

void apply_current_to_v(double *x, double *v, double *w, double *t, double *kappa, unsigned *nnt, unsigned ns, unsigned nn, unsigned k);

void construct_M_hat(double *x, double *M_hat, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, organic_flags *mode);
void construct_QM1P(double *x, double *M_hat, double *QM1P, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, int time, organic_flags *mode);
void construct_all_QM1P(double *x, double *M_hat, double *QM1P, double *w, int *perm, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, organic_flags *mode);
void refine_QM1P(double *x, double *QM1P, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int time);
void refine_QM1P_approx(double *x, double *QM1P, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int time);

void cc_correlator(double *x, double *M_hat, double *R_t0, double *R_0t, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, double *greens, double logDetM, organic_flags *mode);
void phonon_correlator(double *x, unsigned ns, unsigned nn, unsigned nt, double *greens);

double logDetM_and_inv_M_hat(double *M_hat, double *w, int *p, unsigned ns, unsigned flavors, organic_flags *mode);

#endif
