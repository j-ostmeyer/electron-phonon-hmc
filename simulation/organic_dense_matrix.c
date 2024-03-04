#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <omp.h>

#ifndef MKL_LAPACKE
#include <cblas.h>
#include <lapacke.h>
#else
#include <mkl.h>
#include <mkl_lapacke.h>
#endif

#include "organic_flags.h"
#include "organic_aux.h"
#include "organic_lattice.h"

#include "organic_dense_matrix.h"

#include "taylor-coeff_pairs.h" // defines the pairs in the Taylor product expansion

// ATTENTION: dense matrix operators not parallelised!
void mat_mul(double *m, double *x, double *y, unsigned n){
	for(unsigned i = 0; i < n; i++){
		y[i] = 0;
		for(unsigned k = 0; k < n; k++, m++) y[i] += (*m) * x[k];
	}
}

void mat_mul_tr(double *m, double *x, double *y, unsigned n){
	for(unsigned k = 0; k < n; k++) y[k] = 0;
	for(unsigned i = 0; i < n; i++)
		for(unsigned k = 0; k < n; k++, m++) y[k] += (*m) * x[i];
}

void mat_mat_mul(double *m1, double *m2, double *y, unsigned n){
	for(unsigned i = 0; i < n; i++, m2 += n){
		mat_mul(m1, m2, y, n);
		for(unsigned j = 0; j < n; j++) m2[j] = y[j];
	}
}

void construct_id(double *m, unsigned n){
	// construct identity
	const unsigned n2 = n*n;
	const unsigned n1 = n+1;
#pragma omp parallel for
	for(unsigned i = 0; i < n2; i++) m[i] = 0;
//#pragma omp parallel for
	for(unsigned i = 0; i < n; i++) m[i*n1] = 1;
}

void transpose(double *m, unsigned n){
#pragma omp parallel for //schedule(dynamic)
	for(unsigned i = 0; i < n; i++){
		const unsigned shift = i*n;

		for(unsigned j = i+1; j < n; j++){
			const unsigned pos = shift + j;
			const unsigned posT = j*n + i;

			const double tmp = m[pos];
			m[pos] = m[posT];
			m[posT] = tmp;
		}
	}
}

void transpose_rect(double *x, double *y, unsigned n1, unsigned n2){
#pragma omp parallel for //schedule(dynamic)
	for(unsigned i = 0; i < n2; i++){
		const unsigned shift = i*n1;

		for(unsigned j = 0; j < n1; j++){
			const unsigned pos = shift + j;
			const unsigned posT = j*n2 + i;

			y[posT] = x[pos];
		}
	}
}

double trace(double *m, unsigned n){
	double tr = 0;
	const unsigned n1 = n+1;
//#pragma omp parallel for
	for(unsigned i = 0; i < n; i++) tr += m[i*n1];
	return tr;
}

//Applying the single-particle Hamiltonian to v, w is the output, T are hoppings, 
void apply_mu_K_to_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn){
	const unsigned nn2 = nn/2;

//#pragma omp parallel for
	for(unsigned i = 0; i < ns; i++) w[i] = mu * v[i];
	
//#pragma omp parallel for
	for(unsigned i = 0; i < ns; i++){
		const unsigned shift = i*nn, shift2 = i*nn2;

		for(unsigned k = 0; k < nn2; k++){
			const unsigned k2 = 2*k;
			const unsigned neighbour = nnt[shift+k2];
			const double fac = t[k] - kappa[k]*x[shift2+k];

			const double here = fac * v[neighbour];
			const double there = fac * v[i];

//#pragma omp atomic
			w[i] += here;
//#pragma omp atomic
			w[neighbour] += there;
		}
	}
}

//Multiplying M with K on the left
void apply_mu_K_to_M(double *x, double *M, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn){
#pragma omp parallel for
	for(unsigned i = 0; i < ns; i++){
		double *M_cp = M + i*ns;
		double *w_cp = w + i*ns;
		apply_mu_K_to_v(x, M_cp, w_cp, t, kappa, mu, nnt, ns, nn);
		for(unsigned j = 0; j < ns; j++) M_cp[j] = w_cp[j];
	}
}

void set_taylor_pair_order(){
	const unsigned n_tuples = NUM_TAYLOR_PAIRS/2;
	unsigned tuple;

	for(unsigned i = 0, down = 0, up = n_tuples-1; i < n_tuples; i++){
		// group coefficients into tuples of two pairs each and average stride close to 1
		// sort them so that the overall stride matches the actual position
		if(i % 3 == 1){
			tuple = down;
			down++;
		}else{
			tuple = up;
			up--;
		}

		for(unsigned j = 0; j < 2; j++){
			const unsigned pos = j? NUM_TAYLOR_PAIRS-1 - tuple : tuple;

			taylor_pair_order[2*i + j] = pos;
		}
	}

	if(NUM_TAYLOR_PAIRS % 2) taylor_pair_order[2*n_tuples] = n_tuples;

//	for(unsigned i = 0; i < NUM_TAYLOR_PAIRS; i++)
//		printf("%d:\t%g + %g I\n", taylor_pair_order[i], creal(taylor_pair_coeffs[taylor_pair_order[i]]), cimag(taylor_pair_coeffs[taylor_pair_order[i]]));
}

// 1 + h dt
void apply_K_plus_1_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn){
	apply_mu_K_to_v(x, v, w, t, kappa, mu, nnt, ns, nn);
//#pragma omp parallel for
	for(unsigned i = 0; i < ns; i++) v[i] += w[i];
}

void apply_2nd_order_mu_K_to_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, double complex a){
	const double a_lin = 2*creal(a), a_quad = pow(cabs(a), 2);
	double *z = w + ns;

	apply_mu_K_to_v(x, v, w, t, kappa, mu, nnt, ns, nn);
	if(cimag(a))
		apply_mu_K_to_v(x, w, z, t, kappa, mu, nnt, ns, nn);

	if(cimag(a)){
		//#pragma omp parallel for
		for(unsigned i = 0; i < ns; i++) v[i] += a_lin * w[i] + a_quad * z[i];
	}else{
		//#pragma omp parallel for
		for(unsigned i = 0; i < ns; i++) v[i] += a_lin * w[i];
	}
}


void apply_pair_expansion_mu_K_to_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int revert, int stop){
	if(revert){
		for(int i = NUM_TAYLOR_PAIRS-1; i > stop; i--)
			apply_2nd_order_mu_K_to_v(x, v, w, t, kappa, mu, nnt, ns, nn, .5*taylor_pair_coeffs[taylor_pair_order[i]]/NUM_TAYLOR_PAIRS);
	}else{
		for(int i = 0; i < stop; i++)
			apply_2nd_order_mu_K_to_v(x, v, w, t, kappa, mu, nnt, ns, nn, .5*taylor_pair_coeffs[taylor_pair_order[i]]/NUM_TAYLOR_PAIRS);
	}
}

//Is this exp(h*dt)???
void apply_exp_mu_K_to_v(double *x, double *v, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int revert, int stop){
	// calculating exp(mu) exactly
	// increases precision for small Nt, but decreases acceptance
	//const double exp_mu = exp(mu);
	//for(unsigned i = 0; i < ns; i++) v[i] *= exp_mu;
	//apply_pair_expansion_mu_K_to_v(x, v, w, t, kappa, 0, nnt, ns, nn, revert, stop);
	
	// treat mu and K in the same way in the expansion
	// better approximation of the largest and most important contributions
	apply_pair_expansion_mu_K_to_v(x, v, w, t, kappa, mu, nnt, ns, nn, revert, stop);
}

void apply_exp_mu_K_to_M(double *x, double *M, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int revert, int stop){
#pragma omp parallel for
	for(unsigned i = 0; i < ns; i++){
		double *M_cp = M + i*ns;
		double *w_cp = w + 2*i*ns;
		apply_exp_mu_K_to_v(x, M_cp, w_cp, t, kappa, mu, nnt, ns, nn, revert, stop);
	}
}

void construct_hamilton(double *x, double *ham, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn){
	// construct identity
	construct_id(ham, ns);

	// assuming only one time slice
	apply_mu_K_to_M(x, ham, w, t, kappa, mu, nnt, ns, nn);
}

void apply_current_to_v(double *x, double *v, double *w, double *t, double *kappa, unsigned *nnt, unsigned ns, unsigned nn, unsigned k){
	// apply current operator I_k (in direction k) for all sites to vector
	const unsigned nn2 = nn/2;

//#pragma omp parallel for
	for(unsigned i = 0; i < ns; i++) w[i] = 0;
	
//#pragma omp parallel for
	for(unsigned i = 0; i < ns; i++){
		const unsigned shift = i*nn, shift2 = i*nn2;
		const unsigned k2 = 2*k;
		const unsigned neighbour = nnt[shift+k2];
		const double fac = t[k] - kappa[k]*x[shift2+k];

		// anti-symmetric version of hopping
		const double here  =  fac * v[neighbour];
		const double there = -fac * v[i];

		//#pragma omp atomic
		w[i] += here;
		//#pragma omp atomic
		w[neighbour] += there;
	}
}

void construct_M_hat(double *x, double *M_hat, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, organic_flags *mode){
	// construct identity
	construct_id(M_hat, ns);

	// multiply all M_t from right to left starting with the last one
	for(int tau = nt-1; tau >= 0; tau--){
		apply_exp_mu_K_to_M(x + tau*ns*nn/2, M_hat, w, t, kappa, mu, nnt, ns, nn, 0, NUM_TAYLOR_PAIRS);
	}

	if(!mode->zero_filling && !mode->single_particle){
		// add one
		//#pragma omp parallel for
		for(unsigned i = 0; i < ns; i++) M_hat[i*(ns+1)] += 1;
	}
}

void construct_QM1P(double *x, double *M_hat, double *QM1P, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, int time, organic_flags *mode){
	const unsigned mat_dim = ns*ns;

	if(mode->schur_form_inv){
		// exact analytic derivative using exponential for 0th time slice
		if(!mode->single_particle){
#pragma omp parallel for
			for(unsigned i = 0; i < mat_dim; i++) QM1P[i] = -M_hat[i];
			for(unsigned i = 0; i < ns; i++) QM1P[i*(ns+1)] += 1;
		}else{
#pragma omp parallel for
			for(unsigned i = 0; i < mat_dim; i++) QM1P[i] = M_hat[i];
		}
	}else{
		// exact derivative of polynomially approximated exponential for given time slice
		if(!mode->single_particle){
			// start inside out with M_hat^-1
#pragma omp parallel for
			for(unsigned i = 0; i < mat_dim; i++) QM1P[i] = M_hat[i];
		}else
			construct_id(QM1P, ns);

		// multiply right to left M_{t+1} *...* M_{Nt-1} * M_hat^-1
		for(int tau = nt-1; tau > time; tau--){
			apply_exp_mu_K_to_M(x + tau*ns*nn/2, QM1P, w, t, kappa, mu, nnt, ns, nn, 0, NUM_TAYLOR_PAIRS);
		}

		// transpose (Q_t * M_hat^-1), so P_t^t can be applied from the left
		transpose(QM1P, ns);

		// multiply right to left (Q_t * M_hat^-1 * P_t)^t
		// = (M_0 * M_1 *...* M_{t-1})^t * (Q_t * M_hat^-1)^t
		// = M_{t-1} *...* M_1 * M_0 * (Q_t * M_hat^-1)^t
		for(int tau = 0; tau < time; tau++){
			apply_exp_mu_K_to_M(x + tau*ns*nn/2, QM1P, w, t, kappa, mu, nnt, ns, nn, 0, NUM_TAYLOR_PAIRS);
		}

		// revert transposition
		transpose(QM1P, ns);

		refine_QM1P(x, QM1P, w, t, kappa, mu, nnt, ns, nn, time);
	}
}

void construct_all_QM1P(double *x, double *M_hat, double *QM1P, double *w, int *perm, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, organic_flags *mode){
	const unsigned mat_dim = ns*ns;
	const unsigned iterations = ceil_log2(nt);

	double *QM1P_left, *QM1P_right;
	unsigned chunk, remaining = nt-1;

	if(!mode->single_particle){
		// start inside out with M_hat^-1
#pragma omp parallel for
		for(unsigned i = 0; i < mat_dim; i++) QM1P[i] = M_hat[i];
	}else
		construct_id(QM1P, ns);

	// store number of matrices already accumulated
	// from the left (first half) and the right (second half)
	for(unsigned i = 0; i < 2*nt; i++) perm[i] = 0;

	for(unsigned r = 0; r < iterations; r++){
		chunk = ceil_half(remaining);

		for(int i = nt-1; i >= 0; i--){ // going backward to avoid re-iteration of same entry
			if(i == 0 || perm[i] || perm[nt+i]){
				QM1P_left = QM1P + i*mat_dim;

				if(i + chunk >= nt) printf("ERROR: Try to access time-slice %d out of %d!\n", i + chunk, nt);

				if(perm[i+chunk] == 0 && perm[nt+i+chunk] == 0){
					// copy all the information from the previous block
					QM1P_right = QM1P + (i+chunk)*mat_dim;
#pragma omp parallel for
					for(unsigned j = 0; j < mat_dim; j++) QM1P_right[j] = QM1P_left[j];

					perm[i+chunk] = perm[i];
					perm[nt+i+chunk] = perm[nt+i];

					// apply the chunk from the right

					// transpose (Q_t * M_hat^-1), so P_t^t can be applied from the left
					transpose(QM1P_right, ns);

					// multiply right to left (Q_t * M_hat^-1 * P_t)^t
					// = (M_0 * M_1 *...* M_{t-1})^t * (Q_t * M_hat^-1)^t
					// = M_{t-1} *...* M_1 * M_0 * (Q_t * M_hat^-1)^t
					for(int tau = perm[nt+i+chunk]; tau < perm[nt+i+chunk] + chunk; tau++){
						apply_exp_mu_K_to_M(x + tau*ns*nn/2, QM1P_right, w, t, kappa, mu, nnt, ns, nn, 0, NUM_TAYLOR_PAIRS);
					}

					// revert transposition
					transpose(QM1P_right, ns);

					perm[nt+i+chunk] += chunk;
				}

				// multiply right to left M_{t+1} *...* M_{Nt-1} * M_hat^-1
				for(int tau = nt-1 - perm[i]; tau > nt-1 - perm[i] - chunk; tau--){
					apply_exp_mu_K_to_M(x + tau*ns*nn/2, QM1P_left, w, t, kappa, mu, nnt, ns, nn, 0, NUM_TAYLOR_PAIRS);
				}

				perm[i] += chunk;
			}
		}

		remaining -= chunk;
	}

	if(remaining) printf("ERROR: Construction of QM1P matrix has %d un-accounted time-slices!\n", remaining);

	for(unsigned i = 0; i < nt; i++)
		if(perm[i] + perm[nt+i] != nt-1) printf("ERROR: QM1P[%d] has %d + %d time-slices instead of %d!\n", i, perm[i], perm[nt+i], nt-1);
}

void refine_QM1P(double *x, double *QM1P, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int time){
#ifndef TAYLOR_ORDER_1
	const unsigned mat_dim = ns*ns;
	double *QM1P_cp = w + 2*mat_dim; // mat_dim
	double *QM1P_loc = QM1P_cp + mat_dim; // mat_dim
	double *QM1P_loc2 = QM1P_loc + mat_dim; // mat_dim

	// QM1P: input and final result, QM1P_cp: copy of input, QM1P_loc: each local term
#pragma omp parallel for
	for(unsigned i = 0; i < mat_dim; i++){
		QM1P_cp[i] = QM1P[i];
		QM1P[i] = 0;
	}

	//print_mat(QM1P_cp, ns);

	for(unsigned pair = 0; pair < NUM_TAYLOR_PAIRS; pair++){
		const unsigned pos = taylor_pair_order[pair];
		const double complex a = .5*taylor_pair_coeffs[pos] / NUM_TAYLOR_PAIRS;

		// start inside out with preliminary QM1P
#pragma omp parallel for
		for(unsigned i = 0; i < mat_dim; i++) QM1P_loc[i] = QM1P_cp[i];

		// multiply right to left
		apply_exp_mu_K_to_M(x + time*ns*nn/2, QM1P_loc, w, t, kappa, mu, nnt, ns, nn, 0, pair);

		// transpose, so the rest can be applied from the left
		transpose(QM1P_loc, ns);

		apply_exp_mu_K_to_M(x + time*ns*nn/2, QM1P_loc, w, t, kappa, mu, nnt, ns, nn, 1, pair);

		if(cimag(a)){
#pragma omp parallel for
			for(unsigned i = 0; i < mat_dim; i++) QM1P_loc2[i] = QM1P_loc[i];

			apply_mu_K_to_M(x + time*ns*nn/2, QM1P_loc2, w, t, kappa, mu, nnt, ns, nn);

			// revert transposition
			transpose(QM1P_loc2, ns);
		}
		transpose(QM1P_loc, ns);

		const double a_lin = 2*creal(a);
#pragma omp parallel for
		for(unsigned i = 0; i < mat_dim; i++) QM1P[i] += a_lin * QM1P_loc[i];

		if(cimag(a)){
			apply_mu_K_to_M(x + time*ns*nn/2, QM1P_loc, w, t, kappa, mu, nnt, ns, nn);

			const double a_quad = pow(cabs(a), 2);
#pragma omp parallel for
			for(unsigned i = 0; i < mat_dim; i++) QM1P[i] += a_quad * (QM1P_loc[i] + QM1P_loc2[i]);
		}

		//print_mat(QM1P_loc, ns);
		//print_mat(QM1P_loc2, ns);
	}
#endif
}

void refine_QM1P_approx(double *x, double *QM1P, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, int time){
	// Calculates the polynomial approximation to the derivative exp(mu + K)
	// instead of the exact derivative of the polynomial approximation.
	// Works well for high Taylor orders >= 6, but bad for lower orders.
	apply_mu_K_to_M(x + time*ns*nn/2, QM1P, w, t, kappa, mu, nnt, ns, nn);
}

void cc_correlator(double *x, double *M_hat, double *R_t0, double *R_0t, double *w, double *t, double *kappa, double mu, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, double *greens, double logDetM, organic_flags *mode){
	const unsigned nn2 = nn/2, mat_dim = ns*ns;

	if(!mode->zero_filling && !mode->single_particle){
		// building the transfer matrices 0 -> t, t -> 0
#pragma omp parallel for
		for(unsigned j = 0; j < mat_dim; j++){
			R_t0[j] = M_hat[j];
			R_0t[j] = M_hat[j];
		}
		transpose(R_0t, ns); // only ever use transposed version of R_0t!
		//for(unsigned i = 0; i < ns; i++) R_0t[i*(ns+1)] += 1;
	}else{
		construct_id(R_t0, ns);
		construct_id(R_0t, ns);
	}

	// multiply right to left M_{t} *...* M_{Nt-1} * M_hat^-1
	for(int tau = nt-1; tau >= 0; tau--){
		apply_exp_mu_K_to_M(x + tau*ns*nn2, R_t0, w, t, kappa, mu, nnt, ns, nn, 0, NUM_TAYLOR_PAIRS);
		if(tau){
			double *R_tmp = R_t0 + tau*mat_dim;
#pragma omp parallel for
			for(unsigned j = 0; j < mat_dim; j++) R_tmp[j] = R_t0[j];
		}
	}

	for(unsigned tau = 0; tau < nt; tau++, R_t0 += mat_dim){
		const unsigned shift = tau*ns*nn2;

		for(unsigned a = 0; a < nn2; a++){
			for(unsigned b = 0; b < nn2; b++){
				double g = 0;

#pragma omp parallel for reduction(+:g)
				for(unsigned i = 0; i < ns; i++){
					for(unsigned j = 0; j < ns; j++){

						// the different neighbours
						const unsigned ia = nnt[i*nn+2*a];
						const unsigned jb = nnt[j*nn+2*b];

						// the effective hoppings
						const double chi_ia = t[a] - kappa[a] * x[i*nn2 + a];
						const double chi_jb = t[b] - kappa[b] * x[shift + j*nn2 + b];

						double corr = R_t0[j*ns + i] * R_0t[jb*ns + ia];
						corr += R_t0[jb*ns + ia] * R_0t[j*ns + i];
						corr -= R_t0[jb*ns + i] * R_0t[j*ns + ia];
						corr -= R_t0[j*ns + ia] * R_0t[jb*ns + i];

						// disconnected part, appears to be irrelevant
						//corr -= M_hat[i*ns + ia] * M_hat[j*ns + jb];
						//corr -= M_hat[ia*ns + i] * M_hat[jb*ns + j];
						//corr += M_hat[i*ns + ia] * M_hat[jb*ns + j];
						//corr += M_hat[ia*ns + i] * M_hat[j*ns + jb];

						g += chi_ia*chi_jb * corr;
					}
				}

				const unsigned pos = tau*nn2*nn2 + a*nn2 + b;

				if(mode->single_particle)
					greens[pos] = g*exp(-logDetM);
				else
					greens[pos] = g/ns;
			}
		}

		if(tau < nt-1)
			apply_exp_mu_K_to_M(x + shift, R_0t, w, t, kappa, mu, nnt, ns, nn, 0, NUM_TAYLOR_PAIRS);
	}
}

void phonon_correlator(double *x, unsigned ns, unsigned nn, unsigned nt, double *greens){
	const unsigned nn2 = nn/2;

	for(unsigned tau = 0; tau < nt; tau++){
		const unsigned shift = tau*ns*nn2;

		for(unsigned a = 0; a < nn2; a++){
			for(unsigned b = 0; b < nn2; b++){
				double g = 0;

				for(unsigned i = 0; i < ns; i++)
					g += x[i*nn2 + a] * x[shift + i*nn2 + b];

				const unsigned pos = tau*nn2*nn2 + a*nn2 + b;

				greens[pos] = g/ns;
			}
		}
	}
}

double logDetM_and_inv_M_hat(double *M_hat, double *w, int *p, unsigned ns, unsigned flavors, organic_flags *mode){
	double log_det = 0;

	if(mode->single_particle){
		log_det = log(trace(M_hat, ns)); 
	}else{
		//print_mat(M_hat, ns);
		const unsigned mat_dim = ns*ns;

		if(mode->schur_form_inv){
			for(unsigned i = 0; i < ns; i++) M_hat[i*(ns+1)] -= 1;
			int sdim;
			double *wr = w + mat_dim, *wi = wr + ns;
			LAPACKE_dgees(LAPACK_COL_MAJOR, 'V', 'N', NULL, ns, M_hat, ns, &sdim, wr, wi, w, ns);

			double min_ev = wr[0], max_ev = wr[0];
			for(unsigned i = 1; i < ns; i++){
				if(wr[i] < min_ev) min_ev = wr[i];
				if(wr[i] > max_ev) max_ev = wr[i];
			}

			const double eps = max_ev*DBL_EPSILON/1e-3;
			if(min_ev < -1 || eps > 1){
				for(unsigned i = 0; i < mat_dim; i++){
					if(fabs(M_hat[i]) < eps) M_hat[i] = 0;
				}
			}

			for(unsigned i = 0; i < ns; i++) M_hat[i*(ns+1)] += 1;

			//print_mat(M_hat, ns);
		}

		LAPACKE_dgetrf(LAPACK_COL_MAJOR, ns, ns, M_hat, ns, p);

		int info = 1;
		for(unsigned i = 0; i < ns; i++){
			const double diag = M_hat[i*(ns+1)];
			log_det += log(fabs(diag));
			info *= diag > 0? 1: -1;
			info *= p[i] == i+1? 1: -1;
		}
		if(info < 0 && flavors % 2) printf("Negative determinant!\n");

		//for(unsigned i = 0; i < ns; i++) printf("%d, ", p[i]);
		//for(unsigned i = 0; i < ns; i++) printf("%g, ", M_hat[i*(ns+1)]);
		//printf("\n");

		LAPACKE_dgetri(LAPACK_COL_MAJOR, ns, M_hat, ns, p);

		//print_mat(M_hat, ns);

		if(mode->schur_form_inv){
			double *mat_prod = w + 2*mat_dim;
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, ns, ns, ns, 1, M_hat, ns, w, ns, 0, mat_prod, ns);
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ns, ns, ns, 1, w, ns, mat_prod, ns, 0, M_hat, ns);
		}
	}

	return log_det;
}

#ifdef EXACT_INVERSE_MAIN
int main(){
	//const unsigned ns = 4;
	double det;
	const unsigned l = 3;
	unsigned ns, nn;

	double t[3] = {1,2,3};
	double kappa[3] = {1,1,1};

	unsigned *nnt = construct_triangular(l, l, &ns, &nn);

	double *x = malloc(ns*nn/2 * sizeof(double));
	double *M_hat = malloc(ns*ns * sizeof(double));
	double *w = malloc(ns * sizeof(double));
	int *p = malloc(ns * sizeof(int));

	for(unsigned i = 0; i < ns; i++){
		for(unsigned k = 0; k < nn/2; k++) x[i*nn/2 + k] = 1.*i/(k+1);
	}
	//print_vec(x, ns*nn/2);

	construct_M_hat(x, M_hat, w, t, kappa, 0, nnt, ns, nn, 1, mode);
	print_mat(M_hat, ns);

	//double mat[16] = {4.022820, -3.146379, 1.028209, 1.917387, -3.146379, 3.396020, -3.827133, -1.976714, 3.028209, -3.827133, 5.015514, 1.980140, 1.917387, -1.976714, 1.980140, 1.730459};

	det = exp(logDetM_and_inv_M_hat(M_hat, w, p, ns, mode));
	printf("det = %.8g\n", det);

	free(p);
	free(w);
	free(M_hat);
	free(nnt);

	return 0;
}
#endif
