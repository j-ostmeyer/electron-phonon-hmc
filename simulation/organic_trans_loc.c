#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <float.h>
#include <time.h>

#ifndef MKL_LAPACKE
#include <lapacke.h>
#else
#include <mkl.h>
#include <mkl_lapacke.h>
#endif

#include "mt19937-64.h"
#include "organic_random.h"

#include "organic_flags.h"
#include "organic_aux.h"
#include "organic_lattice.h"
#include "organic_dense_matrix.h"

#include "organic_trans_loc.h"

void hist_trans_loc(double *x, double *t, double *kappa, double mu, double *contr_mat, unsigned *nnt, unsigned ns, unsigned nn, unsigned nt, double width, double *results, double *greens, FILE *res_out){
	const unsigned nn2 = nn/2, dim = ns*nn2, mat_dim = ns*ns;

	double partition_sum = 0;

	double *ham = x + dim; // mat_dim, Hamiltonian matrix
	double *w = x + dim + mat_dim; // mat_dim, buffer
	double *ev = x + dim + 2*mat_dim; // ns, eigen-energies
	double *correl = x + dim + 2*mat_dim + ns; // nn2*nn2, current-current-correlator elements
	double *proj = x + dim + 2*mat_dim + ns + nn2*nn2; // nn2, projection <m|I_k|n> for each k

	//double *all_I_nm = calloc(mat_dim, sizeof(double));

	random_vector(x, dim);
	construct_hamilton(x, ham, w, t, kappa, -mu, nnt, ns, nn);
	LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', ns, ham, ns, ev);

	//printf("mu = %g, E_min = %g, E_max = %g\n", mu, ev[0], ev[ns-1]);

	//fwrite_config("energies.bin", ev, ns, "a");

	for(unsigned i = 0; i < ns; i++){
		for(unsigned k = 0; k < nn2; k++)
			apply_current_to_v(x, ham + i*ns, w + k*ns, t, kappa, nnt, ns, nn, k);

		for(unsigned j = 0; j < ns; j++){
			if(ev[j] < ev[i]){ // only collect positive frequencies
				//if(j < i) all_I_nm[i*ns + j] = all_I_nm[j*ns + i];
				//else printf("ERROR: E[%d] = %g < %g = E[%d]\n", j, ev[j], ev[i], i);
				continue;
			}

			//const double boltzmann = exp(-.5*(ev[i] + ev[j]));
			const double boltzmann = exp(-ev[i]) + exp(-ev[j]);
			const double freq = ev[j] - ev[i];
			const unsigned bin = (unsigned) (freq / width + .5);

			for(unsigned k = 0; k < nn2; k++)
				proj[k] = scalar_dot(ham + j*ns, w + k*ns, ns);

			for(unsigned a = 0; a < nn2; a++){
				for(unsigned b = 0; b < nn2; b++) correl[a*nn2 + b] = proj[a]*proj[b];
			}

			const double I_nm = contr_gr_mat(proj, contr_mat, nn2);
			//all_I_nm[i*ns + j] = I_nm;

			greens[bin] += boltzmann * I_nm / (width * ns);
		}

		partition_sum += exp(-ev[i]);
	}

	//fwrite_config("I_nm_elements.bin", all_I_nm, mat_dim, "a");
	//free(all_I_nm);

	results[0] = partition_sum / ns;
	fprint_results(res_out, results, 1, 1);
}

unsigned run_trans_loc(double *greens, unsigned l1, unsigned l2, unsigned nt, double *t, double *kappa, double mu, double *contr_mat, unsigned meas, const char *res_name, const char *corr_name){
	unsigned ns, nn;

	unsigned *nnt = construct_triangular(l1, l2, &ns, &nn);

	const unsigned nn2 = nn/2, dim = ns*nn2, mat_dim = ns*ns;

	double *results = greens + nt;
	double *x = malloc((dim + 2*mat_dim + ns + (nn2+1)*nn2) * sizeof(double));

	FILE *res_out = res_name? fopen(res_name, "a"): NULL;

	init_genrand64(time(NULL));

	double width = -2*mu/nt;
	for(unsigned i = 0; i < nn2; i++) width += 8*fabs(kappa[i])/nt;

	for(unsigned i = 0; i < meas; i++)
		hist_trans_loc(x, t, kappa, mu, contr_mat, nnt, ns, nn, nt, width, results + i, greens, res_out);

	for(unsigned k = 0; k < nt; k++) greens[k] /= meas;

	fwrite_sigma(corr_name, greens, width, nt);

	free(nnt);
	free(x);
	if(res_out) fclose(res_out);

	return meas + nt; // dimension of results array for sanity check
}
