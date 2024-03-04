#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "organic_aux.h"

#include "organic_lattice.h"

unsigned *construct_triangular(unsigned l1, unsigned l2, unsigned *ns, unsigned *nn){
	// constructs a parallelogram with periodic boundaries
	unsigned *nnt;
	int i1, i2, shift;

	*nn = 6;
	*ns = l1*l2;

	const unsigned size = *ns;
	const unsigned neighbours = *nn, links = neighbours/2;
	nnt = malloc(2 * neighbours * size * sizeof(unsigned));

	for(unsigned i = 0; i < size; i++){
		i1 = i%l1;
		i2 = i/l1;

		// nearest sites
		shift = i*neighbours;
		nnt[shift] = i + (i1+1)%l1 - i1; // right
		nnt[shift+1] = i + (i1+l1-1)%l1 - i1; // left
		nnt[shift+2] = i + ((i2+1)%l2 - i2)*l1; // top right
		nnt[shift+3] = i + ((i2+l2-1)%l2 - i2)*l1; // bottom left
		nnt[shift+4] = nnt[shift+2] + (i1+l1-1)%l1 - i1; // top left
		nnt[shift+5] = nnt[shift+3] + (i1+1)%l1 - i1; // bottom right

		// nearest links
		shift = (size + i)*neighbours;
		nnt[shift] = i*links; // right
		nnt[shift+1] = nnt[i*neighbours+1]*links; // left
		nnt[shift+2] = i*links + 1; // top right
		nnt[shift+3] = nnt[i*neighbours+3]*links + 1; // bottom left
		nnt[shift+4] = i*links + 2; // top left
		nnt[shift+5] = nnt[i*neighbours+5]*links + 2; // bottom right
	}

	return nnt;
}

double *construct_contraction_mat(double *a, double *b, unsigned n){
	// project lattice vectors onto Cartesian basis
	const unsigned d = 2;
	double *x = malloc(d*n * sizeof(double));

	for(unsigned i = 0; i < d; i++){ // here we assume triangular structure and n=3
		x[i*n]   = a[i]; //right
		x[i*n+1] = .5 * ( a[i] + b[i]); // top right
		x[i*n+2] = .5 * (-a[i] + b[i]); // top left
	}
	
	return x;
}

double contr_gr_mat(double *gr, double *mat, unsigned n){
	// calculating the trace in Cartesian coordinates
	const unsigned d = 2;
	double contr = 0;

	for(unsigned c = 0; c < d; c++){ // Cartesian trace
		double sum = 0;
		for(unsigned a = 0; a < n; a++){ // lattice row projection
			//for(unsigned b = 0; b < n; b++){ // lattice column projection
			//	contr += mat[c*n + a] * gr[a*n + b] * mat[c*n + b];
			//}
			sum += mat[c*n + a] * gr[a];
		}
		contr += sum*sum;
	}

	return contr;
}

double contr_gr_triangle(double *gr, unsigned n){
	double contr = 0;

	// summing symmetrised version with diag 1 and all off-diag entries 1/2
	for(unsigned a = 0; a < n; a++){
		for(unsigned b = 0; b <= a; b++) contr += gr[a*n + b] + gr[b*n + a];
	}
	contr *= .5;

	// subtracting corner terms
	contr -= gr[n-1] + gr[(n-1)*n];

	return contr;
}

double contract_greens(double *gr, double *mat, unsigned nn, unsigned tau){
	const unsigned nn2 = nn/2;
	double contraction;

	if(mat)
		contraction = contr_gr_mat(gr + nn2*nn2*tau, mat, nn2);
	else
		contraction = contr_gr_triangle(gr + nn2*nn2*tau, nn2);

	return contraction;
}
