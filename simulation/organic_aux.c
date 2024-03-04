#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <omp.h>

#include "organic_lattice.h"
#include "organic_aux.h"

unsigned ceil_half(unsigned n){
	return n%2? n/2 + 1: n/2;
}

unsigned ceil_log2(unsigned n){
	unsigned l = 0;
	for(l = 0; n > 1; l++) n = ceil_half(n);
	return l;
}

double sinhc(double x){
	return x? sinh(x)/x: 1;
}

double tanhc(double x){
	return x? tanh(x)/x: 1;
}

double norm(double complex x){
	const double a = creal(x), b = cimag(x);
	return a*a + b*b;
}

void print_vec(double *vec, unsigned n){
	unsigned i;
	for(i = 0; i < n; i++) printf("%g\n", vec[i]);
}

void print_mat(double *m, unsigned long n){
	unsigned i, k;

	printf("{");
	for(i = 0; i < n; i++){
		printf("{");
		for(k = 0; k < n; k++, m++) printf("%.15g,\t", *m);
		printf("},\n");
	}
	printf("}\n");
}

void fprint_corr(FILE *out, double *corr, double *contr_mat, unsigned nn, unsigned nt){
	if(out){
		const unsigned nt2 = nt*nt;

		for(unsigned i = 0; i < nt; i++)
			fprintf(out, "%d\t%.15g\n", i, contract_greens(corr, contr_mat, nn, i) * nt2);
		//fprintf(out, "\n");
		fflush(out);
	}
}

void fprint_results(FILE *out, double *res, unsigned n, unsigned r){
	if(out){
		for(unsigned i = 0; i < n; i++){
			for(unsigned k = 0; k < r; k++, res++) fprintf(out, "%.15g\t", *res);
			fprintf(out, "\n");
		}
		fflush(out);
	}
}

void fprint_sigma(const char *name, double *greens, double width, unsigned n){
	FILE *out = fopen(name, "a");

	for(unsigned i = 0; i < n; i++){
		fprintf(out, "%.15g\t%.15g\n", i*width, greens[i]);
	}

	fclose(out);
}

void fwrite_sigma(const char *name, double *greens, double width, unsigned n){
	FILE *out = fopen(name, "w");

	fwrite(&n, sizeof(unsigned), 1, out);
	fwrite(&width, sizeof(double), 1, out);
	fwrite(greens, sizeof(double), n, out);

	fclose(out);
}

void print_config(double *x, unsigned ns, unsigned nn, unsigned nt){
	for(unsigned tau = 0; tau < nt; tau++){
		for(unsigned i = 0; i < ns; i++){
			for(unsigned k = 0; k < nn/2; k++, x++) printf("%g, ", *x);
			printf("\t");
		}
		printf("\n");
	}
}

void fprint_config(const char *name, double *x, unsigned ns, unsigned nn, unsigned nt, const char *mode){
	const unsigned nn2 = nn/2;
	FILE *out = fopen(name, mode);
	for(unsigned tau = 0; tau < nt; tau++){
		for(unsigned i = 0; i < ns; i++){
			for(unsigned k = 0; k < nn2; k++, x++) fprintf(out, "%.16g\t", *x);
			fprintf(out, "\n");
		}
	}
	fprintf(out, "\n");
	fclose(out);
}

void fwrite_config(const char *name, double *x, unsigned n, const char* mode){
	FILE *out = fopen(name, mode);
	fwrite(x, sizeof(double), n, out);
	fclose(out);
}

double *fscan_array(const char *in_name, unsigned *length, double *beta, unsigned *nt){
	unsigned i, n;
	FILE *in = fopen(in_name, "r");

	if(beta) fscanf(in, "%lg ", beta);
	if(nt) fscanf(in, "%u ", nt);

	n = 1024;
	double *x = malloc(n * sizeof(double));
	for(i = 0; fscanf(in, "%lg ", x+i) > 0;){
		if(++i == n){
			n *= 2;
			x = realloc(x, n*sizeof(double));
		}
	}
	*length = i;
	
	fclose(in);

	return x;
}

double *fscan_greens(const char *in_name, double *h, unsigned *length){
	unsigned i, n;
	double dummy;
	FILE *in = fopen(in_name, "r");

	n = 1024;
	double *x = malloc(n * sizeof(double));
	fscanf(in, "%lg %lg ", h, x);
	fscanf(in, "%lg %lg ", h, x+1);
	for(i = 2; fscanf(in, "%lg %lg ", &dummy, x+i) > 0;){
		if(++i == n){
			n *= 2;
			x = realloc(x, n*sizeof(double));
		}
	}
	*length = i;
	
	fclose(in);

	return x;
}

double *fread_greens(const char *in_name, double *h, unsigned *length){
	FILE *in = fopen(in_name, "r");

	fread(length, sizeof(unsigned), 1, in);
	fread(h, sizeof(double), 1, in);

	const unsigned n = *length;

	double *x = malloc(n * sizeof(double));
	fread(x, sizeof(double), n, in);
	
	fclose(in);

	return x;
}

double average(double *x, unsigned n){
	double res = 0;

#pragma omp parallel for reduction(+:res)
	for(unsigned i = 0; i < n; i++) res += x[i];

	return res/n;
}

double average_sq(double *x, unsigned n){
	double res = 0;

#pragma omp parallel for reduction(+:res)
	for(unsigned i = 0; i < n; i++) res += pow(x[i], 2);

	return res/n;
}

void average_links(double *x, double *scale, unsigned ns, unsigned nn, unsigned nt, double *res){
	const unsigned n = ns*nt;

	for(unsigned k = 0; k < nn; k++) res[k] = 0;

	for(unsigned i = 0; i < n; i++){
		for(unsigned k = 0; k < nn; k++) res[k] += x[i*nn+k];
	}

	//for(unsigned k = 0; k < nn; k++) res[k] *= scale[k]/ns;
	for(unsigned k = 0; k < nn; k++) res[k] /= ns;
}

double scalar_dot(double *x, double *y, unsigned d){
	double total = 0;

	for(unsigned i = 0; i < d; i++) total += x[i]*y[i];

	return total;
}
