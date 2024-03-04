#ifndef ORGANIC_FLAGS
#define ORGANIC_FLAGS

typedef struct OrganicFlags{
	int schur_form_inv;
	int zero_filling;
	int static_disorder;
	int single_particle;
	int constant_force;
	int no_fourier_acc;
	int classical_fourier_acc;
	unsigned cfa_steps;
	double cfa_mass;
} organic_flags;

void set_flags(FILE *in, organic_flags *mode);

#endif
