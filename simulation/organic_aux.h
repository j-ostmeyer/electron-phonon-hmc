#ifndef ORGANIC_AUX
#define ORGANIC_AUX

#define NUM_RES 12

#ifndef M_PI
#define M_PI 3.141592653589793
#endif
#define M_2PI 6.283185307179586
#define HALF_M_PI 1.570796326794897

#define max(a,b)              \
	({                        \
	 __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b;       \
	 })

#define min(a,b)              \
	({                        \
	 __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a < _b ? _a : _b;       \
	 })

unsigned ceil_half(unsigned n);
unsigned ceil_log2(unsigned n);

double sinhc(double x);
double tanhc(double x);

double norm(double complex x);

void print_vec(double *vec, unsigned n);
void print_mat(double *m, unsigned long n);

void fprint_corr(FILE *out, double *corr, double *contr_mat, unsigned nn, unsigned nt);
void fprint_results(FILE *out, double *res, unsigned n, unsigned r);
void fprint_sigma(const char *name, double *greens, double width, unsigned n);
void fwrite_sigma(const char *name, double *greens, double width, unsigned n);

void print_config(double *x, unsigned ns, unsigned nn, unsigned nt);
void fprint_config(const char *name, double *x, unsigned ns, unsigned nn, unsigned nt, const char *mode);
void fwrite_config(const char *name, double *x, unsigned n, const char* mode);

double *fscan_array(const char *in_name, unsigned *length, double *beta, unsigned *nt);
double *fscan_greens(const char *in_name, double *h, unsigned *length);
double *fread_greens(const char *in_name, double *h, unsigned *length);

double average(double *x, unsigned n);
double average_sq(double *x, unsigned n);
void average_links(double *x, double *scale, unsigned ns, unsigned nn, unsigned nt, double *res);

double scalar_dot(double *x, double *y, unsigned d);

#endif
