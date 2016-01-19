#ifndef INIT
#define INIT

#include <vector>
#include <complex>
#include <random>

struct para{
    double mu;
    double gam;
    double sigma_j;
    double rho;
    double k;
    double v_p;
    double sigma_v;
    double weight;
    double norm_weight;
    double cum_norm_weight;
    double* v;
    double* v_star;
    double* z;
    double* w;
    double u1;
    double u2;
    double u3;
    para* next;
};

extern double max_v;
extern double max_w;
extern double min_v;
extern double min_w;
extern para* head;

extern std::vector<double> price;
extern double delta;
extern double dmax;
extern int total;
extern int no_particles;
extern int flag;
extern double ep;
extern double TOL;

//init.cpp
void init();
double normal();
double uniform();
void print(para* p);
//
//posterior.cpp
double prior(para* p);
double posterior(double zeta,para* p);
double post_no_prior(double zeta,double z,double vt,double vu,double v_star,double yt,double yu,double w,para* p);
double likelihood(double zeta,double z,double vt,double vu,double v_star,double yt,double yu,para* p);
double variance_gamma(double z,para* p);
double transition(double vt,double vu,para* p);
double aux_g(double v,double w,double vt,double vu,para* p);
std::complex<double> phi(double v,double w,double vt,double vu,para* p);
//
//resample_update.cpp
double ESS(double zeta,double prev_zeta);
double find_new_zeta(double prev_zeta,double u_zeta,double l_zeta,double curr_ESS);
void update_norm_weights();
void resample();
void copy_particle(para* new_curr,para* old,int i);
void destroy(para* p);
void resample_tester();
//
//kernel.cpp
para* set_kernel();
//
//update_para.cpp
void update_para(double zeta,para* kernel);
void update_mu(double zeta,para* curr,double sd);
void update_gam(double zeta,para* curr,double sd);
void update_sj(double zeta,para* curr,double sd);
void update_rho(double zeta,para* curr,double sd);
void update_k(double zeta,para* curr,double sd);
void update_vp(double zeta,para* curr,double sd);
void update_sv(double zeta,para* curr,double sd);
//
//bessel.cpp
double chebev(double a, double b, double c[], int m, double x);
void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi);
double bessik(double x, double xnu);
//
//rnglib
void advance_state ( int k );
bool antithetic_get ( );
void antithetic_memory ( int i, bool &value );
void antithetic_set ( bool value );
void cg_get ( int g, int &cg1, int &cg2 );
void cg_memory ( int i, int g, int &cg1, int &cg2 );
void cg_set ( int g, int cg1, int cg2 );
int cgn_get ( );
void cgn_memory ( int i, int &g );
void cgn_set ( int g );
void get_state ( int &cg1, int &cg2 );
int i4_uni ( );
void ig_get ( int g, int &ig1, int &ig2 );
void ig_memory ( int i, int g, int &ig1, int &ig2 );
void ig_set ( int g, int ig1, int ig2 );
void init_generator ( int t );
void initialize ( );
bool initialized_get ( );
void initialized_memory ( int i, bool &initialized );
void initialized_set ( );
void lg_get ( int g, int &lg1, int &lg2 );
void lg_memory ( int i, int g, int &lg1, int &lg2 );
void lg_set ( int g, int lg1, int lg2 );
int multmod ( int a, int s, int m );
float r4_uni_01 ( );
double r8_uni_01 ( );
void set_initial_seed ( int ig1, int ig2 );
void set_seed (  int cg1, int cg2 );
void timestamp ( );
// ranlib
char ch_cap ( char ch );
float genbet ( float aa, float bb );
float genchi ( float df );
float genexp ( float av );
float genf ( float dfn, float dfd );
float gengam ( float a, float r );
float *genmn ( float parm[] );
int *genmul ( int n, float p[], int ncat );
float gennch ( float df, float xnonc );
float gennf ( float dfn, float dfd, float xnonc );
float gennor ( float av, float sd );
void genprm ( int iarray[], int n );
float genunf ( float low, float high );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int ignbin ( int n, float pp );
int ignnbn ( int n, float p );
int ignpoi ( float mu );
int ignuin ( int low, int high );
int lennob ( char *s );
void phrtsd ( char *phrase, int &seed1, int &seed2 );
void prcomp ( int maxobs, int p, float mean[], float xcovar[], float answer[] );
float r4_exp ( float x );
float r4_exponential_sample ( float lambda );
float r4_max ( float x, float y );
float r4_min ( float x, float y );
float r4vec_covar ( int n, float x[], float y[] );
double r8_exponential_sample ( double lambda );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8vec_covar ( int n, double x[], double y[] );
int s_eqi ( char *s1, char *s2 );
float sdot ( int n, float dx[], int incx, float dy[], int incy );
float *setcov ( int p, float var[], float corr );
void setgmn ( float meanv[], float covm[], int p, float parm[] );
float sexpo ( );
float sgamma ( float a );
float snorm ( );
int spofa ( float a[], int lda, int n );
void stats ( float x[], int n, float &av, float &var, float &xmin, float &xmax );
void trstat ( std::string pdf, float parin[], float &av, float &var );
//

#endif