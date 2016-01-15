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
    double* v;
    double* v_star;
    double* z;
    double* w;
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
void initialize();
double normal();
double uniform();
void print();
double prior(para* p);
double posterior(double zeta,para* p);
double post_no_prior(double zeta,double z,double vt,double vu,double v_star,double yt,double yu,double w,para* p);
double likelihood(double zeta,double z,double vt,double vu,double v_star,double yt,double yu,para* p);
double variance_gamma(double z,para* p);
double transition(double vt,double vu,para* p);

/*
double post_no_prior(double z,double vt,double vu,double v_star,double yt,double yu,double w);
double posterior();
double prior();
double likelihood(double z,double vt,double vu,double v_star,double yt,double yu);
double variance_gamma(double z);
double aux_g(double v,double w,double vt,double vu);
std::complex<double> phi(double v,double w,double vt,double vu);
double transition(double vt,double vu);

double sd(std::vector<double> history);
double mean(std::vector<double> history);
double normal(double s1,double s2);
double uniform(unsigned s);
double accept(std::vector<double> hist);

double update_gam();
double update_mu();
double update_sigma_j();
double update_rho();
double update_k();
double update_v_p();
double update_sigma_v();

void update_latent_w(int iter,double** hist);
void update_latent_z(int iter,double** hist);
void update_latent_v_s(int iter,double** hist);
void update_latent_v(int iter,double** hist);

double mean(int iter,double* history);
double sd(int iter,double* history);
void accept_latent(double** hist);
double scaler(std::vector<double> history);

double final_mean(std::vector<double> history,int cutoff);
double accept_latent(int iter,double* hist);
double latent_mean(double ** hist);*/
#endif