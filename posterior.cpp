#include "include.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <complex>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <complex_bessel.h>
#define PII 3.141592654

using namespace std;

void update_all(){
    para* p = head;
    double sum,val;
    int count = 1;
    while(p != NULL){
        cout<<count<<endl;
        sum = 0;
        for(int i=1;i<total;i++){
            val = variance_gamma(p->z[i],p);
            if(abs(val-p->post_z[i])>TOL){
                cout<<val<<" wrong z "<<p->post_z[i]<<endl;
                p->post_z[i] = val;
            }
            val = likelihood(1.0,p->z[i],p->v[i],p->v[i-1],price[i],price[i-1],p);
            if(abs(val- p->post_lik[i])>TOL){
                cout<<val<<" wrong lik "<<p->post_lik[i]<<endl;
                p->post_lik[i] = val;
            }
            val = p->post_z[i]+p->post_lik[i];
            if(abs(val - p->posterior[i])>TOL){
                cout<<"id:"<<i<<" corr:"<<val<<" wrong posterior "<<"wrong:"<<p->posterior[i]<<endl;
                p->posterior[i] = val;
            }
            sum += p->posterior[i];
        }
        if(abs(sum+prior(p)- p->sum_post)>TOL){
            cout<<sum<<" wrong sum "<<p->sum_post<<endl;
            p->sum_post = sum;
        }
        count++;
        p = p->next;
    }
    return ;
}

void print_part_post(){
    para* p = head;
    while(p != NULL){
        for(int i=1;i<total;i++){
            p->post_z[i] = variance_gamma(p->z[i],p);
            p->post_lik[i] = likelihood(0.005,p->z[i],p->v[i],p->v[i-1],price[i],price[i-1],p);
            if(p->post_z[i]>0||p->post_lik[i]>0)
                cout<<"z: "<<p->post_z[i]<<" lik: "<<p->post_lik[i]<<endl;
        }
        p = p->next;
    }
}
void update_untempered_lik(){
    int j =0;
    para* p = head;
    while(p != NULL){
        untempered_lik[j] = lik(1.0,p);
        p = p->next;
        j++;
    }
    return ;
    
}
double prior(para* p){
    double prior_mu,prior_gam,prior_sigj,prior_rho,prior_sigv,prior_k,prior_vp,w_v,phi_v;
    //int I = 1;
    prior_mu = -pow(p->mu,2)/2-log(sqrt(2*PII));
    prior_gam = -pow(p->gam,2)/2-log(sqrt(2*PII));
    phi_v = p->rho*p->sigma_v;
    w_v = p->sigma_v*p->sigma_v*(1-p->rho*p->rho);
    prior_rho = log(exp(-0.5/w_v)*pow(w_v,-1.0-1)*pow(0.5,1)/tgamma(1));
    //double rhob = (p->rho+1)/2;
    //prior_rho = (1/7-1)*log(rhob)+log(tgamma(1.0+1/7)/tgamma(0.14285714));
    prior_sigj = log(exp(-5/p->sigma_j)*pow(p->sigma_j,-2.5-1)*pow(5,2.5)*tgamma(2.5));
    prior_sigv = -pow(phi_v,2)/(2*pow(0.5/w_v,2))-log(sqrt(2*PII*pow(0.5/w_v,2)));
    prior_k = -pow(p->k,2)/2-log(sqrt(2*PII));
    prior_vp = -pow(p->v_p,2)/2-log(sqrt(2*PII));
    //double d = 4.0*p->k*p->v_p/pow(p->sigma_v,2);
    //cout<<w_v<<" "<<phi_v<<" sigv:"<<p->sigma_v<<" rho:"<<p->rho<<" prior_sigv:"<<prior_sigv<<endl;
    //cout<<prior_mu+prior_gam+prior_rho+prior_sigv+prior_sigj+prior_k+prior_vp<<endl;
    //cout<<prior_mu<<" "<<prior_gam<<" "<<prior_rho<<" "<<prior_sigv<<" "<<prior_sigj<<" "<<prior_k<<" "<<prior_vp<<endl;
    return prior_mu+prior_gam+prior_rho+prior_sigv+prior_sigj+prior_k+prior_vp;
}
double mcmc_posterior(para* p,double* lat_z,double* lat_lik,double* lat_posterior){   //untempered posterior
    double val = 0;
    lat_lik[0] = 0;
    lat_z[0] = 0;
    lat_posterior[0]=0;
    for(int i =1;i<total;i++){
        lat_lik[i] = likelihood(1.0,p->z[i],p->v[i],p->v[i-1],price[i],price[i-1],p);
        lat_z[i] = variance_gamma(p->z[i],p);
        lat_posterior[i] = lat_z[i]+lat_lik[i];
        val = val+ lat_posterior[i];
    }
    //cout<<"prior1:"<<prior(p)<<endl;
    //cout<<"mcmc_fn: ";
    return val+prior(p);
}
double posterior(double zeta,para* p){
    //cout<<"posterior_fn: ";
    double val1 = prior(p);
    double sum= 0 ;
    for(int i =1;i<total;i++){
        sum = sum+post_no_prior(zeta,p->z[i],p->v[i],p->v[i-1],price[i],price[i-1],p);
    }
    //cout<<sum<<" "<<val1<<endl;
    //cout<<"prior: "<<val1<<" transition: "<<sum<<endl;
    //cout<<"prior2m:"<<val1<<endl;
    return val1+sum;
}


double post_no_prior(double zeta,double z,double vt,double vu,double yt,double yu,para* p){
    double val2 = likelihood(zeta,z,vt,vu,yt,yu,p);
    double val3 = variance_gamma(z,p);
    //if(val2>0) cout<<"like>1";
    //if(val3>0) cout<<"vg>1";
    //cout<<endl;
    return val2+val3;
}

double lik(double zeta,para* p){
    double sum = 0;
    for(int i=1;i<total;i++){
        sum += likelihood(zeta,p->z[i],p->v[i],p->v[i-1],price[i],price[i-1],p);
        //cout<<likelihood(zeta,p->z[i],p->v[i],p->v[i-1],price[i],price[i-1],p)<<endl;
    }
    return sum*zeta;
}
double likelihood(double zeta,double z,double vt,double vu,double yt,double yu,para* p){
    double** cov;double** inv;
    cov = new double*[2];
    inv = new double*[2];
    for(int i =0;i<2;i++){
        cov[i] = new double[2];
        inv[i] = new double[2];
    }
    cov[0][0] = vu*delta;
    cov[0][1] = vu*delta*p->rho*p->sigma_v;
    cov[1][0] = cov[0][1];
    cov[1][1] = cov[0][0]*pow(p->sigma_v,2);
    inverse(inv,cov);

    double var = (1-pow(p->rho,2))*vu*delta;
    double norm_v = (vt-(vu-p->k*(p->v_p-vu)*delta))/sqrt(vu*delta);
    double norm_y = (yt-(yu + p->mu*delta + z))/(p->sigma_v*sqrt(vu*delta));
    double val = -0.5*(norm_y*norm_y+norm_v*norm_v-2*p->rho*norm_v*norm_y)/(1-p->rho*p->rho)+log(1/(sqrt(2*PII*(1-p->rho*p->rho))*p->sigma_v*vu));
    //cout<<-0.5*(norm_y*norm_y+norm_v*norm_v-2*p->rho*norm_v*norm_y)/(1-p->rho*p->rho)<<" "<<log(1/(sqrt(2*PII*(1-p->rho*p->rho))*p->sigma_v*vu))<<endl;
    //double val = -0.5*mult(yt-mn_y,vt-mn_v,inv)+log(1/sqrt(determinant(cov)*2*PII));cout<<determinant(cov)<<endl;
    //cout<<-log(sqrt(determinant(cov)*2*PII))<<" "<<-0.5*mult(yt-mn_y,vt-mn_v,inv)<<" "<<val<<endl;
    for(int i=0;i<2;i++){
        delete cov[i];
        delete inv[i];
    }
    delete cov;delete inv;
    return val*zeta;
}



double variance_gamma(double z,para* p){

    double alpha = delta*1.0;
    double sigma = p->sigma_j*sqrt(delta);
    double val = -log(sigma);
    sigma  = sigma*sigma;
    val += p->gam*z/sigma;
    val -= delta/alpha*log(alpha);
    //cout<<"1:"<<val<<endl;
    val -= tgamma(1.0*delta/alpha);
    //cout<<"2:"<<val<<endl;
    val += (0.5*delta/alpha-0.25)*log(z*z/(p->gam*p->gam+2*sigma/alpha));
    
    double va = (sqrt(p->gam*p->gam+2*sigma/alpha)*abs(z))/sigma;
    //cout<<va<<" "<<val<<endl;
    if(abs(alpha-delta)<TOL)
    {
        val += -0.5*log(va) - va;
    }
    else
    {
        val += log(bessik(va,1.0*delta/alpha - 0.5));
    }
    /*
    double val1 = log(2)+p->gam*pow(z,2)/sigma-log(alpha)*delta/alpha-log(tgamma(1.0*delta/alpha))-log(sqrt(2*PII));
    double val2 = (0.5*delta/alpha-0.25)*log(z*z/(p->gam*p->gam + 2.0*sigma*sigma/alpha));
    double va = sqrt(pow(z,2)*(pow(p->gam,2)+2*pow(sigma,2)/alpha))/pow(sigma,2);
    double val3 = bessik(va,1.0*delta/alpha-0.5);*/
    //cout.clear();
    return val;//val1+log(val3)+val2;
}

void swap(double* val1, double* val2){
    double temp = *val1;
    *val1 = *val2;
    *val2 = temp;
    return ;
}

void print_mat(double** mat,int r,int c){
    cout<<"start"<<endl;
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            cout<<mat[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<"end"<<endl;
    return ;
}

double determinant(double** mat){
    return mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1];
}

void inverse(double** inv,double** mat){
    double det = determinant(mat);
    inv[0][0] = mat[1][1]/det;
    inv[1][1] = mat[0][0]/det;
    inv[1][0] = -mat[1][0]/det;
    inv[0][1] = -mat[0][1]/det;
    return ;
}

void multiply(double** res,double** mat1,double** mat2,int r1,int c1,int r2,int c2){
    double sum;
    for(int i=0;i<r1;i++){
        for(int j=0;j<c2;j++){
            sum = 0;
            for(int q=0;q<c1;q++)
                sum += mat1[i][q]*mat2[q][j];
            res[i][j] = sum;
        }
    }
    return ;
}
double mult(double norm_y,double norm_v,double** mat){
    double val[2];
    val[0] = mat[0][0]*norm_y+mat[1][0]*norm_v;
    val[1] = mat[0][1]*norm_y+mat[1][1]*norm_v;
    double ans = norm_y*val[0]+norm_v*val[1];
    //cout<<norm_y<<" "<<norm_v<<" "<<val[0]<<" "<<val[1]<<" "<<ans<<endl;
    return ans;
}
/*
double aux_g(double v,double w,double vt,double vu,para* p){
    int I = 0;
    if(v>min_v)
        if(v<max_v)
            if(w>min_w)
                if(w<max_w)
                    I=1;
    std::complex<double> val = phi(v,w,vt,vu,p);
    double g = cos(v*w)*val.real()+sin(v*w)*val.imag();
    double aux = abs(g)*I+abs(g)*exp(-v-w)*(1-I);
    //cout<<aux<<endl;
    return log(aux);
}

std::complex<double> phi(double v,double w,double vt,double vu,para* p){
    double de = 1.0*delta;
    double com = -2*pow(p->sigma_v,2)*w;
    double lambda = 4.0*p->k/(pow(p->sigma_v,2)*(1-exp(-p->k*delta)));
    complex<double> nu(pow(p->k,2),com);
    nu = sqrt(nu);
    double d = 0.5*4*p->k*p->v_p/pow(p->sigma_v,2)-1;
    complex<double> v1 = sqrt(vt*vu)*4.0*nu*exp(-0.5*nu*de)/(pow(p->sigma_v,2)*(1.0-exp(-nu*de)));
    complex<double> v2 = sqrt(vt/lambda*vu/lambda)*exp(-0.5*p->k*de)/(1.0-exp(-p->k*de));
    complex<double> den = besselI(d,v2);
    //cout<<"prob"<<endl;
    complex<double> val1 = besselI(d,v1);
    //cout<<" v2: "<<v2<<" v1: "<<v1<<" den: "<<den<<" val1: "<<val1<<" divide:"<<val1/den<<endl;
    
    val1 = val1/den*nu*exp(-0.5*(nu-p->k)*de)*(1.0-exp(-p->k*de))/(p->k*(1.0-exp(-nu*de)));
    val1 = val1*exp((vu+vt)/pow(p->sigma_v,2)*(p->k*(1.0+exp(-p->k*de))/(1.0-exp(-p->k*de))-nu*(1.0+exp(-nu*de))/(1.0-exp(-nu*de))));

    return val1;
}
 
 double transition(double vt,double vu,para* p){
 double lambda = 4.0*p->k/(pow(p->sigma_v,2)*(1-exp(-p->k*delta)));
 double  d = 4*p->k*p->v_p/pow(p->sigma_v,2);
 
 double val;
 boost::math::non_central_chi_squared_distribution<long double> myNonCentralChiSquared(d, lambda*exp(-p->k*delta)*vu);
 val = log(pdf(myNonCentralChiSquared,vt/lambda)/lambda+ep);
 //cout<<"chisq: "<<val<<" "<<pdf(myNonCentralChiSquared,vt/lambda)<<endl;
 
 return val;
 }
*/
