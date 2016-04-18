#include "include.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <complex>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <complex_bessel.h>
#define PII 3.141592654

using namespace std;

//posterior distribution
double likelihood(double zeta,double z,double vt,double vu,double yt,double yu,para p){
    double y1 = p.k;
    double y2 = p.v_p;
    double y3 = p.sigma_v;
    double k = 0.5*exp(y3)/exp(y1);
    double v_p = 0.5*exp(y3)/exp(y2);
    double sigma_v = sqrt(2*p.k*p.v_p-exp(y3));
    
    double var = (1-pow(p.rho,2))*vu*1.0/250;
    double norm_v = (vt-(vu+k*(v_p-vu)*1.0/250))/(sigma_v*sqrt(vu*1.0/250));
    double norm_y = (yt-(yu + p.mu*1.0/250 + z))/sqrt(vu*1.0/250);
    double val = -0.5*(norm_y*norm_y+norm_v*norm_v-2*p.rho*norm_v*norm_y)/(1-p.rho*p.rho)+log(1/(sqrt(2*PII*(1-p.rho*p.rho))*sigma_v*vu*1.0/250));
    
    
    return val*zeta;
}

double variance_gamma(double z,para p){
    
    double alpha = delta*1.0;
    double sigma = p.sigma_j*sqrt(1.0/250);
    double val = -log(sigma);
    sigma  = sigma*sigma;
    //cout<<"be:"<<sigma<<" "<<p.sigma_j<<" "<<sqrt(1.0/250)<<" "<<val<<endl;
    val += p.gam*z/sigma;
    val -= delta/alpha*log(1.0/250);
    //cout<<"1:"<<val<<endl;
    val -= tgamma(1.0*delta/alpha);
    //cout<<"2:"<<val<<endl;
    val += (0.5*delta/alpha-0.25)*log(z*z/(p.gam*p.gam+2*sigma/250));
    
    double va = (sqrt(p.gam*p.gam+2*sigma/250)*abs(z))/sigma;
    //cout<<va<<" "<<val<<endl;
    if(abs(alpha-delta)<TOL)
    {
        val += -0.5*log(va) - va;
    }
    else
    {
        val += log(bessik(va,1.0*delta/alpha - 0.5));
    }
    
    return val;//val1+log(val3)+val2;
}

double prior(para p){
    double prior_mu,prior_gam,prior_sigj,prior_rho,prior_sigv,prior_k,prior_vp,w_v,phi_v;
    
    double y1 = p.k;
    double y2 = p.v_p;
    double y3 = p.sigma_v;
    double k = 0.5*exp(y3)/exp(y1);
    double v_p = 0.5*exp(y3)/exp(y2);
    double sigma_v = sqrt(2*k*v_p-exp(y3));
    
    prior_mu = -pow(p.mu,2)/2-log(sqrt(2*PII));
    prior_gam = -pow(p.gam,2)/2-log(sqrt(2*PII));
    
    phi_v = p.rho*sigma_v;
    w_v = sigma_v*sigma_v*(1-p.rho*p.rho);
    prior_rho = log(exp(-0.5/w_v)*pow(w_v,-1.0-1)*pow(0.5,1)/tgamma(1));
    prior_sigv = -pow(phi_v,2)/(2*pow(0.5/w_v,2))-log(sqrt(2*PII*pow(0.5/w_v,2)));
    
    prior_sigj = log(exp(-5/p.sigma_j)*pow(p.sigma_j,-2.5-1)*pow(5,2.5)*tgamma(2.5));
    prior_k = -pow(k,2)/2-log(sqrt(2*PII));
    prior_vp = -pow(v_p,2)/2-log(sqrt(2*PII));
    
    //cout<<prior_mu<<" "<<prior_gam<<" "<<prior_rho<<" "<<prior_sigv<<" "<<prior_sigj<<" "<<prior_k<<" "<<prior_vp<<endl;
    return prior_mu+prior_gam+prior_rho+prior_sigv+prior_sigj+prior_k+prior_vp;
}

//Checking of the posterior values
void update_all(){
    para* p = head;
    double sum,like,val;
    for(int j=0;j<no_particles;j++){
        sum = 0;
        like = 0;
        for(int i=1;i<total;i++){
            val = variance_gamma(p[j].z[i],p[j]);
            if(abs(val-p[j].post_z[i])>TOL){
                cout<<val<<" wrong z "<<p[j].post_z[i]<<endl;
                p[j].post_z[i] = val;
            }
            val = likelihood(1.0,p[j].z[i],p[j].v[i],p[j].v[i-1],price[i],price[i-1],p[j]);
            if(abs(val- p[j].post_lik[i])>TOL){
                cout<<val<<" wrong lik "<<p[j].post_lik[i]<<endl;
                p[j].post_lik[i] = val;
            }
            like += val;
            val = p[j].post_z[i]+p[j].post_lik[i];
            if(abs(val - p[j].posterior[i])>TOL){
                cout<<"id:"<<i<<" corr:"<<val<<" wrong posterior "<<"wrong:"<<p[j].posterior[i]<<endl;
                p[j].posterior[i] = val;
            }
            sum += p[j].posterior[i];
        }
        if(abs(like-p[j].sum_lik)>TOL){
            cout<<like<<" wrong sum_lik "<<p[j].sum_lik<<endl;
            p[j].sum_lik = like;
        }
        if(abs(sum+prior(p[j])- p[j].sum_post)>TOL){
            cout<<sum<<" "<<sum+prior(p[j])<<" wrong sum "<<p[j].sum_post<<endl;
            p[j].sum_post = sum;
        }
    }
    return ;
}

//Liklihood Values
double lik(double zeta,para p){
    double sum = 0;
    for(int i=1;i<total;i++){
        sum += likelihood(zeta,p.z[i],p.v[i],p.v[i-1],price[i],price[i-1],p);
        sum += variance_gamma(p.z[i],p);
    }
    return sum*zeta;
}


//For the MCMC calculation
double mcmc_posterior(para p,double* lat_z,double* lat_lik,double* lat_posterior,double* sum_lik){   //untempered posterior
    double val = 0;
    lat_lik[0] = 0;
    lat_z[0] = 0;
    lat_posterior[0]=0;
    
    for(int i =1;i<total;i++){
        lat_lik[i] = likelihood(1.0,p.z[i],p.v[i],p.v[i-1],price[i],price[i-1],p);
        lat_z[i] = variance_gamma(p.z[i],p);
        lat_posterior[i] = lat_z[i]+lat_lik[i];
        val = val+ lat_posterior[i];
        sum_lik[0] += lat_lik[i];
        //if(lat_lik[i] != lat_lik[i])
        //  cout<<"lik "<<i<<endl;
    }
    
    return val+prior(p);
}
