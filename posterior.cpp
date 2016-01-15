#include "include.h"
#include <vector>
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/bessel.hpp>
#include <complex>
#include <complex_bessel.h>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include "rnglib.h"
#include "ranlib.h"

using namespace std;

double prior(para* p){
    double prior_mu,prior_gam,prior_sigj,prior_rho,prior_com;
    int I = 1;
    prior_mu = -pow(p->mu,2)/2;
    prior_gam = -pow(p->gam,2)/2;
    double rhob = (p->rho+1)/2;
    prior_rho = (1/7-1)*log(rhob);
    prior_sigj = log(exp(-p->sigma_j)/tgamma(1.0));
    double d = 4.0*p->k*p->v_p/pow(p->sigma_v,2);
    if(0.5*d>dmax || 0.5*d<1)
        I = 0;
    prior_com = -p->k+log(I)-p->sigma_v-p->v_p;
    //cout<<prior_mu<<" "<<prior_gam<<" "<<p->rho<<" "<<prior_rho<<" "<<prior_sigj<<" "<<prior_com<<endl;
    return prior_mu+prior_gam+prior_rho+prior_sigj+prior_com;
}

double posterior(double zeta,para* p){
    flag = 1;
    double val1 = prior(p);
    double sum= 0 ;
    for(int i =1;i<total;i++){
        sum = sum+post_no_prior(zeta,p->z[i],p->v[i],p->v[i-1],p->v_star[i],price[i],price[i-1],p->w[i],p);
        //if(!flag)
          //  break;
        
    }
    return val1+sum;
}


double post_no_prior(double zeta,double z,double vt,double vu,double v_star,double yt,double yu,double w, para* p){
    double val2 = likelihood(zeta,z,vt,vu,v_star,yt,yu,p);
    double val3 = variance_gamma(z,p);//error
    //double val4 = aux_g(v_star,w,vt,vu);
    double val5 = transition(vt,vu,p);//error
    return val3;
}

double likelihood(double zeta,double z,double vt,double vu,double v_star,double yt,double yu,para* p){
    double var = (1-pow(p->rho,2))*v_star;
    double v_bar = (vt-vu-p->k*p->v_p*delta+p->k*v_star)/p->sigma_v;
    double mn = yu + p->mu*1/7 + z + p->rho*v_bar;
    double val = -pow(yt-mn,2)/(2*var)-0.5*log(var);
    return val*zeta;
}



double variance_gamma(double z,para* p){
    
    double alpha = delta*1.0;
    double sigma = p->sigma_j*sqrt(delta);
    double val1 = p->gam*pow(z,2)/sigma-log(alpha)*delta/alpha-log(tgamma(1.0*delta/alpha));
    double val2 = (2*log(z)-log(pow(p->gam,2)+2*pow(sigma,2)/alpha))*(delta/(2*alpha)-0.5);
    double va = sqrt(pow(z,2)*(pow(p->gam,2)+2*pow(sigma,2)/alpha))/pow(sigma,2);
    //complex<double> a(va,0);
    complex<double> temp;
    double temp2 = boost::math::cyl_bessel_k(delta/alpha-0.5,va);
    temp = sp_bessel::besselK(delta/alpha-0.5,va);
    cout<<val1<<" "<<val2<<" "<<va<<" "<<temp<<" "<<temp2<<endl;
    if(norm(temp)==0||isnan(va)||isinf(va))
        flag=0;
    double val3 = log(temp.real());
    return val1+val2+val3;
}
/*
double aux_g(double v,double w,double vt,double vu){
    int I = 0;
    if(v>min_v)
        if(v<max_v)
            if(w>min_w)
                if(w<max_w)
                    I=1;
    std::complex<double> val = phi(v,w,vt,vu);
    if(norm(val) == 0)
        flag=0;
    double g = cos(v*w)*val.real()+sin(v*w)*val.imag();
    double aux = abs(g)*I+abs(g)*exp(-v-w)*(1-I);
    
    return(log(aux));
}

std::complex<double> phi(double v,double w,double vt,double vu){
    double de = 1.0*delta;
    double com = -2*pow(sigma_v,2)*w;
    complex<double> nu(pow(k,2),com);
    nu = sqrt(nu);
    double d = 0.5*4*k*v/pow(sigma_v,2)-1;
    complex<double> v1 = sqrt(vt*vu)*4.0*nu*exp(-0.5*nu*de)/(pow(sigma_v,2)*(1.0-exp(-nu*de)));
    complex<double> v2 = sqrt(vt*vu)*4.0*k*exp(-0.5*k*de)/(pow(sigma_v,2)*(1.0-exp(-k*de)));
    complex<double> den = sp_bessel::besselI(d,v2);
    complex<double> val1(0,0);
    if(norm(den)==0||isinf(d)||isnan(d)||isinf(v1.real())||isinf(v1.imag())||isinf(v2.real())||isinf(v2.imag())){
        flag = 0;
        return val1;
    }
    else{
        val1 = sp_bessel::besselI(d,v1)/den;
        val1 = val1*nu*exp(-0.5*(nu-k)*de)*(1.0-exp(-k*de))/(k*(1.0-exp(-nu*de)));
        val1 = val1*exp((vu+vt)/pow(sigma_v,2)*(k*(1.0+exp(-k*de))/(1.0-exp(-k*de))-nu*(1.0+exp(-nu*de))/(1.0-exp(-nu*de))));
    }
    return(val1);
}
*/
double transition(double vt,double vu,para* p){
    double ncp = 4.0*p->k*exp(-p->k*delta)*vu/(pow(p->sigma_v,2)*(1-exp(-p->k*delta)));
    double  d = 4*p->k*p->v_p/pow(p->sigma_v,2);

    double val;
    boost::math::non_central_chi_squared_distribution<long double> myNonCentralChiSquared(d, ncp);
    val = log(pdf(myNonCentralChiSquared,vt))+2*log(p->sigma_v)+log(1-exp(-p->k*delta))-log(p->k);
    /*if(isinf(vt)||isinf(vu)||isnan(vt)||isnan(ncp)){
        val = 0;
        flag = 0;
    }
    else{
        boost::math::non_central_chi_squared_distribution<long double> myNonCentralChiSquared(d, ncp);
        val = pdf(myNonCentralChiSquared,vt);
    }*/
    return val;

}



