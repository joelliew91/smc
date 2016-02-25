#include <iostream>
#include "include.h"

using namespace std;

void update_para(double zeta){
    for(int i =0;i<no_particles;i++){
        
        update_mu(head,kernel->mu,acc,i,zeta);
        update_gam(head,kernel->gam,acc,i,zeta);
        update_sj(head,kernel->sigma_j,acc,i,zeta);
        update_rho(head,kernel->rho,acc,i,zeta);
        update_k(head,kernel->k,acc,i,zeta);
        update_vp(head,kernel->v_p,acc,i,zeta);
        update_sv(head,kernel->sigma_v,acc,i,zeta);
        update_lat_z(head,kernel,acc,i);
        update_lat_v(head,kernel,acc,i);
    }
    return ;
}

void update_lat_v(para* curr,para* kernel,para* acc,int j){
    for(int i=1;i<total;i++){
        double old_post = curr[j].post_lik[i];
        if(i !=1)
            old_post +=curr[j].post_lik[i-1];
        if(i != total-1)
            old_post += curr[j].post_lik[i+1];
        double temp = curr[j].v[i];
        curr[j].v[i] = exp(log(curr[j].v[i]) + normal()*kernel->v[0]);
        
        double new_post = likelihood(1.0,curr[j].z[i],curr[j].v[i],curr[j].v[i-1],price[i],price[i-1],curr[j]);
        
        if(i != 1)
            new_post +=likelihood(1.0,curr[j].z[i-1],curr[j].v[i-1],curr[j].v[i-2],price[i-1],price[i-2],curr[j]);
        
        if(i != total-1)
            new_post += likelihood(1.0,curr[j].z[i+1],curr[j].v[i+1],curr[j].v[i],price[i+1],price[i],curr[j]);
        
        double R = new_post - old_post + curr[j].v[i] - temp;
        double u = uniform();
        //cout<<"R:"<<R<<endl;
        if(exp(R)>u){
            curr[j].sum_post += new_post - old_post;
            curr[j].sum_lik += new_post - old_post;
            
            curr[j].post_lik[i] = likelihood(1.0,curr[j].z[i],curr[j].v[i],curr[j].v[i-1],price[i],price[i-1],curr[j]);
            curr[j].posterior[i] = curr[j].post_z[i] + curr[j].post_lik[i];
            
            if(i !=1){
                curr[j].post_lik[i-1] = likelihood(1.0,curr[j].z[i-1],curr[j].v[i-1],curr[j].v[i-2],price[i-1],price[i-2],curr[j]);
                curr[j].posterior[i-1] = curr[j].post_z[i-1] + curr[j].post_lik[i-1];
            }
            if(i != total-1){
                curr[j].post_lik[i+1] = likelihood(1.0,curr[j].z[i+1],curr[j].v[i+1],curr[j].v[i],price[i+1],price[i],curr[j]);
                curr[j].posterior[i+1] = curr[j].post_z[i+1] + curr[j].post_lik[i+1];
            }
            acc->v[0] += 1.0;
        }
        else{
            curr[j].v[i] = temp;
        }
    }
    return ;
}


void update_lat_z(para* curr,para* kernel,para* acc,int j){
    for(int i=1;i<total;i++){
        double old_post = curr[j].posterior[i];
        double temp = curr[j].z[i];
        curr[j].z[i] += normal()*kernel->z[0];
        double new_post = variance_gamma(curr[j].z[i],curr[j]) + likelihood(1.0,curr[j].z[i],curr[j].v[i],curr[j].v[i-1],price[i],price[i-1],curr[j]);
        double R = new_post - old_post;
        double u = uniform();
        //cout<<"R:"<<R<<endl;
        if(exp(R)>u){
            curr[j].sum_post = curr[j].sum_post + new_post - old_post;
            curr[j].post_z[i] = variance_gamma(curr[j].z[i],curr[j]);
            curr[j].sum_lik -= curr[j].post_lik[i];
            curr[j].post_lik[i] = likelihood(1.0,curr[j].z[i],curr[j].v[i],curr[j].v[i-1],price[i],price[i-1],curr[j]);
            curr[j].sum_lik += curr[j].post_lik[i];
            curr[j].posterior[i] = curr[j].post_z[i]+curr[j].post_lik[i];
            acc->z[0] += 1.0;
        }
        else{
            curr[j].z[i] = temp;
        }
    }
    
    return ;
}


void update_sv(para* curr,double sd,para* acc,int i,double zeta){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double* sum_lik = new double[1];
    sum_lik[0] = 0;
    double old_post = curr[i].sum_post+log(0.5*exp(curr[i].sigma_v)/sqrt(2*0.5*exp(curr[i].sigma_v)/exp(curr[i].k)*0.5*exp(curr[i].sigma_v)/exp(curr[i].v_p)-exp(curr[i].sigma_v)))+(zeta-1)*curr[i].sum_lik;
    
    double temp = curr[i].sigma_v;
    curr[i].sigma_v += normal()*sd;

    double new_post = mcmc_posterior(curr[i],new_post_z,new_post_lik,new_posterior,sum_lik);
    
    
    double R = new_post+(zeta-1)*sum_lik[0] - old_post+log(0.5*exp(curr[i].sigma_v)/sqrt(2*0.5*exp(curr[i].sigma_v)/exp(curr[i].k)*0.5*exp(curr[i].sigma_v)/exp(curr[i].v_p)-exp(curr[i].sigma_v)));
    
    double u = uniform();
    
    if(exp(R)>u){
        curr[i].sum_post = new_post;
        delete curr[i].post_z;
        delete curr[i].post_lik;
        delete curr[i].posterior;
        curr[i].post_z = new_post_z;
        curr[i].post_lik = new_post_lik;
        curr[i].posterior = new_posterior;
        curr[i].sum_lik = sum_lik[0];
        acc->sigma_v += 1.0;
    }
    else{
        curr[i].sigma_v = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }
    delete sum_lik;
    
    return ;
}


void update_vp(para* curr,double sd,para* acc,int i,double zeta){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double* sum_lik = new double[1];
    sum_lik[0]=0;
    
    double old_post = curr[i].sum_post+exp(curr[i].v_p)+(zeta-1)*curr[i].sum_lik;
    double temp = curr[i].v_p;
    curr[i].v_p += normal()*sd;
    
    double new_post = mcmc_posterior(curr[i],new_post_z,new_post_lik,new_posterior,sum_lik);
    double R = new_post +(zeta-1)*sum_lik[0] - old_post +exp(curr[i].v_p);
    double u = uniform();
    if(exp(R)>u){
        curr[i].sum_post = new_post;
        delete curr[i].post_z;
        delete curr[i].post_lik;
        delete curr[i].posterior;
        curr[i].post_z = new_post_z;
        curr[i].post_lik = new_post_lik;
        curr[i].posterior = new_posterior;
        curr[i].sum_lik = sum_lik[0];
        acc->v_p += 1.0;
    }
    else{
        curr[i].v_p = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }
    
    delete sum_lik;
    return ;
}


void update_k(para* curr,double sd,para* acc,int i,double zeta){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double* sum_lik = new double[1];
    sum_lik[0] = 0;
    
    double old_post = curr[i].sum_post+exp(curr[i].k)+(zeta-1)*curr[i].sum_lik;
    double temp = curr[i].k;
    curr[i].k += normal()*sd;
    
    double new_post = mcmc_posterior(curr[i],new_post_z,new_post_lik,new_posterior,sum_lik);
    double R = new_post +(zeta-1)*sum_lik[0]- old_post+exp(curr[i].k);
    double u = uniform();
    if(exp(R)>u){
        curr[i].sum_post = new_post;
        delete curr[i].post_z;
        delete curr[i].post_lik;
        delete curr[i].posterior;
        curr[i].post_z = new_post_z;
        curr[i].post_lik = new_post_lik;
        curr[i].posterior = new_posterior;
        curr[i].sum_lik = sum_lik[0];
        acc->k += 1.0;
    }
    else{
        curr[i].k = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }
    delete sum_lik;
    return ;
}

void update_rho(para* curr,double sd,para* acc,int i,double zeta){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double* sum_lik = new double[1];
    sum_lik[0] = 0;
    
    double old_post = curr[i].sum_post+log(curr[i].rho+1)+(zeta-1)*curr[i].sum_lik;
    double temp = curr[i].rho;
    curr[i].rho = exp(log((1+curr[i].rho)/(1-curr[i].rho)) + normal()*sd);
    curr[i].rho = (curr[i].rho-1)/(curr[i].rho+1);
    double new_post = mcmc_posterior(curr[i],new_post_z,new_post_lik,new_posterior,sum_lik);
    double R = new_post+(zeta-1)*sum_lik[0] - old_post+log(curr[i].rho+1);
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(exp(R)>u){
        curr[i].sum_post = new_post;
        delete curr[i].post_z;
        delete curr[i].post_lik;
        delete curr[i].posterior;
        curr[i].post_z = new_post_z;
        curr[i].post_lik = new_post_lik;
        curr[i].posterior = new_posterior;
        curr[i].sum_lik = sum_lik[0];
        acc->rho += 1.0;
    }
    else{
        curr[i].rho = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }
    
    delete sum_lik;
    return ;
}

void update_sj(para* curr,double sd,para* acc,int i,double zeta){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double* sum_lik = new double[1];
    sum_lik[0] = 0;
    
    double old_post = curr[i].sum_post+curr[i].sigma_j+(zeta-1)*curr[i].sum_lik;
    
    double temp = curr[i].sigma_j;
    curr[i].sigma_j = exp(log(curr[i].sigma_j) + normal()*sd);
    double new_post = mcmc_posterior(curr[i],new_post_z,new_post_lik,new_posterior,sum_lik);
    double R = new_post +(zeta-1)*sum_lik[0] - old_post +curr[i].sigma_j;
    
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(exp(R)>u){
        curr[i].sum_post = new_post;
        delete curr[i].post_z;
        delete curr[i].post_lik;
        delete curr[i].posterior;
        curr[i].post_z = new_post_z;
        curr[i].post_lik = new_post_lik;
        curr[i].posterior = new_posterior;
        curr[i].sum_lik = sum_lik[0];
        acc->sigma_j += 1.0;
    }
    else{
        curr[i].sigma_j = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }
    delete sum_lik;
    return ;
}


void update_gam(para* curr,double sd,para* acc,int i,double zeta){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double* sum_lik = new double[1];
    sum_lik[0] = 0;
    
    double old_post = curr[i].sum_post+(zeta-1)*curr[i].sum_lik;

    double temp = curr[i].gam;
    curr[i].gam += normal()*sd;
    double new_post = mcmc_posterior(curr[i],new_post_z,new_post_lik,new_posterior,sum_lik);
    double R = new_post+(zeta-1)*sum_lik[0] - old_post;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(exp(R)>u){
        curr[i].sum_post = new_post;
        delete curr[i].post_z;
        delete curr[i].post_lik;
        delete curr[i].posterior;
        curr[i].post_z = new_post_z;
        curr[i].post_lik = new_post_lik;
        curr[i].posterior = new_posterior;
        curr[i].sum_lik = sum_lik[0];
        acc->gam += 1.0;
    }
    else{
        curr[i].gam = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;

    }
    
    delete sum_lik;
    return ;
}


void update_mu(para* curr,double sd,para* acc,int i,double zeta){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double* sum_lik = new double[1];
    sum_lik[0] = 0;
    
    double old_post = curr[i].sum_post+(zeta-1)*curr[i].sum_lik;

    double temp = curr[i].mu;
    curr[i].mu += normal()*sd;
    double new_post = mcmc_posterior(curr[i],new_post_z,new_post_lik,new_posterior,sum_lik);
    double R = new_post+(zeta-1)*sum_lik[0] - old_post;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(exp(R)>u){
        //cout<<"acc ";
        curr[i].sum_post = new_post;
        delete curr[i].post_z;
        delete curr[i].post_lik;
        delete curr[i].posterior;
        curr[i].post_z = new_post_z;
        curr[i].post_lik = new_post_lik;
        curr[i].posterior = new_posterior;
        curr[i].sum_lik = sum_lik[0];
        acc->mu += 1.0;
        
    }
    else{
        //cout<<"rej ";
        curr[i].mu = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }
    delete sum_lik;
    
    return ;
}





