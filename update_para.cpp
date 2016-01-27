#include <iostream>
#include "include.h"

using namespace std;

void update_para(double zeta,para* kernel,para* acc){
    para* temp = head;
    //print(head);
    while(temp != NULL){
        update_mu(zeta,temp,kernel->mu,acc);
        //cout<<"mu"<<endl;
        update_gam(zeta,temp,kernel->gam,acc);
        //cout<<"gam"<<endl;
        update_sj(zeta,temp,kernel->sigma_j,acc);
        //cout<<"sj"<<endl;
        update_rho(zeta,temp,kernel->rho,acc);
        //cout<<"rho"<<endl;
        update_k(zeta,temp,kernel->k,acc);
        //cout<<"k"<<endl;
        update_vp(zeta,temp,kernel->v_p,acc);
        //cout<<"vp"<<endl;
        update_sv(zeta,temp,kernel->sigma_v,acc);
        //update s_v
        update_lat_z(zeta,temp,kernel,acc);
        update_lat_v(zeta,temp,kernel,acc);
        //update lat_v
        //update lat_vs
        //update lat_z
        //update lat_w
        //cout<<temp->k-temp->u1<<" "<<temp->v_p - temp->u2<<" "<<temp->u3-temp->sigma_v<<endl;
        //cout<<temp->mu<<" "<<temp->gam<<" "<<temp->sigma_j<<" "<<temp->rho<<" "<<temp->k<<" "<<temp->v_p<<endl;
        temp = temp->next;
    }
    return ;
}

void update_lat_v(double zeta,para* curr,para* kernel,para* acc){
    for(int i=1;i<total-1;i++){
        double old_post = curr->post_lik[i]+curr->post_lik[i+1];
        if(i>1)
            old_post += curr->post_lik[i-1];
        double temp = curr->v[i];
        curr->v[i] = exp(log(curr->v[i]) + normal()*kernel->v[i]);
        double new_post = likelihood(zeta,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr);
        new_post += likelihood(zeta,curr->z[i+1],curr->v[i+1],curr->v[i],price[i+1],price[i],curr);
        if(i>1)
            new_post += likelihood(zeta,curr->z[i-1],curr->v[i-1],curr->v[i-2],price[i-1],price[i-2],curr);
        double R = new_post - old_post + curr->v[i] - temp;
        double u = uniform();
        //cout<<"R:"<<R<<endl;
        if(R<log(u)||R!=R)
            curr->v[i] = temp;
        else{
            curr->sum_post -= curr->post_lik[i];
            curr->post_lik[i] = likelihood(zeta,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr);
            curr->sum_post += curr->post_lik[i];
            acc->v[i] += 1.0;
        }
    }
    
    return ;
}

void update_lat_z(double zeta,para* curr,para* kernel,para* acc){
    for(int i=1;i<total;i++){
        double old_post = curr->post_z[i] + curr->post_lik[i];
        double temp = curr->z[i];
        curr->z[i] = curr->z[i] + normal()*kernel->z[i];
        double new_post = variance_gamma(curr->z[i],curr) + likelihood(zeta,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr);
        double R = new_post - old_post;
        double u = uniform();
    //cout<<"R:"<<R<<endl;
        if(R<log(u)||R!=R)
            curr->z[i] = temp;
        else{
            curr->sum_post += new_post - old_post;
            curr->post_z[i] = new_post;
            curr->post_lik[i] = likelihood(zeta,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr);
            acc->z[i] += 1.0;
        }
    }
    
    return ;
}

void update_sv(double zeta,para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];

    double old_post = curr->sum_post+curr->sigma_v;
    double temp = curr->sigma_v;
    curr->sigma_v = exp(log(curr->sigma_v/(curr->u3-curr->sigma_v)) + normal()*sd);
    curr->u3 = sqrt(2*curr->k*curr->v_p);
    curr->sigma_v = curr->sigma_v*curr->u3/(curr->sigma_v+1);
    double new_post = mcmc_posterior(zeta,curr,new_post_z,new_post_lik);
    double R = new_post - old_post+curr->sigma_v;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(R<log(u)||R!=R)
        curr->sigma_v = temp;
    else{
        curr->sum_post = new_post;
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }
        acc->sigma_v += 1.0;
    }
    delete new_post_lik;delete new_post_z;
    
    return ;
}

void update_vp(double zeta,para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    
    double old_post = curr->sum_post+curr->v_p;
    double temp = curr->v_p;
    curr->v_p = exp(log(curr->v_p-curr->u2) + normal()*sd);
    curr->u2 = pow(curr->sigma_v,2)/(2*curr->k);
    curr->v_p += curr->u2;
    double new_post = mcmc_posterior(zeta,curr,new_post_z,new_post_lik);
    double R = new_post - old_post +curr->v_p;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(R<log(u)||R!=R)
        curr->v_p = temp;
    else{
        curr->sum_post = new_post;
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }
        acc->v_p += 1.0;
    }
    delete new_post_lik;delete new_post_z;
    
    return ;
}

void update_k(double zeta,para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    
    double old_post = curr->sum_post+curr->k;
    double temp = curr->k;
    curr->k = exp(log(curr->k-curr->u1) + normal()*sd);
    curr->u1 = pow(curr->sigma_v,2)/(2*curr->v_p);
    curr->k += curr->u1;
    double new_post = mcmc_posterior(zeta,curr,new_post_z,new_post_lik);
    double R = new_post - old_post+curr->k;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(R<log(u)||R!=R)
        curr->k = temp;
    else{
        curr->sum_post = new_post;
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }
        acc->k += 1.0;
    }
    delete new_post_lik;delete new_post_z;
    
    return ;
}

void update_rho(double zeta,para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];

    double old_post = curr->sum_post+log(curr->rho+1);
    double temp = curr->rho;
    curr->rho = exp(log((1+curr->rho)/(1-curr->rho)) + normal()*sd);
    curr->rho = (curr->rho-1)/(curr->rho+1);
    double new_post = mcmc_posterior(zeta,curr,new_post_z,new_post_lik);
    double R = new_post - old_post+log(curr->rho+1);
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(R<log(u)||R!=R)
        curr->rho = temp;
    else{
        curr->sum_post = new_post;
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }
        acc->rho += 1.0;
    }
    delete new_post_lik;delete new_post_z;
    
    return ;
}

void update_sj(double zeta,para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    
    double old_post = curr->sum_post+curr->sigma_j;
    double temp = curr->sigma_j;
    curr->sigma_j = exp(log(curr->sigma_j) + normal()*sd);
    double new_post = mcmc_posterior(zeta,curr,new_post_z,new_post_lik);
    double R = new_post - old_post+curr->sigma_j;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(R<log(u)||R!=R)
        curr->sigma_j = temp;
    else{
        curr->sum_post = new_post;
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }
        acc->sigma_j += 1.0;
    }
    delete new_post_lik;delete new_post_z;
    
    return ;
}

void update_gam(double zeta,para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    
    double old_post = curr->sum_post;
    double temp = curr->gam;
    curr->gam = curr->gam + normal()*sd;
    double new_post = mcmc_posterior(zeta,curr,new_post_z,new_post_lik);
    double R = new_post - old_post;
    double u = uniform();
   //cout<<"R:"<<R<<endl;
    if(R<log(u)||R!=R)
        curr->gam = temp;
    else{
        curr->sum_post = new_post;
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }
        acc->gam += 1.0;
    }
    delete new_post_lik;delete new_post_z;
    
    return ;
}

void update_mu(double zeta,para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    
    double old_post = curr->sum_post;
    double temp = curr->mu;
    curr->mu = curr->mu + normal()*sd;
    double new_post = mcmc_posterior(zeta,curr,new_post_z,new_post_lik);
    double R = new_post - old_post;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(R<log(u)||R!=R){
        curr->mu = temp;
        //cout<<"mu:"<<curr->mu<<" old_post:"<<old_post<<" new_post"<<new_post<<endl;
    }
    else{
        //cout<<"mu:"<<curr->mu<<" old_post:"<<old_post<<" new_post"<<new_post<<endl;
        curr->sum_post = new_post;
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }
        acc->mu += 1.0;
    }
    delete new_post_lik;delete new_post_z;
    
    return ;
}




