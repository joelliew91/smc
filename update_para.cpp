#include <iostream>
#include "include.h"

using namespace std;
/*
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
*/
void update_para(double zeta,para* kernel){
    para* temp = head;
    //print(head);
    while(temp != NULL){
        update_mu(zeta,temp,kernel->mu);
        //cout<<"mu"<<endl;
        update_gam(zeta,temp,kernel->gam);
        //cout<<"gam"<<endl;
        update_sj(zeta,temp,kernel->sigma_j);
        //cout<<"sj"<<endl;
        update_rho(zeta,temp,kernel->rho);
        //cout<<"rho"<<endl;
        update_k(zeta,temp,kernel->k);
        //cout<<"k"<<endl;
        update_vp(zeta,temp,kernel->v_p);
        //cout<<"vp"<<endl;
        //update_sv(zeta,temp,kernel->sigma_v);
        //update s_v
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

void update_sv(double zeta,para* curr,double sd){
    double old_post = posterior(zeta,curr)+curr->sigma_v;
    double temp = curr->sigma_v;
    curr->sigma_v = exp(log(curr->sigma_v/(curr->u3-curr->sigma_v)) + normal()*sd);
    curr->u3 = sqrt(2*curr->k*curr->v_p);
    curr->sigma_v = curr->sigma_v*curr->u3/(curr->sigma_v+1);
    double new_post = posterior(zeta,curr)+curr->sigma_v;
    double R = exp(new_post - old_post);
    double u = uniform();
    if(R<u)
        curr->sigma_v = temp;
    
    return ;
}

void update_vp(double zeta,para* curr,double sd){
    double old_post = posterior(zeta,curr)+curr->v_p;
    double temp = curr->v_p;
    curr->v_p = exp(log(curr->v_p-curr->u2) + normal()*sd);
    curr->u2 = pow(curr->sigma_v,2)/(2*curr->k);
    curr->v_p += curr->u2;
    double new_post = posterior(zeta,curr)+curr->v_p;
    double R = exp(new_post - old_post);
    double u = uniform();
    if(R<u)
        curr->v_p = temp;
    
    return ;
}

void update_k(double zeta,para* curr,double sd){
    double old_post = posterior(zeta,curr)+curr->k;
    double temp = curr->k;
    curr->k = exp(log(curr->k-curr->u1) + normal()*sd);
    curr->u1 = pow(curr->sigma_v,2)/(2*curr->v_p);
    curr->k += curr->u1;
    double new_post = posterior(zeta,curr)+curr->k;
    double R = exp(new_post - old_post);
    double u = uniform();
    if(R<u)
        curr->k = temp;
    
    return ;
}

void update_rho(double zeta,para* curr,double sd){
    double old_post = posterior(zeta,curr)+log(curr->rho+1);
    double temp = curr->rho;
    curr->rho = exp(log((1+curr->rho)/(1-curr->rho)) + normal()*sd);
    curr->rho = (curr->rho-1)/(curr->rho+1);
    double new_post = posterior(zeta,curr)+log(curr->rho+1);
    double R = exp(new_post - old_post);
    double u = uniform();
    if(R<u)
        curr->rho = temp;
    
    return ;
}

void update_sj(double zeta,para* curr,double sd){
    double old_post = posterior(zeta,curr)+log(curr->sigma_j);
    double temp = curr->sigma_j;
    curr->sigma_j = exp(log(curr->sigma_j) + normal()*sd);
    double new_post = posterior(zeta,curr)+log(curr->sigma_j);
    double R = exp(new_post - old_post);
    double u = uniform();
    if(R<u)
        curr->sigma_j = temp;
    
    return ;
}

void update_gam(double zeta,para* curr,double sd){
    double old_post = posterior(zeta,curr);
    double temp = curr->gam;
    curr->gam = curr->gam + normal()*sd;
    double new_post = posterior(zeta,curr);
    double R = exp(new_post - old_post);
    double u = uniform();
    if(R<u)
        curr->gam = temp;

    return ;
}

void update_mu(double zeta,para* curr,double sd){
    double old_post = posterior(zeta,curr);
    double temp = curr->mu;
    curr->mu = curr->mu + normal()*sd;
    double new_post = posterior(zeta,curr);
    double R = exp(new_post - old_post);
    double u = uniform();
    if(R<u)
        curr->mu = temp;
    
    return ;
}




