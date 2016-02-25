#include "include.h"
#include <iostream>
#include <random>
#include <vector>


using namespace std;

void update_norm_weights(double zeta,double prev_zeta){
    para* p = head;
    double sum=0;double cum_sum = 0;
    double zeta_diff = zeta - prev_zeta;
    
    for(int i=0;i<no_particles;i++){
        p[i].weight = exp(p[i].sum_lik*zeta_diff)*p[i].norm_weight;
        sum += p[i].weight;
    }
    p = head;
    for(int i=0;i<no_particles;i++){
        p[i].norm_weight = p[i].weight/sum;
        cum_sum += p[i].norm_weight;
        p[i].cum_norm_weight = cum_sum;
    }

    return ;
}

double ESS_0(){
    para* temp = head;
    double num = 0;
    double den = 0;
    for(int i=0;i<no_particles;i++){
        num += temp[i].weight;
        den += pow(temp[i].weight,2);
    }
    return num*num/den;
}

void copy_particle(para* new_curr,para old,int j){
    new_curr[j].mu = old.mu;
    new_curr[j].gam = old.gam;
    new_curr[j].sigma_j = old.sigma_j;
    new_curr[j].rho = old.rho;
    new_curr[j].k = old.k;
    new_curr[j].v_p = old.v_p;
    new_curr[j].sigma_v = old.sigma_v;
    new_curr[j].sum_post = old.sum_post;
    new_curr[j].sum_lik = old.sum_lik;
    new_curr[j].weight = 1;
    new_curr[j].norm_weight = 1.0/no_particles;
    new_curr[j].cum_norm_weight = 1.0*(j+1)/no_particles;
    new_curr[j].v = new double[total];
    new_curr[j].z = new double[total];
    new_curr[j].posterior = new double[total];
    new_curr[j].post_lik = new double[total];
    new_curr[j].post_z = new double[total];
    double sum = 0;
    for(int i=0;i<total;i++){
        new_curr[j].v[i] = old.v[i];
        new_curr[j].z[i] = old.z[i];
        new_curr[j].posterior[i] = old.posterior[i];
        new_curr[j].post_lik[i] = old.post_lik[i];
        new_curr[j].post_z[i] = old.post_z[i];
    }
    return ;
}

void destroy(para* p){
    for(int i=0;i<no_particles;i++){
        delete p[i].z;
        delete p[i].v;
        delete p[i].post_lik;
        delete p[i].post_z;
        delete p[i].posterior;
    }
    delete p;
    p = NULL;
    return ;
}

void resample(){

    para* new_head = new para[no_particles];
    para* temp = head;
    double u;int j;
    for(int i=0;i<no_particles;i++){
        j=0;
        u = uniform();
        while(u>temp[j].cum_norm_weight && j<no_particles){
            j++;
        }
        //cout<<u<<" hey "<<temp->cum_norm_weight<<endl;
        copy_particle(new_head,temp[j],i);
        //cout<<"done copy"<<endl;
        
    }
    destroy(head);
    head = new_head;
    
    return ;
}

double ESS(double zeta,double prev_zeta){ // to be used after initialization ie iteration>1
    double zeta_diff = zeta - prev_zeta;
    para* temp = head;
    double num = 0;double temp_w;
    double den = 0;
    double like;

    for(int i=0;i<no_particles;i++){
        temp_w = exp(temp[i].sum_lik*zeta_diff)*temp[i].norm_weight;
        num += temp_w;
        den += temp_w*temp_w;
    }
    return num*num/den;
}

double find_new_zeta(double prev_zeta,double u_zeta,double l_zeta,double curr_ESS){
    double curr_zeta = 0.5*(u_zeta+l_zeta);
    double new_ESS = ESS(curr_zeta,prev_zeta);
    //cout<<new_ESS<<" curr_zeta: "<<curr_zeta<<endl;
    double diff = new_ESS - 0.9*curr_ESS;
    if(abs(curr_zeta-1.0)<TOL) return 1.0;
    if(abs(curr_zeta-l_zeta)<TOL) return l_zeta;
    if(abs(diff)<10)
        return curr_zeta;
    if(diff>0)
        return find_new_zeta(prev_zeta,u_zeta,curr_zeta,curr_ESS);
    return find_new_zeta(prev_zeta,curr_zeta,l_zeta,curr_ESS);
    
}

/*


void resample_tester(para* p,double zeta,double prev_zeta){
    para* temp = head;
    while(temp != NULL){
        temp->weight = uniform();
        temp = temp->next;
    }
    update_norm_weights(zeta,prev_zeta);
    return ;
}*/




