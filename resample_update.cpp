#include "include.h"
#include <iostream>
#include <random>
#include <vector>

using namespace std;

void update_norm_weights(){
    para* temp = head;
    double sum = 0;
    while(temp != NULL){
        sum += temp->weight;
        temp = temp->next;
    }
    temp = head;
    double cum_sum = 0;
    while(temp != NULL){
        temp->norm_weight = temp->weight/sum;
        cum_sum += temp->norm_weight;
        temp->cum_norm_weight = cum_sum;
        temp = temp->next;
    }
    
    return ;
}

double ESS(double zeta,double prev_zeta){
    double zeta_diff = zeta - prev_zeta;
    para* temp = head;
    double num = 0;
    double den = 0;
    
    while(temp != NULL){
        double lik = 0;
        for(int i=1;i<total;i++)
            lik += likelihood(zeta_diff,temp->z[i],temp->v[i],temp->v[i-1],temp->v_star[i],price[i],price[i-1],temp);

        temp->weight = exp(lik)*temp->norm_weight;
        //cout<<lik<<" "<<temp->weight<<endl;
        num += temp->weight;
        den += pow(temp->weight,2);
        temp = temp->next;
    }
    return pow(num,2)/den;
}

double find_new_zeta(double prev_zeta,double u_zeta,double l_zeta,double curr_ESS){
    double curr_zeta = 0.5*(u_zeta+l_zeta);
    //cout<<"curr_zeta: "<<curr_zeta<<endl;
    double new_ESS = ESS(curr_zeta,prev_zeta);
    double diff = new_ESS - 0.95*curr_ESS;
    if(abs(diff)<TOL)
        return curr_zeta;
    if(diff>0)
        return(find_new_zeta(prev_zeta,u_zeta,curr_zeta,curr_ESS));
    return find_new_zeta(prev_zeta,curr_zeta,l_zeta,curr_ESS);
    
}

void copy_particle(para* new_curr,para* old,int i){
    new_curr->mu = old->mu;
    new_curr->gam = old->gam;
    new_curr->sigma_j = old->sigma_j;
    new_curr->rho = old->rho;
    new_curr->k = old->k;
    new_curr->v_p = old->v_p;
    new_curr->sigma_v = old->sigma_v;
    new_curr->u1 = old->u1;
    new_curr->u2 = old->u2;
    new_curr->u3 = old->u3;
    new_curr->weight = 1;
    new_curr->norm_weight = 1.0/no_particles;
    new_curr->cum_norm_weight = 1.0*(i+1)/no_particles;
    new_curr->v = new double[total];
    new_curr->v_star = new double[total];
    new_curr->w = new double[total];
    new_curr->z = new double[total];
    for(int i=0;i<total;i++){
        new_curr->v[i] = old->v[i];
        new_curr->v_star[i] = old->v_star[i];
        new_curr->z[i] = old->z[i];
        new_curr->w[i] = old->w[i];
    }
    return ;
}

void destroy(para* p){
    para* next;
    next = p->next;
    while(next != NULL){
        delete p;
        p = next;
        next = next->next;
    }
    delete p;
    p = NULL;
    return ;
}
void resample(){
    para* curr = NULL;
    para* new_head = NULL;
    para* temp_new;
    para* temp = head;
    double u;
    for(int i=0;i<no_particles;i++){
        u = uniform();
        while(u>temp->cum_norm_weight && temp->next != NULL){
            temp = temp->next;
        }
        //cout<<u<<" hey "<<temp->cum_norm_weight<<endl;
        temp_new = new para;
        copy_particle(temp_new,temp,i);
        //cout<<"done copy"<<endl;
        if(new_head==NULL){
            curr = temp_new;
            new_head = curr;
            curr->next = NULL;
        }
        else{
            curr->next = temp_new;
            curr = curr->next;
            curr->next = NULL;
        }
        temp = head;
    }
    destroy(head);
    head = new_head;
    
    return ;
}

void resample_tester(){
    para* temp = head;
    while(temp != NULL){
        temp->weight = uniform();
        temp = temp->next;
    }
    update_norm_weights();
    return ;
}




