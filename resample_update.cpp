#include "include.h"
#include <iostream>
#include <random>
#include <vector>


using namespace std;

void update_norm_weights(para* p){
    para* temp = p;
    double sum = 0;
    double min = 1000000000;
    while(temp != NULL){
        if(sum/temp->weight<1000000){
            sum += temp->weight;
            if(min>temp->weight)
                min = temp->weight;
        }
        else{
            temp->weight = 0;
        }
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
    ADAPT_WEIGHT = min;
    return ;
}
double ESS_0(){
    para* temp = head;
    double num = 0;
    double den = 0;
    while(temp!=NULL){
        num += temp->weight;
        den += pow(temp->weight,2);
        temp = temp->next;
    }
    return pow(num,2)/den;
}
double ESS(double zeta,double prev_zeta){ // to be used after initialization ie iteration>1
    double zeta_diff = zeta - prev_zeta;
    para* temp = head;
    double num = 0;
    double den = 0;
    double like;
    int j =0;
    while(temp != NULL){
        num += exp(untempered_lik[j]*zeta_diff+log(temp->norm_weight));
        den += exp((untempered_lik[j]*zeta_diff+log(temp->norm_weight))*2);
        temp = temp->next;
        j++;
        cout<<"num: "<<num<<" den: "<<den<<endl;
    }

    return exp(log(num*2))/den;
}

double find_new_zeta(double prev_zeta,double u_zeta,double l_zeta,double curr_ESS){
    double curr_zeta = 0.5*(u_zeta+l_zeta);
    double new_ESS = ESS(curr_zeta,prev_zeta);
    cout<<new_ESS<<" curr_zeta: "<<curr_zeta<<endl;
    double diff = new_ESS - 0.95*curr_ESS;
    if(abs(diff)<0.00001)
        return curr_zeta;
    if(diff>0)
        return find_new_zeta(prev_zeta,curr_zeta,l_zeta,curr_ESS);
    return find_new_zeta(prev_zeta,u_zeta,curr_zeta,curr_ESS);
    
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
    new_curr->sum_post = old->sum_post;
    new_curr->weight = 1;
    new_curr->norm_weight = 1.0/no_particles;
    new_curr->cum_norm_weight = 1.0*(i+1)/no_particles;
    new_curr->v = new double[total];
    new_curr->z = new double[total];
    new_curr->posterior = new double[total];
    new_curr->post_lik = new double[total];
    new_curr->post_z = new double[total];
    for(int i=0;i<total;i++){
        new_curr->v[i] = old->v[i];
        new_curr->z[i] = old->z[i];
        new_curr->posterior[i] = old->posterior[i];
        new_curr->post_lik[i] = old->post_lik[i];
        new_curr->post_z[i] = old->post_z[i];
    }
    return ;
}

void destroy(para* p){
    para* next;
    next = p->next;
    while(next != NULL){
        delete p->z;
        delete p->v;
        delete p->post_lik;
        delete p->post_z;
        delete p->posterior;
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

void resample_tester(para* p){
    para* temp = head;
    while(temp != NULL){
        temp->weight = uniform();
        temp = temp->next;
    }
    update_norm_weights(p);
    return ;
}




