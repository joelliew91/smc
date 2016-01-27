#include <iostream>
#include "include.h"

using namespace std;

void adapt_kernel(para* kernel,para* acc){
    kernel->mu = kernel->mu*check_mult(acc->mu);
    kernel->gam = kernel->gam*check_mult(acc->gam);
    kernel->sigma_j = kernel->sigma_j*check_mult(acc->sigma_j);
    kernel->rho = kernel->rho*check_mult(acc->rho);
    kernel->k = kernel->k*check_mult(acc->k);
    kernel->v_p = kernel->v_p*check_mult(acc->v_p);
    kernel->sigma_v = kernel->sigma_v*check_mult(acc->sigma_v);
    for(int i=1;i<total;i++){
        kernel->z[i] = kernel->z[i]*check_mult(acc->z[i]);
        kernel->v[i] = kernel->v[i]*check_mult(acc->v[i]);
    }
    return ;
}

double check_mult(double acc_rate){
    double val = 1.0;
    if(acc_rate>0.85)
        val = val*5.0;
    else
        if(acc_rate<0.15)
            val = val/5;
    return val;
}
para* acc_init(){
    para* acc = new para;
    acc->mu = 0;
    acc->gam = 0;
    acc->sigma_j = 0;
    acc->rho = 0;
    acc->k = 0;
    acc->v_p = 0;
    acc->sigma_v = 0;
    acc->v = new double[total];
    acc->z = new double[total];
    for(int i=0;i<total;i++){
        acc->v[i]=0;
        acc->z[i]=0;
    }
    return acc;
}

void acc_end(para* acc){
    delete acc->v;
    delete acc->z;
    delete acc;
    return ;
}

para* reset_kernel(para* ker){
    delete ker->v;
    delete ker->z;
    delete ker;
    return set_kernel();
}

para* set_kernel(){
    para* kernel = new para;
    kernel->mu = 0;
    kernel->gam = 0;
    kernel->sigma_j = 0;
    kernel->rho = 0;
    kernel->k = 0;
    kernel->v_p = 0;
    kernel->sigma_v = 0;
    kernel->v = new double[total];
    kernel->z = new double[total];
    for(int i=0;i<total;i++){
        kernel->v[i]=0;
        kernel->z[i]=0;
    }
    para* curr = head;
    while(curr != NULL){
        kernel->mu += curr->mu;
        kernel->gam += curr->gam;
        kernel->sigma_j += log(curr->sigma_j);
        kernel->rho += log((1+curr->rho)/(1-curr->rho));
        kernel->k += log(curr->k - curr->u1);
        kernel->v_p += log(curr->v_p - curr->u2);
        kernel->sigma_v += log(curr->sigma_v/(curr->u3 - curr->sigma_v));
        //cout<<kernel->k<<" "<<kernel->v_p<<" "<<kernel->sigma_v<<endl;
        for(int i=0;i<total;i++){
            kernel->v[i] += log(curr->v[i]);
            kernel->z[i] += curr->z[i];
        }
        curr = curr->next;
    }
    
    curr = head;
    para* m = kernel;
    
    kernel = new para;
    kernel->mu = 0;
    kernel->gam = 0;
    kernel->sigma_j = 0;
    kernel->rho = 0;
    kernel->k = 0;
    kernel->v_p = 0;
    kernel->sigma_v = 0;
    kernel->v = new double[total];
    kernel->z = new double[total];
    for(int i=0;i<total;i++){
        kernel->v[i]=0;
        kernel->z[i]=0;
    }

    while(curr != NULL){
        kernel->mu += pow(curr->mu - m->mu/no_particles,2);
        kernel->gam += pow(curr->gam - m->gam/no_particles,2);
        kernel->sigma_j += pow(log(curr->sigma_j) - m->sigma_j/no_particles,2);
        kernel->rho += pow(log((1+curr->rho)/(1-curr->rho)) - m->rho/no_particles,2);
        kernel->k += pow(log(curr->k-curr->u1) - m->k/no_particles,2);
        kernel->v_p += pow(log(curr->v_p-curr->u2) - m->v_p,2);
        kernel->sigma_v += pow(log(curr->sigma_v/(curr->u3-curr->sigma_v)) - m->sigma_v,2);
        
        for(int i=0;i<total;i++){
            kernel->v[i] += pow(log(curr->v[i])-m->v[i]/no_particles,2);
            kernel->z[i] += pow(curr->z[i]-m->z[i]/no_particles,2);
        }
        curr = curr->next;
    }
    kernel->mu = sqrt(kernel->mu/(no_particles-1));
    kernel->gam = sqrt(kernel->gam/(no_particles-1));
    kernel->sigma_j = sqrt(kernel->sigma_j/(no_particles-1));
    kernel->rho = sqrt(kernel->rho/(no_particles-1));
    kernel->k = sqrt(kernel->k/(no_particles-1));
    kernel->v_p = sqrt(kernel->v_p/(no_particles-1));
    kernel->sigma_v = sqrt(kernel->sigma_v/(no_particles-1));
    for(int i=0;i<total;i++){
        kernel->v[i] = sqrt(kernel->v[i]/(no_particles-1));
        kernel->z[i] = sqrt(kernel->z[i]/(no_particles-1));
    }
    delete m->v;
    delete m->z;
    delete m;
    kernel->next = NULL;
	return kernel;
}
