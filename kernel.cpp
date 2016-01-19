#include <iostream>
#include "include.h"

using namespace std;

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
    kernel->v_star = new double[total];
    kernel->w = new double[total];
    kernel->z = new double[total];
    for(int i=0;i<total;i++){
        kernel->v[i]=0;
        kernel->v_star[i] = 0;
        kernel->w[i] = 0;
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
            kernel->v_star[i] += log(curr->v_star[i]);
            kernel->w[i] += log(curr->w[i]);
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
    kernel->v_star = new double[total];
    kernel->w = new double[total];
    kernel->z = new double[total];
    for(int i=0;i<total;i++){
        kernel->v[i]=0;
        kernel->v_star[i] = 0;
        kernel->w[i] = 0;
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
            kernel->v_star[i] += pow(log(curr->v_star[i])-m->v_star[i]/no_particles,2);
            kernel->w[i] += pow(log(curr->w[i])-m->w[i]/no_particles,2);
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
        kernel->v_star[i] = sqrt(kernel->v_star[i]/(no_particles-1));
        kernel->w[i] = sqrt(kernel->w[i]/(no_particles-1));
        kernel->z[i] = sqrt(kernel->z[i]/(no_particles-1));
    }
    delete m;
    kernel->next = NULL;
	return kernel;
}
