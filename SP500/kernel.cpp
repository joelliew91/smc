#include <iostream>
#include "include.h"

using namespace std;

void print_kernel(para* ker){
    cout<<ker->mu<<" "<<ker->gam<<" "<<ker->sigma_j<<" "<<ker->rho<<" "<<ker->k<<" "<<ker->v_p<<" "<<ker->sigma_v<<" "<<ker->v[0]<<" "<<ker->z[0]<<endl;
    return ;
}
void set_kernel(){
    
    para* m = new para;
    
    m->mu = 0;
    m->gam = 0;
    m->sigma_j = 0;
    m->rho = 0;
    m->k = 0;
    m->v_p = 0;
    m->sigma_v = 0;
    m->v = new double[1];
    m->z = new double[1];
    for(int i=0;i<1;i++){
        m->v[i]=0;
        m->z[i]=0;
    }
    
    


    for(int i=0;i<no_particles;i++){
        m->mu += head[i].mu;
        m->gam += head[i].gam;
        m->sigma_j += log(head[i].sigma_j);
        m->rho += log((1+head[i].rho)/(1-head[i].rho));
        m->k += head[i].k;
        m->v_p += head[i].v_p;
        m->sigma_v += head[i].sigma_v;
        //cout<<curr->k - curr->u1<<" "<<curr->v_p - curr->u2<<" "<<(curr->u3 - curr->sigma_v)<<endl;
        for(int i=0;i<total;i++){
            m->v[0] += log(head[i].v[i]);
            m->z[0] += head[i].z[i];
        }
    }
    
    kernel->mu = 0;
    kernel->gam = 0;
    kernel->sigma_j = 0;
    kernel->rho = 0;
    kernel->k = 0;
    kernel->v_p = 0;
    kernel->sigma_v = 0;
    for(int i=0;i<1;i++){
        kernel->v[i]=0;
        kernel->z[i]=0;
    }
    
    for(int i=0;i<no_particles;i++){
        kernel->mu += pow(head[i].mu - m->mu/no_particles,2);
        kernel->gam += pow(head[i].gam - m->gam/no_particles,2);
        kernel->sigma_j += pow(log(head[i].sigma_j) - m->sigma_j/no_particles,2);
        kernel->rho += pow(log((1+head[i].rho)/(1-head[i].rho)) - m->rho/no_particles,2);
        kernel->k += pow(head[i].k - m->k/no_particles,2);
        kernel->v_p += pow(head[i].v_p - m->v_p/no_particles,2);
        kernel->sigma_v += pow(head[i].sigma_v - m->sigma_v/no_particles,2);
        
        for(int j=0;j<total;j++){
            kernel->v[0] += pow(log(head[i].v[j])-m->v[0]/(total*no_particles),2);
            kernel->z[0] += pow(head[i].z[j]-m->z[0]/(total*no_particles),2);
        }
    }
    
    
    kernel->mu = sqrt(kernel->mu/(no_particles-1));
    kernel->gam = sqrt(kernel->gam/(no_particles-1));
    kernel->sigma_j = sqrt(kernel->sigma_j/(no_particles-1));
    kernel->rho = sqrt(kernel->rho/(no_particles-1));
    kernel->k = sqrt(kernel->k/(no_particles-1));
    kernel->v_p = sqrt(kernel->v_p/(no_particles-1));
    kernel->sigma_v = sqrt(kernel->sigma_v/(no_particles-1));
    
    for(int i=0;i<1;i++){
        kernel->v[i] = sqrt(kernel->v[0]/(total*no_particles-1));
        kernel->z[i] = sqrt(kernel->z[0]/(total*no_particles-1));
    }
    //cout<<m.mu<<" "<<m.gam<<" "<<m.sigma_j<<" "<<m.rho<<" "<<m.k<<" "<<m.v_p<<" "<<m.sigma_v<<" "<<m.v[0]<<" "<<m.z[0]<<endl;
    delete m->v;
    delete m->z;
    delete m;
    return ;
}

void acc_init(){
    acc->mu = 0;
    acc->gam = 0;
    acc->sigma_j = 0;
    acc->rho = 0;
    acc->k = 0;
    acc->v_p = 0;
    acc->sigma_v = 0;
    for(int i=0;i<1;i++){
        acc->v[i]=0;
        acc->z[i]=0;
    }
    return ;
}

void adapt_kernel(){
    if(kernel->sigma_v<TOL) kernel->sigma_v = 1;
    if(kernel->k<TOL) kernel->k = 1;
    if(kernel->v_p<TOL) kernel->v_p= 1;
    if(kernel->rho<TOL) kernel->rho = 1;
    if(kernel->sigma_j<TOL) kernel->sigma_j =1;
    if(kernel->gam<TOL) kernel->gam = 1;
    if(kernel->mu<TOL) kernel->mu = 1;
    if(kernel->v[0]<TOL) kernel->v[0] = 1;
    if(kernel->z[0]<TOL) kernel->z[0] = 1;
    kernel->mu = kernel->mu*check_mult(acc->mu);
    kernel->gam = kernel->gam*check_mult(acc->gam);
    kernel->sigma_j = kernel->sigma_j*check_mult(acc->sigma_j);
    kernel->rho = kernel->rho*check_mult(acc->rho);
    kernel->k = kernel->k*check_mult(acc->k);
    kernel->v_p = kernel->v_p*check_mult(acc->v_p);
    kernel->sigma_v = kernel->sigma_v*check_mult(acc->sigma_v);
    
    for(int i=0;i<1;i++){
        double acc_z = acc->z[0] / total;
        double acc_v = acc->v[0] / total;
        kernel->z[i] = kernel->z[i]*check_mult(acc_z);
        kernel->v[i] = kernel->v[i]*check_mult(acc_v);
    }
    return ;
}


double check_mult(double acc_rate){
    double val = 1.0;
    if((acc_rate/no_particles)>0.85)
        val = val*5.0;
    else{
        if((acc_rate/no_particles)<0.15)
            val = val/5;
        if(val==0)
            val = 0.0001;
    }
    return val;
}
/*




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
    kernel->u1 = 0;
    kernel->u2 = 0;
    kernel->u3 = 0;
    kernel->v = new double[total];
    kernel->z = new double[total];
    kernel->next = NULL;
    for(int i=0;i<1;i++){
        kernel->v[i]=0;
        kernel->z[i]=0;
    }
    para* curr = head;
    while(curr != NULL){
        kernel->mu += curr->mu;
        kernel->gam += curr->gam;
        kernel->sigma_j += log(curr->sigma_j);
        kernel->rho += log((1+curr->rho)/(1-curr->rho));
        kernel->k += curr->k;
        kernel->v_p += curr->v_p;
        //kernel->sigma_v += log(curr->sigma_v/(curr->u3 - curr->sigma_v));
        kernel->sigma_v += curr->sigma_v;
        //cout<<curr->k - curr->u1<<" "<<curr->v_p - curr->u2<<" "<<(curr->u3 - curr->sigma_v)<<endl;
        for(int i=0;i<total;i++){
            kernel->v[0] += log(curr->v[i]);
            kernel->z[0] += curr->z[i];
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
    for(int i=0;i<1;i++){
        kernel->v[i]=0;
        kernel->z[i]=0;
    }

    while(curr != NULL){
        kernel->mu += pow(curr->mu - m->mu/no_particles,2);
        kernel->gam += pow(curr->gam - m->gam/no_particles,2);
        kernel->sigma_j += pow(log(curr->sigma_j) - m->sigma_j/no_particles,2);
        kernel->rho += pow(log((1+curr->rho)/(1-curr->rho)) - m->rho/no_particles,2);
        kernel->k += pow(curr->k - m->k/no_particles,2);
        kernel->v_p += pow(curr->v_p - m->v_p/no_particles,2);
        kernel->sigma_v += pow(curr->sigma_v - m->sigma_v/no_particles,2);
        
        for(int i=0;i<total;i++){
            kernel->v[0] += pow(log(curr->v[i])-m->v[0]/(total*no_particles),2);
            kernel->z[0] += pow(curr->z[i]-m->z[0]/(total*no_particles),2);
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
    for(int i=0;i<1;i++){
        kernel->v[i] = sqrt(kernel->v[0]/(total*no_particles-1));
        kernel->z[i] = sqrt(kernel->z[0]/(total*no_particles-1));
    }
    delete m->v;
    delete m->z;
    delete m;
    kernel->next = NULL;
	return kernel;
}
*/