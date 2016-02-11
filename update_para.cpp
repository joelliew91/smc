#include <iostream>
#include "include.h"

using namespace std;

void update_para(para* kernel,para* acc){
    para* temp = head;
    int count = 1;
    while(temp != NULL){
        update_mu(temp,kernel->mu,acc);
        //cout<<"mu"<<endl;
        update_gam(temp,kernel->gam,acc);
        //cout<<"gam"<<endl;
        update_sj(temp,kernel->sigma_j,acc);
        //cout<<"sj"<<endl;
        update_rho(temp,kernel->rho,acc);
        //cout<<"rho"<<endl;*/
        update_k(temp,kernel->k,acc);
        //cout<<"k"<<endl;
        update_vp(temp,kernel->v_p,acc);
        //cout<<"vp"<<endl;
        update_sv(temp,kernel->sigma_v,acc);
        //update s_v
        update_lat_z(temp,kernel,acc);
        update_lat_v(temp,kernel,acc);
        update_untempered_lik();
        temp = temp->next;
        cout<<count<<endl;
        count++;
    }
    return ;
}

void update_lat_v(para* curr,para* kernel,para* acc){
    int count = 0;
    for(int i=1;i<total;i++){
        double old_post = curr->post_lik[i]+curr->post_lik[i-1];
        if(i != total-1)
            old_post += curr->post_lik[i+1];
        double temp = curr->v[i];
        curr->v[i] = exp(log(curr->v[i]) + normal()*kernel->v[i]);
        double new_post = likelihood(1.0,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr)+likelihood(1.0,curr->z[i-1],curr->v[i-1],curr->v[i-2],price[i-1],price[i-2],curr);
        if(i != total-1)
           new_post += likelihood(1.0,curr->z[i+1],curr->v[i+1],curr->v[i],price[i+1],price[i],curr);

        double R = new_post - old_post + curr->v[i] - temp;
        double u = uniform();
        //cout<<"R:"<<R<<endl;
        if(exp(R)>u){
            curr->sum_post += new_post-old_post;
            curr->post_lik[i] = likelihood(1.0,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr);
            curr->post_lik[i-1] = likelihood(1.0,curr->z[i-1],curr->v[i-1],curr->v[i-2],price[i-1],price[i-2],curr);
            curr->posterior[i] = curr->post_z[i] + curr->post_lik[i];
            curr->posterior[i-1] = curr->post_z[i-1] + curr->post_lik[i-1];
            if(i != total-1){
                curr->post_lik[i+1] = likelihood(1.0,curr->z[i+1],curr->v[i+1],curr->v[i],price[i+1],price[i],curr);
                curr->posterior[i+1] = curr->post_z[i+1] + curr->post_lik[i+1];
            }
            //cout<<curr->post_z[i]<<" "<<variance_gamma(curr->z[i],curr)<<endl;
            //cout<<curr->posterior[i]<<" "<<post_no_prior(1.0,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr)<<" "<<variance_gamma(curr->z[i],curr)+likelihood(1.0,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr)<<endl;
            acc->v[i] += 1.0;
        }
        else{
            count++;
            curr->v[i] = temp;
        }
    }
   return ;
}

void update_lat_z(para* curr,para* kernel,para* acc){
    for(int i=1;i<total;i++){
        double old_post = curr->posterior[i];
        double temp = curr->z[i];
        curr->z[i] = curr->z[i] + normal()*kernel->z[i];
        double new_post = variance_gamma(curr->z[i],curr) + likelihood(1.0,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr);
        double R = new_post - old_post;
        double u = uniform();
        //cout<<"R:"<<R<<endl;
        if(exp(R)>u){
            curr->sum_post = curr->sum_post + new_post - old_post;
            if(abs(posterior(1.0,curr)-curr->sum_post)>TOL)
                cout<<posterior(1.0,curr)<<" "<<curr->sum_post<<endl;
            curr->post_z[i] = variance_gamma(curr->z[i],curr);
            curr->post_lik[i] = likelihood(1.0,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr);
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
            acc->z[i] += 1.0;
        }
        else{
            curr->z[i] = temp;
            if(abs(posterior(1.0,curr)-curr->sum_post)>TOL)
                cout<<posterior(1.0,curr)<<" "<<curr->sum_post<<endl;
        }
    }
    
    return ;
}

void update_sv(para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double old_post = curr->sum_post+log(0.5*exp(curr->sigma_v)/sqrt(2*0.5*exp(curr->sigma_v)/exp(curr->k)*0.5*exp(curr->sigma_v)/exp(curr->v_p)-exp(curr->sigma_v)));
    double temp = curr->sigma_v;

    //curr->sigma_v = exp(log(curr->sigma_v/(curr->u3-curr->sigma_v)) + normal()*sd);
    //curr->u3 = sqrt(2*curr->k*curr->v_p);
    //curr->sigma_v = curr->u3+curr->u3/(curr->sigma_v-1);
    //curr->u3 = sqrt(2*curr->k*curr->v_p);
    curr->sigma_v = curr->sigma_v + normal()*sd;
    //cout<<"temp:"<<temp<<" "<<curr->sigma_v/(curr->u3-curr->sigma_v)<<" "<<curr->u3<<" "<<sd<<endl;

    double new_post = mcmc_posterior(curr,new_post_z,new_post_lik,new_posterior);
    //cout<<"new post:"<<new_post<<endl;
    double R = new_post - old_post+log(0.5*exp(curr->sigma_v)/sqrt(2*0.5*exp(curr->sigma_v)/exp(curr->k)*0.5*exp(curr->sigma_v)/exp(curr->v_p)-exp(curr->sigma_v)));

    double u = uniform();
    
   /* double y1 = curr->k;
    double y2 = curr->v_p;
    double y3 = curr->sigma_v;
    curr->k = 0.5*exp(y3)/exp(y1);
    curr->v_p = 0.5*exp(y3)/exp(y2);
    curr->sigma_v = sqrt(2*curr->k*curr->v_p-exp(y3));
    cout<<"k:"<<curr->k<<" vp:"<<curr->v_p<<" sv:"<<curr->sigma_v<<" "<<R<<endl;
    curr->k = y1;
    curr->v_p = y2;
    curr->sigma_v = y3;
    
    cout<<new_post<<" "<<log(0.5*exp(curr->sigma_v)/sqrt(2*0.5*exp(curr->sigma_v)/exp(curr->k)*0.5*exp(curr->sigma_v)/exp(curr->v_p)-exp(curr->sigma_v)))<<endl;*/
    if(exp(R)>u){
        //cout<<"accepted sv"<<endl;
        curr->sum_post = new_post;
        delete curr->post_z;
        delete curr->post_lik;
        delete curr->posterior;
        curr->post_z = new_post_z;
        curr->post_lik = new_post_lik;
        curr->posterior = new_posterior;
        //curr->u1 = 0.5*curr->sigma_v*curr->sigma_v/curr->v_p;
        //curr->u2 = 0.5*curr->sigma_v*curr->sigma_v/curr->k;
        //cout<<"R:"<<R<<" :"<<old_post<<" new_post"<<new_post<<endl;
        //cout<<curr->u3-curr->sigma_v<<endl;
        acc->sigma_v += 1.0;
    }
    else{
        curr->sigma_v = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
        //cout<<curr->u3-curr->sigma_v<<endl;
    }

    
    return ;
}

void update_vp(para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    double old_post = curr->sum_post+exp(curr->v_p);
    double temp = curr->v_p;
    curr->v_p = curr->v_p + normal()*sd;

    double new_post = mcmc_posterior(curr,new_post_z,new_post_lik,new_posterior);
    double R = new_post - old_post +exp(curr->v_p);
    //cout<<"new post:"<<new_post<<endl;
    double u = uniform();
    //cout<<"R:"<<R<<" log u:"<<u<<endl;
    if(exp(R)>u){
        //cout<<"acc v"<<endl;
        curr->sum_post = new_post;
        delete curr->post_z;
        delete curr->post_lik;
        delete curr->posterior;
        curr->post_z = new_post_z;
        curr->post_lik = new_post_lik;
        curr->posterior = new_posterior;
        //curr->u1 = 0.5*curr->sigma_v*curr->sigma_v/curr->v_p;
        //curr->u3 = sqrt(2*curr->v_p*curr->k);
        //cout<<"R:"<<R<<" :"<<old_post<<" new_post"<<new_post<<endl;
        /*
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }*/
        acc->v_p += 1.0;
    }
    else{
        curr->v_p = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }

    
    return ;
}

void update_k(para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    
    double old_post = curr->sum_post+exp(curr->k);
    double temp = curr->k;
    curr->k = curr->k + normal()*sd;
    //curr->u1 = pow(curr->sigma_v,2)/(2*curr->v_p);

    double new_post = mcmc_posterior(curr,new_post_z,new_post_lik,new_posterior);
    double R = new_post - old_post+exp(curr->k);
    //cout<<"new post:"<<new_post<<endl;
    double u = uniform();
    //cout<<"R:"<<R<<" log u:"<<u<<endl;
    if(exp(R)>u){
        //cout<<"acc k"<<endl;
        curr->sum_post = new_post;
        delete curr->post_z;
        delete curr->post_lik;
        delete curr->posterior;
        curr->post_z = new_post_z;
        curr->post_lik = new_post_lik;
        curr->posterior = new_posterior;
        //curr->u2 = 0.5*curr->sigma_v*curr->sigma_v/curr->k;
        //curr->u3 = sqrt(2*curr->v_p*curr->k);
        // cout<<"R:"<<R<<" :"<<old_post<<" new_post"<<new_post<<endl;
        /*for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }*/
        acc->k += 1.0;
    }
    else{
        curr->k = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }
    
    return ;
}

void update_rho(para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    
    double old_post = curr->sum_post+log(curr->rho+1);
    double temp = curr->rho;
    curr->rho = exp(log((1+curr->rho)/(1-curr->rho)) + normal()*sd);
    curr->rho = (curr->rho-1)/(curr->rho+1);
    double new_post = mcmc_posterior(curr,new_post_z,new_post_lik,new_posterior);
    double R = new_post - old_post+log(curr->rho+1);
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(exp(R)>u){
        curr->sum_post = new_post;
        delete curr->post_z;
        delete curr->post_lik;
        delete curr->posterior;
        curr->post_z = new_post_z;
        curr->post_lik = new_post_lik;
        curr->posterior = new_posterior;
        //cout<<"R:"<<R<<" :"<<old_post<<" new_post"<<new_post<<endl;
        /*
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }*/
        acc->rho += 1.0;
    }
    else{
        curr->rho = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }

    
    return ;
}

void update_sj(para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    
    double old_post = curr->sum_post+curr->sigma_j;//curr->sum_post;
    double temp = curr->sigma_j;
    curr->sigma_j = exp(log(curr->sigma_j) + normal()*sd);
    double new_post = mcmc_posterior(curr,new_post_z,new_post_lik,new_posterior);
    double R = new_post - old_post+curr->sigma_j;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(exp(R)>u){
        curr->sum_post = new_post;
        delete curr->post_z;
        delete curr->post_lik;
        delete curr->posterior;
        curr->post_z = new_post_z;
        curr->post_lik = new_post_lik;
        curr->posterior = new_posterior;
        //cout<<"R:"<<R<<" :"<<old_post<<" new_post"<<new_post<<endl;
        /*
        for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }*/
        acc->sigma_j += 1.0;
    }
    else{
        curr->sigma_j = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
    }
    
    return ;
}

void update_gam(para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    
    double old_post = curr->sum_post;
    double temp = curr->gam;
    curr->gam = curr->gam + normal()*sd;
    double new_post = mcmc_posterior(curr,new_post_z,new_post_lik,new_posterior);
    double R = new_post - old_post;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(exp(R)>u){
        curr->sum_post = new_post;
        delete curr->post_z;
        delete curr->post_lik;
        delete curr->posterior;
        curr->post_z = new_post_z;
        curr->post_lik = new_post_lik;
        curr->posterior = new_posterior;
       // cout<<"R:"<<R<<" :"<<old_post<<" new_post"<<new_post<<endl;
        /*
        for(int i=0;i<total;i++){
        curr->post_z[i] = new_post_z[i];
        curr->post_lik[i] = new_post_lik[i];
        curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }*/
            acc->gam += 1.0;
    }
    else{
        curr->gam = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;

    }

    
    return ;
}

void update_mu(para* curr,double sd,para* acc){
    double* new_post_lik = new double[total];
    double* new_post_z = new double[total];
    double* new_posterior = new double[total];
    
    double old_post = curr->sum_post;
    //cout<<"old_post:"<<old_post<<" ";
    double temp = curr->mu;
    curr->mu = curr->mu + normal()*sd;
    double new_post = mcmc_posterior(curr,new_post_z,new_post_lik,new_posterior);
    double R = new_post - old_post;
    double u = uniform();
    //cout<<"R:"<<R<<endl;
    if(exp(R)>u){
        //cout<<"R:"<<R<<" :"<<old_post<<" new_post"<<new_post<<endl;
        //cout<<curr->sum_post<<" ";
        curr->sum_post = new_post;
        //cout<<curr->sum_post<<" "<<posterior(1.0,curr)<<endl;
        int i = 10;
        double temp = curr->posterior[i];
        delete curr->post_z;
        delete curr->post_lik;
        delete curr->posterior;
        curr->post_z = new_post_z;
        curr->post_lik = new_post_lik;
        curr->posterior = new_posterior;
        //cout<<"prev:"<<temp<<" now:"<<curr->post_lik[i]<<" correct:"<<likelihood(1.0,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr)<<endl;
        //cout<<"prev:"<<temp<<" now:"<<curr->posterior[i]<<" correct:"<<variance_gamma(curr->z[i],curr)+likelihood(1.0,curr->z[i],curr->v[i],curr->v[i-1],price[i],price[i-1],curr)<<endl;
        /*for(int i=0;i<total;i++){
            curr->post_z[i] = new_post_z[i];
            curr->post_lik[i] = new_post_lik[i];
            curr->posterior[i] = curr->post_z[i]+curr->post_lik[i];
        }*/
        acc->mu += 1.0;

    }
    else{
        curr->mu = temp;
        delete new_post_lik;
        delete new_post_z;
        delete new_posterior;
        //cout<<curr->sum_post<<" correct:"<<posterior(1.0,curr)<<" rej"<<endl;
        //cout<<"mu:"<<curr->mu<<" old_post:"<<old_post<<" new_post"<<new_post<<endl;
    }

    
    return ;
}




