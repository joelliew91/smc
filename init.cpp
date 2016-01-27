#include "include.h"
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <boost/math/distributions.hpp>


using namespace std;
double max(){
    double v = -100000;
    for(int i=0;i<price.size();i++)
        if(v<price[i])
            v = price[i];
    return v;
}
void print(para* p){
    para* curr = p;int i =0;
    while(curr != NULL){
        cout<<" sum temp post:"<<curr->sum_post<<" posterior:"<<posterior(1.0,curr)<<" mu: "<<curr->mu<<" gam: "<<curr->gam<<" s_j: "<<curr->sigma_j<<" s_v: "<<curr->sigma_v<<" rho: "<<curr->rho<<" k: "<<curr->k<<" v: "<<curr->v_p<<" w: "<<curr->weight<<" norm_w: "<<curr->norm_weight<<" cum_w: "<<curr->cum_norm_weight<<" untempered_lik: "<<untempered_lik[i]<<endl;
        //cout<<"s_v-u3: "<<curr->u3-curr->sigma_v<<" k-u1: "<<curr->k-curr->u1<<" v-u2: "<<curr->v_p-curr->u2<<endl;
        curr = curr->next;
        i++;
    }
    return ;
}

double normal(){
    double s1 = uniform();
    double s2 = uniform();
    if(s1<0.0000001) s1 += 0.0000001;
    return sqrt(-2*log(s1))*cos(2*M_PI*s2);
}


double uniform(){
    return rand()%10000000/pow(10.0,7)+ep;
}

double rgamma(int alpha){
    double sum=0;
    double s;
    for(int i=0;i<alpha;i++){
        s = uniform();
        sum += -log(s);
    }
    return sum;
}


void init(){

    ifstream myfile("data.csv");
    int count = 0;
    string line;
    while(myfile.good()){
        getline(myfile,line);
        int pos = line.find(",");
        if(pos==3){
            price.push_back(stod(line.substr(pos+1)));
        }
        
    }
    myfile.close();
    
    total = 1*price.size(); //check delta*price.size();
    para* temp;
    para* curr;
    for(int i=0;i<no_particles;i++){
        temp = new para;
        temp->norm_weight = 1;
        temp->mu = normal();
        temp->gam = normal();
        temp->sigma_j = 1.0/gengam(1.0,0.5);
        temp->sigma_v = gengam(1.0,1.0);
        temp->v_p = gengam(1.0,1.0);
        temp->rho = 2*genbet(0.14285714,1.0)-1+0.000001;
        temp->k = (uniform()+1)*pow(temp->sigma_v,2)/temp->v_p*0.5;
        temp->weight = 0;
        temp->norm_weight = 0;
        temp->cum_norm_weight = 0;
        temp->z = new double[total];
        temp->v = new double[total];
        temp->posterior = new double[total];
        temp->post_lik = new double[total];
        temp->post_z = new double[total];
        temp->u1 = pow(temp->sigma_v,2)/(2*temp->v_p);
        temp->u2 = pow(temp->sigma_v,2)/(2*temp->k);
        temp->u3 = sqrt(2*temp->k*temp->v_p);
        //cout<<2*temp->k*temp->v_p/pow(temp->sigma_v,2)<<" "<<temp->k - temp->u1<<" "<<temp->v_p - temp->u2<<" "<<temp->u3-temp->sigma_v<<endl;
        double d = 4*temp->k*temp->v_p/pow(temp->sigma_v,2);
        double lambda;
        for(int j=0;j<total;j++){
            
            if(j==0){
                temp->z[j] = 1;
                temp->v[j] = 1;
            }
            else{
                int flag = 1;
                lambda = 4*temp->k/(pow(temp->sigma_v,2)*(1-exp(-temp->k*delta)));
                temp->v[j] = uniform();//gennch(d,lambda*temp->v[j-1]*exp(-temp->k*delta))/lambda;
                //cout<<temp->v[j]<<endl;
                double temp_gamma = gengam(1.0,1.0/delta);
                temp->z[j] = temp->gam*temp_gamma+temp->sigma_j*sqrt(temp_gamma)*normal();
            }
            
        }
        if((0.5*4*temp->k*temp->v_p/pow(temp->sigma_v,2))>dmax) cout<<"out"<<endl;


        if(head==NULL){
            head = temp;
            curr = head;
        }
        else{
            curr->next = temp;
            curr = temp;
        }
    }
    curr->next = NULL;
    return ;
}

void init_post(double zeta){
    double sum;
    para* p = head;
    while(p != NULL){
        sum = 0;
        for(int j=0;j<total;j++){
        
            if(j==0){
                p->posterior[j] = 0;
                p->post_z[j] = 0;
                p->post_lik[j] = 0;
            }
            else{
                p->post_z[j] = variance_gamma(p->z[j],p);
                p->post_lik[j] = likelihood(zeta,p->z[j],p->v[j],p->v[j-1],price[j],price[j-1],p);
                p->posterior[j] = p->post_z[j]+p->post_lik[j];
                sum += p->posterior[j];
            }
        }
        p->sum_post = sum+prior(p);
        p->weight = exp(p->sum_post - posterior(0.0,p));
        if(p->weight<ADAPT_WEIGHT) p->weight = 0.0;
        p = p->next;
    }
    return ;
}
