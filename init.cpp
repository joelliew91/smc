#include "include.h"
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <boost/math/distributions.hpp>


using namespace std;
/*
void get_val(para* a){
    para* p = a;
    para* curr = new para;
    curr->mu = 0;
    curr->gam = 0;
    curr->sigma_j = 0;
    curr->rho = 0;
    curr->k = 0;
    curr->v_p = 0;
    curr->sigma_v = 0;
    double y1,y2,y3;
    double t_w = 0.0;
    while(p != NULL){
        y1 = p->k;
        y2 = p->v_p;
        y3 = p->sigma_v;
        p->k = 0.5*exp(y3)/exp(y1);
        p->v_p = 0.5*exp(y3)/exp(y2);
        p->sigma_v = sqrt(2*p->k*p->v_p-exp(y3));
        
        t_w += p->weight;
        curr->mu += p->mu*p->norm_weight;
        curr->gam += p->gam*p->norm_weight;
        curr->sigma_j += p->sigma_j*p->norm_weight;
        curr->rho += p->rho*p->norm_weight;
        curr->k += p->k*p->norm_weight;
        curr->v_p += p->v_p*p->norm_weight;
        curr->sigma_v += p->sigma_v*p->norm_weight;
        p->k = y1;
        p->v_p = y2;
        p->sigma_v = y3;
        p = p->next;
    }
    cout<<"mu:"<<curr->mu<<" gamma:"<<curr->gam<<" sigma_j:"<<curr->sigma_j<<" rho:"<<curr->rho<<" k:"<<curr->k<<" v:"<<curr->v_p<<" sigma_v:"<<curr->sigma_v<<endl;
}

double min(){
    double v = 100000;
    for(int i=1;i<no_particles;i++)
        if(v>untempered_lik[i])
            v = untempered_lik[i];
    return v;
}
double max(){
    double v = -pow(10,100);
    for(int i=1;i<no_particles;i++)
        if(v<untempered_lik[i])
            v = untempered_lik[i];
    return v;
}
void print_extra(){
    para* p = head;int i=0;
    while(p!=NULL){
        cout<<"sum_post:"<<p->sum_post<<" untempered_lik:"<<untempered_lik[i]<<" W:"<<p->weight<<" Norm_W:"<<p->norm_weight<<" cum_norm_w:"<<p->cum_norm_weight<<endl;
        p = p->next;
        i++;
    }
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
*/

void print_acc(para* curr){
    
    cout<<curr->mu/no_particles<<","<<curr->gam/no_particles<<","<<curr->sigma_j/no_particles<<","<<curr->sigma_v/no_particles<<","<<curr->rho/no_particles<<","<<curr->k/no_particles<<","<<curr->v_p/no_particles<<","<<curr->z[0]/(no_particles*total)<<","<<curr->v[0]/(no_particles*total)<<endl;
    return ;
}

void init_post(){
    double sum,like;
    para* p = head;
    for(int i =0;i<no_particles;i++){
        sum = 0;
        like = 0;
        for(int j=0;j<total;j++){
            
            if(j==0){
                p[i].posterior[j] = 0;
                p[i].post_z[j] = 0;
                p[i].post_lik[j] = 0;
            }
            else{
                p[i].post_z[j] = variance_gamma(p[i].z[j],p[i]);
                p[i].post_lik[j] = likelihood(1.0,p[i].z[j],p[i].v[j],p[i].v[j-1],price[j],price[j-1],p[i]);
                p[i].posterior[j] = p[i].post_z[j]+p[i].post_lik[j];
                sum += p[i].posterior[j];
                like += p[i].post_lik[j];
            }
        }
        p[i].sum_post = sum+prior(p[i]);
        p[i].sum_lik = like;
    }
    return ;
}

void print(){
    ofstream myfile;
    if(file==0)
        myfile.open("output.csv");
    else{
        string a = "output";
        string b = ".csv";
        string name = a+to_string(file)+b;
        myfile.open(name);
    }
    para* p = head;
    double y1,y2,y3,k,v,sv;
    myfile<<"mu,gam,sigma_j,sigma_v,rho,k,v,w,norm_w,sum"<<endl;
    for(int i =0;i<no_particles;i++){
        y1 = p[i].k;
        y2 = p[i].v_p;
        y3 = p[i].sigma_v;
        k = 0.5*exp(y3)/exp(y1);
        v = 0.5*exp(y3)/exp(y2);
        sv = sqrt(2*k*v-exp(y3));
        myfile<<p[i].mu<<","<<p[i].gam<<","<<p[i].sigma_j<<","<<sv<<","<<p[i].rho<<","<<k<<","<<v<<","<<p[i].weight<<","<<p[i].norm_weight<<","<<p[i].sum_post<<endl;
        //cout<<" sum temp post:"<<p[i].sum_post<<" w: "<<p[i].weight<<" norm_w: "<<p[i].norm_weight<<" cum_w: "<<p[i].cum_norm_weight<<" untempered_lik: "<<p[i].sum_lik<<endl;
        //cout<<"s_v-u3: "<<curr->u3-curr->sigma_v<<" k-u1: "<<curr->k-curr->u1<<" v-u2: "<<curr->v_p-curr->u2<<endl;
        //cout<<" mu: "<<curr->mu<<" gam: "<<curr->gam<<" s_j: "<<curr->sigma_j<<" s_v: "<<curr->sigma_v<<" rho: "<<curr->rho<<" k: "<<curr->k<<" v: "<<curr->v_p<<endl;
    }
    myfile.close();
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

void init(){
    kernel = new para;
    kernel->v = new double[1];
    kernel->z = new double[1];
    
    acc = new para;
    acc->v = new double[1];
    acc->z = new double[1];
    
    head = new para[no_particles];
    ifstream myfile;
    if(file==0)
        myfile.open("data.csv");
    else{
        string a = "data";
        string b = ".csv";
        string name = a+to_string(file)+b;
        myfile.open(name);
    }
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
    double w_v,phi_v,u1,u2,u3;
    total = 1*price.size(); //check delta*price.size();
    for(int i=0;i<no_particles;i++){
        
        head[i].mu = normal();
        head[i].gam = normal();
        
        head[i].v_p = (uniform()+0.001)*0.3;
        head[i].k = (uniform()+0.001)*0.3;
        head[i].sigma_j = 1.0/gengam(2.5,5.0);//1.0/gengam(1.0,0.5);
        //temp->sigma_v = gengam(1.0,1.0);
        w_v = 1/gengam(1.0,0.5);
        phi_v = normal()*0.5/w_v;
        head[i].sigma_v = sqrt(phi_v*phi_v+w_v);
        head[i].rho =phi_v/head[i].sigma_v;//2*genbet(0.14285714,1.0)-1+0.000001;
        head[i].sigma_v = sqrt((uniform()+0.001)*2*head[i].k*head[i].v_p/3);
        
        //temp->k = (gengam(1.0,0.5)+1.01)*pow(temp->sigma_v,2)/temp->v_p*0.5;
        //cout<<temp->sigma_v<<" "<<temp->v_p<<" "<<temp->k<<endl;
        head[i].sum_lik = 0;
        head[i].sum_post=0;
        head[i].weight = 1;
        head[i].norm_weight = 1;
        head[i].cum_norm_weight = i*head[i].norm_weight;
        head[i].z = new double[total];
        head[i].v = new double[total];
        head[i].posterior = new double[total];
        head[i].post_lik = new double[total];
        head[i].post_z = new double[total];
        u1 = pow(head[i].sigma_v,2)/(2*head[i].v_p);
        u2 = pow(head[i].sigma_v,2)/(2*head[i].k);
        u3 = 2*head[i].k*head[i].v_p;
        
        head[i].k = log(head[i].k - u1);
        head[i].v_p = log(head[i].v_p - u2);
        head[i].sigma_v = log(u3 - head[i].sigma_v*head[i].sigma_v);
        //cout<<"mu:"<<head[i].mu<<" gam:"<<head[i].gam<<" vp:"<<head[i].v_p<<" sigj:"<<head[i].sigma_j<<" sigv:"<<head[i].sigma_v<<" rho:"<<head[i].rho<<" k:"<<head[i].k<<endl;
        for(int j=0;j<total;j++){
            if(j==0){
                //double temp_gamma = gengam(1.0,1.0);
                head[i].z[j] = 1;
                head[i].v[j] = 0.1;
            }
            else{
                double temp_gamma = gengam(1.0,1.0);
                head[i].z[j] = head[i].gam*temp_gamma+head[i].sigma_j*sqrt(temp_gamma)*normal();
                head[i].v[j] = (uniform()+0.001)*0.2;
                
            }
            
        }
        
    }
    
    return ;
}
