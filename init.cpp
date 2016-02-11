#include "include.h"
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <boost/math/distributions.hpp>


using namespace std;

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
    
    double t_w = 0.0;
    while(p != NULL){
        t_w += p->weight;
        curr->mu += p->mu*p->norm_weight;
        curr->gam += p->gam*p->norm_weight;
        curr->sigma_j += p->sigma_j*p->norm_weight;
        curr->rho += p->rho*p->norm_weight;
        curr->k += p->k*p->norm_weight;
        curr->v_p += p->v_p*p->norm_weight;
        curr->sigma_v += p->sigma_v*p->norm_weight;
        p = p->next;
    }
    cout<<"mu:"<<curr->mu<<" gamma:"<<curr->gam<<" sigma_j:"<<curr->sigma_j<<" rho:"<<curr->rho<<" k:"<<curr->k<<" v:"<<curr->v_p<<" sigma_v:"<<curr->sigma_v<<endl;
}

void print_acc(para* curr){
    double min1 = 10;double min2 = 10;
    double max1 = -1;double max2 = -1;
    for(int i=0;i<total;i++){
        if(min1>curr->z[i])
            min1 = curr->z[i];
        if(max1<curr->z[i])
            max1 = curr->z[i];
        
        if(min2>curr->v[i])
            min2 = curr->v[i];
        if(max1<curr->v[i])
            max2 = curr->v[i];
        
    }

    cout<<" mu:"<<curr->mu<<" gam:"<<curr->gam<<" s_j:"<<curr->sigma_j<<" s_v:"<<curr->sigma_v<<" rho:"<<curr->rho<<" k:"<<curr->k<<" v:"<<curr->v_p<<" min z:"<<min1<<" max z:"<<max1<<" min v:"<<min2<<" max v:"<<max2<<endl;
    return ;
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
void print(para* p){
    para* curr = p;int i =0;
    while(curr != NULL){
        cout<<" sum temp post:"<<curr->sum_post<<" mu: "<<curr->mu<<" gam: "<<curr->gam<<" s_j: "<<curr->sigma_j<<" s_v: "<<curr->sigma_v<<" rho: "<<curr->rho<<" k: "<<curr->k<<" v: "<<curr->v_p<<" w: "<<curr->weight<<" norm_w: "<<curr->norm_weight<<" cum_w: "<<curr->cum_norm_weight<<" untempered_lik: "<<untempered_lik[i]<<endl;
        //cout<<"s_v-u3: "<<curr->u3-curr->sigma_v<<" k-u1: "<<curr->k-curr->u1<<" v-u2: "<<curr->v_p-curr->u2<<endl;
        //cout<<" mu: "<<curr->mu<<" gam: "<<curr->gam<<" s_j: "<<curr->sigma_j<<" s_v: "<<curr->sigma_v<<" rho: "<<curr->rho<<" k: "<<curr->k<<" v: "<<curr->v_p<<endl;
        //if(curr->u3 !=0)
        //    cout<<"u3 - sv:"<<curr->u3-curr->sigma_v<<endl;
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

        temp->v_p = gengam(1.0,1.0);
        
        temp->sigma_j = 1.0/gengam(2.5,5.0);//1.0/gengam(1.0,0.5);
        //temp->sigma_v = gengam(1.0,1.0);
        double w_v = 1/gengam(1.0,0.5);
        double phi_v = normal()*0.5/w_v;
        temp->sigma_v = sqrt(phi_v*phi_v+w_v);
        temp->rho =phi_v/temp->sigma_v;//2*genbet(0.14285714,1.0)-1+0.000001;

        
        temp->k = (gengam(1.0,0.5)+1.01)*pow(temp->sigma_v,2)/temp->v_p*0.5;
        temp->weight = 0;
        temp->norm_weight = 1;
        temp->cum_norm_weight = 0;
        temp->z = new double[total];
        temp->v = new double[total];
        temp->posterior = new double[total];
        temp->post_lik = new double[total];
        temp->post_z = new double[total];
        temp->u1 = pow(temp->sigma_v,2)/(2*temp->v_p);
        temp->u2 = pow(temp->sigma_v,2)/(2*temp->k);
        temp->u3 = 2*temp->k*temp->v_p;

        temp->k = log(temp->k - temp->u1);
        temp->v_p = log(temp->v_p - temp->u2);
        temp->sigma_v = log(temp->u3 - temp->sigma_v*temp->sigma_v);
        //if(2*temp->k*temp->v_p-temp->sigma_v*temp->sigma_v<0)
          //  cout<<"mu:"<<temp->mu<<" gam:"<<temp->gam<<" vp:"<<temp->v_p<<" sigj:"<<temp->sigma_j<<" sigv:"<<temp->sigma_v<<" rho:"<<temp->rho<<" k:"<<temp->k<<endl;
        double d = 4*temp->k*temp->v_p/pow(temp->sigma_v,2);
        double norm_y;
        for(int j=0;j<total;j++){
            
            if(j==0){
                temp->z[j] = 1;
                temp->v[j] = 1;
            }
            else{
                double temp_gamma = gengam(1.0,1.0/delta);
                temp->z[j] = temp->gam*temp_gamma+temp->sigma_j*sqrt(temp_gamma)*normal();
                //norm_y = (price[j]-price[j-1]-temp->mu*delta-temp->z[j])/sqrt(temp->v[j-1]*delta);
                temp->v[j] = uniform();//temp->v[j-1]*exp(temp->k*(temp->v_p-temp->v[j-1])*delta+temp->sigma_v*sqrt(temp->v[j-1])*uniform());//gennch(d,lambda*temp->v[j-1]*exp(-temp->k*delta))/lambda;

            }
            
        }
        //if((0.5*4*temp->k*temp->v_p/pow(temp->sigma_v,2))>dmax) cout<<"out"<<endl;


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

void init_post(){
    double sum,like;int i=0;
    para* p = head;
    while(p != NULL){
        sum = 0;
        like = 0;
        for(int j=0;j<total;j++){
        
            if(j==0){
                p->posterior[j] = 0;
                p->post_z[j] = 0;
                p->post_lik[j] = 0;
            }
            else{
                p->post_z[j] = variance_gamma(p->z[j],p);
                p->post_lik[j] = likelihood(1.0,p->z[j],p->v[j],p->v[j-1],price[j],price[j-1],p);
                p->posterior[j] = p->post_z[j]+p->post_lik[j];
                sum += p->posterior[j];
                like += p->post_lik[j];
            }
        }
        p->sum_post = sum+prior(p);
        untempered_lik[i] = like;
        i++;
        p = p->next;
    }
    return ;
}
