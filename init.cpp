#include "include.h"
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <boost/math/distributions.hpp>
#include "rnglib.h"
#include "ranlib.h"

using namespace std;

void print(){
    para* curr = head;
    while(curr != NULL){
        cout<<"mu: "<<curr->mu<<" gam: "<<curr->gam<<" s_j: "<<curr->sigma_j<<" s_v: "<<curr->sigma_v<<" rho: "<<curr->rho<<" k: "<<curr->k<<" v: "<<curr->v_p<<endl;
        curr = curr->next;
    }
    return ;
}

double normal(){
    double s1 = uniform();
    double s2 = uniform();
    return sqrt(-2*log(s1))*cos(2*M_PI*s2);
}


double uniform(){
    return rand()%10000/10000.0;
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


void initialize(){
    max_v = 10000;
    min_v = 0;
    min_w = 0;
    max_w = 10000;
    
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
    default_random_engine gen(rand());
    gamma_distribution<double> dist(0.142857, 1.0);
    gamma_distribution<double> dist1(1.0,1.0);
    for(int i=0;i<no_particles;i++){
        temp = new para;
        temp->mu = normal();
        temp->gam = normal();
        temp->sigma_j = rgamma(1);
        temp->sigma_v = rgamma(1);
        temp->v_p = rgamma(1);
        double a = dist(gen);
        temp->rho = 2*(ep+a/(a+dist1(gen)))-1;
        temp->k = uniform()*(0.5*pow(temp->sigma_v,2)/temp->v_p)+pow(temp->sigma_v,2)/temp->v_p*0.5;
        temp->weight = 1/no_particles;
        temp->norm_weight = 1/no_particles;
        temp->z = new double[total];
        temp->w = new double[total];
        temp->v_star = new double[total];
        temp->v = new double[total];
        double d = 4*temp->k*temp->v_p/pow(temp->sigma_v,2);
        double lambda;
        for(int j=0;j<total;j++){
            
            if(j==0){
                temp->z[j] = 1;
                temp->w[j] = 1;
                temp->v_star[j] = 1;
                temp->v[j] = 1;
            }
            else{
                lambda = 4*temp->k*exp(-temp->k*delta)*temp->v[j-1]/(pow(temp->sigma_v,2)*(1-exp(-temp->k*delta)));
                boost::math::non_central_chi_squared_distribution<double> chisq(d,lambda);
                //cout<<ranlib::genchi (d)<<endl;
                temp->v[j] = 1;
                temp->z[j] = 1;
                temp->w[j] = 1;
                temp->v_star[j] = 1;
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

