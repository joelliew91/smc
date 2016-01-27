#include <iostream>
#include <math.h>
#include "include.h"
#include <random>
#include <vector>
#include <chrono>
#include <time.h>
#include <complex>

double ADAPT_WEIGHT = pow(10,-50);
double TOL = 0.000001;
double ep = pow(10,-300);
double untempered_lik[2000];
using namespace std;

para* head = NULL;

int no_particles=2000;
double dmax = 2;
double delta=1;
vector<double> price;
int total;
int flag = 1;

double zeta = 0.005;
double prev_zeta = 0.000;
double ESS_k;
double threshold = 1000.0;

int main(){
    srand(10);
    init();
    init_post(zeta);
    update_norm_weights(head);
    
    para* p = head;
    para* kernel;
    para* acc;
    
    ESS_k = ESS_0();
    
    if(ESS_k<threshold){
        cout<<"RESAMPLE!"<<endl;
        resample();
        ESS_k = no_particles;
    }
    
    kernel = set_kernel();
    int i = 1;
    update_untempered_lik();
    prev_zeta = zeta;
    //print_part_post();
    //print(head);
    while(zeta <0.0051){
        for(int j = 0;j<10;j++){
            acc = acc_init();
            cout<<"accept rate init"<<endl;
            update_para(zeta,kernel,acc);
            cout<<"update done"<<endl;
            kernel = reset_kernel(kernel); // inclusive of delete
            cout<<"delete done"<<endl;
            adapt_kernel(kernel,acc);
            cout<<"adapt done"<<endl;
            acc_end(acc);
            update_untempered_lik();
            cout<<"update untempered done"<<endl;
            cout<<j<<endl;
        }
        prev_zeta = zeta;
        /*cout<<"start new zeta"<<endl;
        zeta = find_new_zeta(prev_zeta,1.0,prev_zeta,ESS_k);
        cout<<"zeta found"<<endl;
        ESS_k = ESS(zeta,prev_zeta);
        cout<<"ess found"<<endl;
        
        if(ESS_k<threshold){
            cout<<"RESAMPLE!"<<endl;
            resample();
            ESS_k = no_particles;
        }
        */

         

        cout<<i<<" "<<zeta<<endl;
        i++;
        zeta += 1;
    }
    update_all();
    update_untempered_lik();
    //lik(1.0,head);
    print(head);
    /*
    para* acc;
    para* kernel = set_kernel();
    for(int i =0;i<20;i++){
        cout<<i<<endl;
        acc = acc_init();
        update_para(zeta,kernel,acc);
        kernel = reset_kernel(kernel);
        adapt_kernel(kernel,acc);
        acc_end(acc);
    }
    print(head);*/
    //resample_tester();
    //resample();
    /*ESS_k = ESS(zeta,prev_zeta);
    update_norm_weights();
    if(ESS_k<threshold)
        resample();
     
    para* temp = head;
    while(zeta<1.0){
        temp = head;
        while(temp != NULL){
            //update para
            temp = temp->next;
        }
        ESS_k = ESS(zeta,prev_zeta);
        update_norm_weights();
        if(ESS_k<1000)
            resample();

        prev_zeta = zeta;
        zeta = find_new_zeta(prev_zeta,1.0,prev_zeta,ESS_k);
        //set kernel
        cout<<zeta<<endl;
    }*/
    
	return 0;
}
