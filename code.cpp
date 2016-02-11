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

int no_particles=100;
double dmax = 2;
double delta=1;
vector<double> price;
int total;
int flag = 1;

double zeta = 0.005;
double prev_zeta = 0.000;
double ESS_k;
double threshold = 50.0;

int main(){

    srand(10);
    init();
    init_post();
    //print(head);
    update_norm_weights(zeta,prev_zeta);
    //print(head);
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
    //cout<<"kernel: ";
    //print(kernel);
    //cout<<"real:";
    //print(head);

    int i = 1;
    int flag = 1;
    while(zeta <=1.00000001){
        //cout<<"kernel: ";
        //print(kernel);
        acc = acc_init();
        cout<<"accept rate init"<<endl;
        update_para(kernel,acc);   //untempered_lik updated already
        cout<<"update done"<<endl;
        kernel = reset_kernel(kernel); // inclusive of delete
        cout<<"delete done"<<endl;
        adapt_kernel(kernel,acc);
        cout<<"adapt done"<<endl;
        print_acc(acc);
        acc_end(acc);
        
        prev_zeta = zeta;
        cout<<"start new zeta"<<endl;
        if(zeta>0.99){
            zeta = 1.0;
            if(flag)
                flag = 0;
            else{
                cout<<"done"<<endl;
            }
        }
        else{
            zeta = find_new_zeta(prev_zeta,1.0,prev_zeta,ESS_k);
        }
        cout<<"zeta found"<<endl;
        update_norm_weights(zeta,prev_zeta);
        cout<<"weights/norm weights updated"<<endl;
        ESS_k = ESS(zeta,prev_zeta);
        cout<<"ESS done"<<endl;
        if(ESS_k<threshold){
            cout<<"RESAMPLE!"<<endl;
            resample();
            ESS_k = no_particles;
        }
        

         
        //zeta = zeta + 1.0;
        cout<<i<<" "<<zeta<<endl;
        print(head);
        i++;
        if(!flag)
            break;
    }
    get_val(head);
    
	return 0;
}
