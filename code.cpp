#include <iostream>
#include <math.h>
#include "include.h"
#include <random>
#include <vector>
#include <chrono>
#include <time.h>
#include <complex>
#include <fstream>
#include <string>


double TOL = 0.000001;
double ep = pow(10,-300);
using namespace std;

para* head;
para* kernel;
para* acc;

int no_particles=2000;
double dmax = 2;
double delta=1;
vector<double> price;
int total;

double zeta = 0.00001;
double prev_zeta = 0.000;
double ESS_k;
double threshold = 1000.0;
int file=1;
string data_folder = "data/";
string output_folder = "output/";


//The codes in this file and SP500 are the same
//Difference is this one is for retrieving the posterior parameters from the simulated SVVG Model
//the code in SP500 is for retrieving the posterior parameters for the actual SP500 data

int main(){
    while(file != 101){
    srand(time(0));
    init();
    init_post();
    update_norm_weights(zeta,prev_zeta);
    ESS_k = ESS_0();
    if(ESS_k<threshold){
        cout<<"resample!"<<endl;
        resample();
        ESS_k = no_particles;
    }
    
    set_kernel();
    acc_init();
    int count = 1;
    int flag = 1;
    while(zeta<1.000001){
        update_para(zeta);
        
        cout<<count<<","<<zeta<<","<<ESS_k<<",";
        print_acc(acc);
    
        prev_zeta = zeta;
        zeta = find_new_zeta(prev_zeta,1.0,prev_zeta,ESS_k);
        if(!flag)
            break;
        if(zeta-0.9999>TOL)
            flag = 0;
        ESS_k = ESS(zeta,prev_zeta);
        update_norm_weights(zeta,prev_zeta);
        if(ESS_k<threshold){
            resample();
            ESS_k = no_particles;
        }
        
        set_kernel();
        adapt_kernel();
        acc_init();
        count++;
    }
        print();
        file++;
        zeta = 0.00001;
        prev_zeta = 0.0;
        reset();
        cout<<"reset"<<endl;
    }
	return 0;
}
