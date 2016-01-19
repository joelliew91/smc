#include <iostream>
#include <math.h>
#include "include.h"
#include <random>
#include <vector>
#include <chrono>
#include <complex_bessel.h>
#include <time.h>
#include <complex>


double TOL = 0.000001;
double ep = pow(10,-300);

using namespace std;

para* head = NULL;
double max_v;
double max_w;
double min_v;
double min_w;
int no_particles=1000;
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
    para* temp = head;
    while(temp!=NULL){
        cout<<posterior(0.005,temp)<<endl;;
        temp = temp->next;
    }
    //para* ker = set_kernel();
    //print(head);
    /*for(int i =0;i<20;i++){
        update_para(zeta,ker);
        cout<<i<<endl;
    }*/
    //resample_tester();
    //resample();
    /*ESS_k = ESS(zeta,prev_zeta);
    update_norm_weights();
    if(ESS_k<threshold)
        resample();
     
        //complex<double> va(700.0,2.0);
        //cout<<sp_bessel::besselI(0.95,va)<<endl;
        //temp = temp->next;temp = temp->next;
        //cout<<posterior(1,temp)<<endl;
        ////cout<<"done init"<<endl;
        //for(int i=0;i<14;i++)
        //   temp = temp->next;
        //cout<<posterior(1,temp)<<endl;
    
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
