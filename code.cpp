#include <iostream>
#include <math.h>
#include "include.h"
#include <random>
#include <vector>
#include <chrono>
#include <complex_bessel.h>
#include <time.h>
#include <complex>
#include "rnglib.h"
#include "ranlib.h"

double TOL = 0.000001;
double ep = 0.0001;

using namespace std;

para* head = NULL;
double max_v;
double max_w;
double min_v;
double min_w;
int no_particles=100;
double dmax = 2;
double delta=1;
vector<double> price;
int total;
int flag = 1;


int main(){
    srand(10);
    initialize();
    para* temp = head;
    //complex<double> a(2800000.0,0.0);
    //complex<double> t = sp_bessel::besselK(delta/delta-0.5,a);
    //cout<<t<<endl;
    //temp = temp->next;temp = temp->next;
    //cout<<posterior(1,temp)<<endl;
    //while(temp != NULL){
     //   cout<<posterior(1,temp)<<endl;
       // temp = temp->next;
    //}
    //print();
    //cout<<log(-1)<<endl;
	return(0);
}
