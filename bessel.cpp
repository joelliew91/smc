#include "include.h"
#include <iostream>
#include <math.h>
#include <complex>
#include <complex_bessel.h>

#define PII 3.141592654
#define MAXIT 10000
#define XMIN 2.0
#define EPS 1.0e-16
#define NUSE1 5
#define NUSE2 5
#define FPMIN 1.0e-30

using namespace std;

double chebev(double a, double b, double c[], int m, double x)
{
    //void nrerror(char error_text[]);
    double d=0.0,dd=0.0,sv,y,y2;
    int j;
    
    if ((x-a)*(x-b) > 0.0) cout<<"x not in range in routine chebev"<<endl;
    y2=2.0*(y=(2.0*x-a-b)/(b-a));
    for (j=m-1;j>=1;j--) {
        sv=d;
        d=y2*d-dd+c[j];
        dd=sv;
    }
    return y*d-dd+0.5*c[0];
}

void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
    double chebev(double a, double b, double c[], int m, double x);
    double xx;
    static double c1[] = {
        -1.142022680371168e0,6.5165112670737e-3,
        3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
        3.67795e-11,-1.356e-13};
    static double c2[] = {
        1.843740587300905e0,-7.68528408447867e-2,
        1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
        2.423096e-10,-1.702e-13,-1.49e-15};
    
    xx=8.0*x*x-1.0;
    *gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
    *gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
    *gampl= *gam2-x*(*gam1);
    *gammi= *gam2+x*(*gam1);
}

double bessik(double x, double xnu)
{
    void beschb(double x, double *gam1, double *gam2, double *gampl,
                double *gammi);
    //void nrerror(char error_text[]);
    int i,l,nl;
    double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
    gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
    ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
    double val,rkp,rip,ri;
    
    if (x <= 0.0 || xnu < 0.0) cout<<"bad arguments in bessik"<<endl;
    nl=(int)(xnu+0.5);
    xmu=xnu-nl;
    xmu2=xmu*xmu;
    xi=1.0/x;
    xi2=2.0*xi;
    h=xnu*xi;
    if (h < FPMIN) h=FPMIN;
    b=xi2*xnu;
    d=0.0;
    c=h;
    for (i=1;i<=MAXIT;i++) {
        b += xi2;
        d=1.0/(b+d);
        c=b+1.0/c;
        del=c*d;
        h=del*h;
        if (fabs(del-1.0) < EPS) break;
    }
    if (i > MAXIT) cout<<"x too large in bessik; try asymptotic expansion"<<endl;
    ril=FPMIN;
    ripl=h*ril;
    ril1=ril;
    rip1=ripl;
    fact=xnu*xi;
    for (l=nl;l>=1;l--) {
        ritemp=fact*ril+ripl;
        fact -= xi;
        ripl=fact*ritemp+ril;
        ril=ritemp;
    }
    f=ripl/ril;
    if (x < XMIN) {
        x2=0.5*x;
        pimu=PII*xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
        d = -log(x2);
        e=xmu*d;
        fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
        beschb(xmu,&gam1,&gam2,&gampl,&gammi);
        ff=fact*(gam1*cosh(e)+gam2*fact2*d);
        sum=ff;
        e=exp(e);
        p=0.5*e/gampl;
        q=0.5/(e*gammi);
        c=1.0;
        d=x2*x2;
        sum1=p;
        for (i=1;i<=MAXIT;i++) {
            ff=(i*ff+p+q)/(i*i-xmu2);
            c *= (d/i);
            p /= (i-xmu);
            q /= (i+xmu);
            del=c*ff;
            sum += del;
            del1=c*(p-i*ff);
            sum1 += del1;
            if (fabs(del) < fabs(sum)*EPS) break;
        }
        if (i > MAXIT) cout<<"bessk series failed to converge"<<endl;
        rkmu=sum;
        rk1=sum1*xi2;
    } else {
        b=2.0*(1.0+x);
        d=1.0/b;
        h=delh=d;
        q1=0.0;
        q2=1.0;
        a1=0.25-xmu2;
        q=c=a1;
        a = -a1;
        s=1.0+q*delh;
        for (i=2;i<=MAXIT;i++) {
            a -= 2*(i-1);
            c = -a*c/i;
            qnew=(q1-b*q2)/a;
            q1=q2;
            q2=qnew;
            q += c*qnew;
            b += 2.0;
            d=1.0/(b+a*d);
            delh=(b*d-1.0)*delh;
            h += delh;
            dels=q*delh;
            s += dels;
            if (fabs(dels/s) < EPS) break;
        }
        if (i > MAXIT) cout<<"bessik: failure to converge in cf2"<<endl;
        h=a1*h;
        rkmu=sqrt(PII/(2.0*x))*exp(-x)/s;
        rk1=rkmu*(xmu+x+0.5-h)*xi;
    }
    rkmup=xmu*xi*rkmu-rk1;
    rimu=xi/(f*rkmu-rkmup);
    ri=(rimu*ril1)/ril;
    rip=(rimu*rip1)/ril;
    for (i=1;i<=nl;i++) {
        rktemp=(xmu+i)*xi2*rk1+rkmu;
        rkmu=rk1;
        rk1=rktemp;
    }
    
    rkp=xnu*xi*rkmu-rk1;
    
    val = rkmu;
    
    return(val);
}
/*
complex<double> besselI(double v, complex<double> z){
    complex<double> val1(0.0,0.0);
    
    if(norm(z)>700.0){
        complex<double> val2,val3;
        complex<double> i(0.0,1.0);
        val3 = pow(z,v)*pow(-pow(z,2),(-2*v+1)/4)/sqrt(2*PII);
        val2 = exp(-i*(PII*(2*v+1)/4-sqrt(-pow(z,2))))+exp(i*((2*v+1)/4-sqrt(-pow(z,2))));
        val1 = val3*val2;
        
        cout<<z<<" "<<val1<<" "<<val2<<" "<<-i*(PII*(2*v+1)/4-sqrt(-pow(z,2)))<<" "<<i*((2*v+1)/4-sqrt(-pow(z,2)))<<endl;
    }
    else{
        val1 = sp_bessel::besselI(v,z);
        //cout<<"sp_bessel"<<endl;
    }
    return val1;
}
*/









