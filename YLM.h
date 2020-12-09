// taken from: github.com/srmnitc/pyscal/blob/2.10.1/src/pyscal/system.cpp
// starting line
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>
#include <list>
const double PI = 3.14159265354;
double dfactorial(int l, int m);
double PLM(int l, int m, double x);
void convert_to_spherical_coordinates(double x, double y, double z, double &r, double &phi, double &theta);
void YLM(int l, int m, double theta, double phi, double &realYLM, double &imgYLM);
void QLM(int l,int m,double theta,double phi,double &realYLM, double &imgYLM );
//---------------------------------------------------
// Methods for q calculation
using namespace std; 
//---------------------------------------------------
double dfactorial(int l,int m){
  //(l-m)!/(l+m)!
    double fac = 1.00;
    for(int i=0;i<2*m;i++){
        fac*=double(l+m-i);
    }
    return (1.00/fac);
}

double PLM(int l, int m, double x){
  //p_{m+1}^{m+1}(x) = -(2*m+1)*sqrt(1-x^2)p_m^m(x)
    double fact,pll,pmm,pmmp1,somx2;
    int i,ll;
    pll = 0.0;
    if (m < 0 || m > l || fabs(x) > 1.0)
        std::cerr << "impossible combination of l and m" << "\n";
    pmm=1.0;
    if (m > 0){
        somx2=sqrt((1.0-x)*(1.0+x));
        fact=1.0;
        for (i=1;i<=m;i++){
            pmm *= -fact*somx2;
            fact += 2.0;
        }
    }

    if (l == m)
        return pmm;
    else{
        pmmp1=x*(2*m+1)*pmm;
        if (l == (m+1))
            return pmmp1;
        else{
            for (ll=m+2;ll<=l;ll++){
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
            pmm=pmmp1;
            pmmp1=pll;
            }
        return pll;
        }
    }
}

void convert_to_spherical_coordinates(double x, double y, double z, double &r, double &phi, double &theta){
    r = sqrt(x*x+y*y+z*z);
    theta = acos(z/r);
    phi = atan2(y,x);
}

void YLM(int l, int m, double theta, double phi, double &realYLM, double &imgYLM){
  //compute real and imag parts of Y_l^m

    double factor;
    double m_PLM;
    m_PLM = PLM(l,m,cos(theta));
    factor = ((2.0*double(l) + 1.0)/ (4.0*PI))*dfactorial(l,m);
    realYLM = sqrt(factor) * m_PLM * cos(double(m)*phi);
    imgYLM  = sqrt(factor) * m_PLM * sin(double(m)*phi);
}
void QLM(int l,int m,double theta,double phi,double &realYLM, double &imgYLM ){
  //computes real and imag parts of Y_l^m where m can be negative
    realYLM = 0.0;
    imgYLM = 0.0;
    if (m < 0) {
        YLM(l, abs(m), theta, phi, realYLM, imgYLM);
        realYLM = pow(-1.0,m)*realYLM;
        imgYLM = pow(-1.0,m)*imgYLM;
    }
    else{
        YLM(l, m, theta, phi, realYLM, imgYLM);
    }
//      cout << "l: " << l << endl; 
//      cout << "m: " << m << endl; 
//      cout << "theta: " << theta << endl; 
//      cout << "phi: " << phi << endl; 
//      cout << "realYLM in QLM: " << realYLM << endl; 
}

