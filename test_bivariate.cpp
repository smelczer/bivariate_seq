#include "bivariate_lin_seq.h"
#include <NTL/BasicThreadPool.h>

#include <NTL/lzz_pX.h>
NTL_CLIENT

int main(){
    // set this to be the prime you want to work over
    zz_p::init(2013265921);
    SetNumThreads(4);
    
    // D should be less than the char of the field due to interpolation
    long D = 8000;
    
    // this are "bivariate" polynomials
    Vec<zz_pX> num;
    Vec<zz_pX> den;
    zz_pX P;
    
    /**** NUMERATOR ********************************/
    // encode num: 1-x
    SetCoeff(P,0,1);
    SetCoeff(P,1,-1);
    num.append(P);
    
    P = zz_pX(0);
    
    /**** DENOMINATOR *******************************/
    // encode den: (1+x)^2 + (-x)y + y^2
    
    // (1+x)^2
    SetCoeff(P,0,1);
    SetCoeff(P,1,2);
    SetCoeff(P,2,1);
    den.append(P); P=zz_pX(0);
    
    // (-x)
    SetCoeff(P,1,-1);
    den.append(P); P=zz_pX(0);
    
    // 1
    SetCoeff(P,0,1);
    den.append(P);
    
    cout << "num: " << num << endl;
    cout << "den: " << den << endl;
    
    // this creates the object
    bivariate_lin_seq bls{num,den,2,2};
    
    zz_pX n,d; // this is num/den for the output
    double t = GetTime();
    bls.find_row(n,d,D);
    cout << "time: " << GetTime() - t << endl;
    //cout << "D: " << D << endl;
    //cout << "n: " << n << endl;
    //cout << "d: " << d << endl;
    
}
