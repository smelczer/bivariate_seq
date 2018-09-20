#include "bivariate_lin_seq.h"
#include <NTL/BasicThreadPool.h>

#include <NTL/lzz_pX.h>
NTL_CLIENT

int main(){
    // set this to be the prime you want to work over
    long p = 2013265921;
    zz_p::init(p);
    SetNumThreads(4);
    
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

    long D,N;
    long DD = 1000; 
    long NN = 1000;
    long m = 5;

    printf("Terms (%ld,%ld) to (%ld,%ld):\n\n",NN-m,DD-m,NN+m,DD+m);

    zz_pX n,d; // this is num/den for the output
    for(D=DD-m; D < DD+m+1; D++){
        bls.find_row(n,d,D);

        Vec<zz_p> init = get_init(deg(d), n, d);
        Vec<zz_p> coeffs;
        coeffs.SetLength(2*m+1);

        for(N=NN-m; N < NN+m+1; N++){
            coeffs[N-NN+m] = get_elem(N, reverse(d), init);
        }
        cout << coeffs << endl;
    }
    
}
