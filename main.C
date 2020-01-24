#include <iostream>
#include <cmath>
#include <complex>
#include <math.h>

using namespace std;

double gamma(complex<double> input_val){
    return 1;
}

double rho_k (double k, double alpha){
    return pow(pow(k, 2) - pow(alpha, 2), 0.5);
}

double Dk (double k, double q){
    double i = 1;
    double alpha = 1;

    double expoRatio = exp(- M_PI * i * k) / (k + i * q);
    double gammaRatioFirstTerm =  gamma(k - i * q) * pow(gamma(k + i * q), -1);
    double firstTerm = expoRatio * gammaRatioFirstTerm;

    double expoRatio = exp(- M_PI * i * rho_k(k, alpha)) / (k + i * q);
    double gammaRatioSecondTerm =  gamma(rho_k(k, alpha) - i * q) * pow(gamma(rho_k(k, alpha) + i * q), -1);
    double secondTerm = expoRatio * gammaRatioFirstTerm;

    double retVal = firstTerm + secondTerm;
    return retVal;
}

complex<double> F0(double q, double theta){
    complex<double> i (0, 1);
    complex<double> one (1, 0);

    complex<double> expoTerm = exp(i * q * log(pow(sin(theta/2), 2)));
    complex<double> retVal = 0.5 * i * expoTerm * gamma(one - i * q);
    return retVal;
};

complex<double> F1(){
    complex<double> i (0, 1);
    return 1;
};

complex<double> F(double q, double theta){
    return F0(q, theta) + F1();
};

complex<double> G0(double q, double theta){
    complex<double> i (0, 1);
    complex<double> retVal = - i * pow(q, 2) * pow(1/tan(theta), 2) * F0(q, theta);
    return retVal;
};

complex<double> G1(){
    complex<double> i (0, 1);
    complex<double> retVal = 0.5 * i ;
    return retVal;
};

complex<double> G(double q, double theta){
    return G0(q, theta) + G1();
};


double crossSection(double q, double beta, double theta) {
    double cs_firstTerm =  pow(q, 2) * (1 - pow(beta, 2)) * pow(norm(F(q, theta)), 2) / pow(sin(theta), 2);
    double cs_secondTerm =  pow(norm(G(q, theta)), 2) / pow(cos(theta), 2);
    double retVal = cs_firstTerm + cs_secondTerm;
    return retVal;
};

double fDk(double *x, double*p){
    return sin(x[0] * p[0]);
}

TF1 fsin("fsin", fDk, -10, 10, 1);

int plot(){
    fsin.SetParameter(0, 1);
    fsin.Draw();

    return 0;
}

int main() {
    double q = 1;
    double beta = 0.99;
    double theta = 1;
    cout << crossSection(q, beta, theta) << endl;

    return 0;
}
