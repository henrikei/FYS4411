#include "jastrow.h"

Jastrow::Jastrow()
{
}

Jastrow::Jastrow(int nPart)
{
    nParticles = nPart;
    nDimensions = 3;
    beta = 1;
}

void Jastrow::setBeta(const double &b){
    beta = b;
}

double Jastrow::getRatio(const int &particleNum, const mat &rNew, const mat &rOld){
    double argument = 0;
    double fNew = 0;
    double fOld = 0;
    double r12 = 0;
    for (int i = 0; i < particleNum; i++){
        r12Vec = rNew.row(particleNum) - rNew.row(i);
        r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
        fNew = f(r12, particleNum, i);

        r12Vec = rOld.row(particleNum) - rOld.row(i);
        r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
        fOld = f(r12, particleNum, i);

        argument += fNew - fOld;
    }
    for (int i = particleNum + 1; i < nParticles; i++){
        r12Vec = rNew.row(particleNum) - rNew.row(i);
        r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
        fNew = f(r12, particleNum, i);

        r12Vec = rOld.row(particleNum) - rOld.row(i);
        r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
        fOld = f(r12, particleNum, i);

        argument += fNew - fOld;
    }
    return exp(argument);
}

mat Jastrow::getGradientRatio(const mat &r){
    gradientRatio = zeros(nParticles, nDimensions);
    for (int k = 0; k < nParticles; k++){
        for (int i = 0; i < k; i++){
            r12Vec = r.row(k) - r.row(i);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            gradientRatio.row(k) += r12Vec*dfdr(r12, k, i)/r12;
        }
        for (int i = k + 1; i < nParticles; i++){
            r12Vec = r.row(i) - r.row(k);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            gradientRatio.row(k) -= r12Vec*dfdr(r12, k, i)/r12;
        }
    }
    return gradientRatio;
}

double Jastrow::getLaplaceRatio(const mat &r){
    double laplaceRatio = 0;
    for (int k = 0; k < nParticles; k++){
        for (int i = 0; i < k; i++){
            r12Vec = r.row(k) - r.row(i);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            laplaceRatio += (nDimensions - 1)*dfdr(r12, k, i)/r12 + d2fdr2(r12, k, i);
        }
        for (int i = k + 1; i < nParticles; i++){
            r12Vec = r.row(k) - r.row(i);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            laplaceRatio += (nDimensions - 1)*dfdr(r12, k, i)/r12 + d2fdr2(r12, k, i);
        }
        laplaceRatio += gradientRatio(k,0)*gradientRatio(k,0) + gradientRatio(k,1)*gradientRatio(k,1)
                        + gradientRatio(k,2)*gradientRatio(k,2);
    }
    return laplaceRatio;
}

double Jastrow::getBetaDerivativeRatio(const mat &r){
    double value = 0;
    for (int k = 0; k < nParticles; k++){
        for (int i = 0; i < k; i++){
            r12Vec = r.row(k) - r.row(i);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            value -= aFactor(k,i)*r12*r12/((1 + beta*r12)*(1 + beta*r12));
        }
    }
    return value;
}


// Gives the factor a in the Jastrow-function
double Jastrow::aFactor(const int &particleNum1, const int &particleNum2){
    double a = 0;
    if (((particleNum1 < nParticles/2) && (particleNum2 < nParticles/2)) || ((particleNum1 >= nParticles/2) && (particleNum2 >= nParticles/2))){
        a = 0.25;
    } else {
        a = 0.50;
    }
    return a;
}

// f is the function in the exponential of the Jastrow-function
double Jastrow::f(const double &r12, const int &particleNum1, const int &particleNum2){
    return aFactor(particleNum1, particleNum2)*r12/(1 + beta*r12);
}

double Jastrow::dfdr(const double &r12, const int &particleNum1, const int &particleNum2){
    return aFactor(particleNum1, particleNum2)/((1 + beta*r12)*(1 + beta*r12));
}

double Jastrow::d2fdr2(const double &r12, const int &particleNum1, const int &particleNum2){
    return -2*aFactor(particleNum1, particleNum2)*beta/((1 + beta*r12)*(1 + beta*r12)*(1 + beta*r12));
}
