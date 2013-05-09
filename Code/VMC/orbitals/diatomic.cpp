#include "diatomic.h"

Diatomic::Diatomic(int nPart, double r, double a)
{
    nParticles = nPart;
    nDimensions = 3;
    rNuclei << r/2 << 0.0 << 0.0;
    alpha = a;
    gradient = zeros<rowvec>(nDimensions);
}

double Diatomic::getValue(int particleNum, int quantumNum, const mat &R){
    double rp1 = 0; // Distance from proton 1
    double rp2 = 0; // Distance from proton 2
    double value = 0;
    for (int i = 0; i < nDimensions; i++){
        rp1 += (R(particleNum,i) + rNuclei(i))*(R(particleNum,i) + rNuclei(i));
        rp2 += (R(particleNum,i) - rNuclei(i))*(R(particleNum,i) - rNuclei(i));
    }

    // calculate the value of the appropriate orbital
    if(quantumNum == 0){
        value = exp(-alpha*rp1) + exp(-alpha*rp2);
    }
    return value;
}

rowvec3 Diatomic::getGradient(int particleNum, int quantumNum, const mat &R){
    double rp1 = 0; // Distance from proton 1
    double rp2 = 0; // Distance from proton 2
    for (int i = 0; i < nDimensions; i++){
        rp1 += (R(particleNum,i) + rNuclei(i))*(R(particleNum,i) + rNuclei(i));
        rp2 += (R(particleNum,i) - rNuclei(i))*(R(particleNum,i) - rNuclei(i));
    }
    if (quantumNum == 0){
        gradient = -alpha*(R.row(particleNum) + rNuclei)*exp(-alpha*rp1)/rp1 - alpha*(R.row(particleNum) - rNuclei)*exp(-alpha*rp2)/rp2;
    }
    return gradient;
}

double Diatomic::getLaplacian(int particleNum, int quantumNum, const mat &R){
    double rp1 = 0; // Distance from proton 1
    double rp2 = 0; // Distance from proton 2
    double laplacian = 0;
    for (int i = 0; i < nDimensions; i++){
        rp1 += (R(particleNum,i) + rNuclei(i))*(R(particleNum,i) + rNuclei(i));
        rp2 += (R(particleNum,i) - rNuclei(i))*(R(particleNum,i) - rNuclei(i));
    }
    if(quantumNum == 0){
        laplacian = alpha*(alpha*rp1 - 2)*exp(-alpha*rp1)/rp1 + alpha*(alpha*rp2 - 2)*exp(-alpha*rp2)/rp2;
    }
    return laplacian;
}
