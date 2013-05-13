#include "diatomic.h"
#include "hydrogenic.h"

Diatomic::Diatomic(int nPart)
{
    nParticles = nPart;
    nDimensions = 3;
    rNuclei = zeros(nParticles, nDimensions); //vectors from origo to protons
    gradient = zeros<rowvec>(nDimensions);
    hydrogenic = new Hydrogenic(nParticles);
}

void Diatomic::setAlpha(double a){
    alpha = a;
    hydrogenic->setAlpha(a);
}

void Diatomic::setR(double R){
    for (int i = 0; i < nParticles; i++){
        rNuclei(i,0) = R/2;
    }
}

double Diatomic::getValue(int particleNum, int quantumNum, const mat &R){
//    double rp1 = 0; // electron distance to proton 1
//    double rp2 = 0; // electron distance to proton 2
//    double value = 0;
//    for (int i = 0; i < nDimensions; i++){
//        rp1 += (R(particleNum,i) + rNuclei(i))*(R(particleNum,i) + rNuclei(i));
//        rp2 += (R(particleNum,i) - rNuclei(i))*(R(particleNum,i) - rNuclei(i));
//    }
//    rp1 = sqrt(rp1);
//    rp2 = sqrt(rp2);

//    // calculate the value of the appropriate orbital
//    if(quantumNum == 0){
//        value = exp(-alpha*rp1) + exp(-alpha*rp2);
//    }

    return (hydrogenic->getValue(particleNum, quantumNum, R + rNuclei)
            + hydrogenic->getValue(particleNum, quantumNum, R - rNuclei));
}

rowvec3 Diatomic::getGradient(int particleNum, int quantumNum, const mat &R){
//    double rp1 = 0; // Distance from proton 1
//    double rp2 = 0; // Distance from proton 2
//    for (int i = 0; i < nDimensions; i++){
//        rp1 += (R(particleNum,i) + rNuclei(0,i))*(R(particleNum,i) + rNuclei(0,i));
//        rp2 += (R(particleNum,i) - rNuclei(0,i))*(R(particleNum,i) - rNuclei(0,i));
//    }
//    rp1 = sqrt(rp1);
//    rp2 = sqrt(rp2);

//    if (quantumNum == 0){
//        gradient = -alpha*(R.row(particleNum) + rNuclei)*exp(-alpha*rp1)/rp1 - alpha*(R.row(particleNum) - rNuclei)*exp(-alpha*rp2)/rp2;
//    }
    return (hydrogenic->getGradient(particleNum, quantumNum, R + rNuclei)
            + hydrogenic->getGradient(particleNum, quantumNum, R - rNuclei));
}

double Diatomic::getLaplacian(int particleNum, int quantumNum, const mat &R){
//    double rp1 = 0; // Distance from proton 1
//    double rp2 = 0; // Distance from proton 2
//    double laplacian = 0;
//    for (int i = 0; i < nDimensions; i++){
//        rp1 += (R(particleNum,i) + rNuclei(0,i))*(R(particleNum,i) + rNuclei(0,i));
//        rp2 += (R(particleNum,i) - rNuclei(0,i))*(R(particleNum,i) - rNuclei(0,i));
//    }
//    rp1 = sqrt(rp1);
//    rp2 = sqrt(rp2);

//    if(quantumNum == 0){
//        laplacian = alpha*(alpha*rp1 - 2)*exp(-alpha*rp1)/rp1 + alpha*(alpha*rp2 - 2)*exp(-alpha*rp2)/rp2;
//    }
    return (hydrogenic->getLaplacian(particleNum, quantumNum, R + rNuclei)
            + hydrogenic->getLaplacian(particleNum, quantumNum, R - rNuclei));
}

double Diatomic::getAlphaDerivative(int particleNum, int quantumNum, const mat &R){
    (void) particleNum;
    (void) quantumNum;
    (void) R;
    return 0;
}
