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

    double value;

    // return the correct orbital
    if (quantumNum % 2 == 0){
        value = hydrogenic->getValue(particleNum, quantumNum/2, R + rNuclei)
                + hydrogenic->getValue(particleNum, quantumNum/2, R - rNuclei);
    } else {
        value = hydrogenic->getValue(particleNum, quantumNum/2, R + rNuclei)
                 - hydrogenic->getValue(particleNum, quantumNum/2, R - rNuclei);
    }
    return value;
//    return (hydrogenic->getValue(particleNum, quantumNum, R + rNuclei)
//            + hydrogenic->getValue(particleNum, quantumNum, R - rNuclei));
}

rowvec3 Diatomic::getGradient(int particleNum, int quantumNum, const mat &R){

     // return correct orbital
    if (quantumNum % 2 == 0){
        gradient = hydrogenic->getGradient(particleNum, quantumNum/2, R + rNuclei)
                    + hydrogenic->getGradient(particleNum, quantumNum/2, R - rNuclei);
    } else {
        gradient = hydrogenic->getGradient(particleNum, quantumNum/2, R + rNuclei)
                    - hydrogenic->getGradient(particleNum, quantumNum/2, R - rNuclei);
    }

    return gradient;

//    return (hydrogenic->getGradient(particleNum, quantumNum, R + rNuclei)
//            + hydrogenic->getGradient(particleNum, quantumNum, R - rNuclei));
}

double Diatomic::getLaplacian(int particleNum, int quantumNum, const mat &R){

    double laplacian;

    // return correct orbital
    if (quantumNum % 2 == 0){
        laplacian =  hydrogenic->getLaplacian(particleNum, quantumNum/2, R + rNuclei)
                      + hydrogenic->getLaplacian(particleNum, quantumNum/2, R - rNuclei);
    } else {
        laplacian =  hydrogenic->getLaplacian(particleNum, quantumNum/2, R + rNuclei)
                      - hydrogenic->getLaplacian(particleNum, quantumNum/2, R - rNuclei);
    }

    return laplacian;

//    return (hydrogenic->getLaplacian(particleNum, quantumNum, R + rNuclei)
//            + hydrogenic->getLaplacian(particleNum, quantumNum, R - rNuclei));
}

double Diatomic::getAlphaDerivative(int particleNum, int quantumNum, const mat &R){

    double dPhidAlpha;

    // return the correct orbital
    if (quantumNum % 2 == 0){
        dPhidAlpha = hydrogenic->getAlphaDerivative(particleNum, quantumNum/2, R + rNuclei)
                + hydrogenic->getAlphaDerivative(particleNum, quantumNum/2, R - rNuclei);
    } else {
        dPhidAlpha = (hydrogenic->getAlphaDerivative(particleNum, quantumNum/2, R + rNuclei)
                - hydrogenic->getAlphaDerivative(particleNum, quantumNum/2, R - rNuclei));
    }

    return dPhidAlpha;

//    return (hydrogenic->getAlphaDerivative(particleNum, quantumNum, R + rNuclei)
//            + hydrogenic->getAlphaDerivative(particleNum, quantumNum, R - rNuclei));
}
