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

    return (hydrogenic->getValue(particleNum, quantumNum, R + rNuclei)
            + hydrogenic->getValue(particleNum, quantumNum, R - rNuclei));
}

rowvec3 Diatomic::getGradient(int particleNum, int quantumNum, const mat &R){

    return (hydrogenic->getGradient(particleNum, quantumNum, R + rNuclei)
            + hydrogenic->getGradient(particleNum, quantumNum, R - rNuclei));
}

double Diatomic::getLaplacian(int particleNum, int quantumNum, const mat &R){

    return (hydrogenic->getLaplacian(particleNum, quantumNum, R + rNuclei)
            + hydrogenic->getLaplacian(particleNum, quantumNum, R - rNuclei));
}

double Diatomic::getAlphaDerivative(int particleNum, int quantumNum, const mat &R){
    return (hydrogenic->getAlphaDerivative(particleNum, quantumNum, R + rNuclei)
            + hydrogenic->getAlphaDerivative(particleNum, quantumNum, R - rNuclei));
}
