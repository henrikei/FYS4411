#include "wavefunction.h"

waveFunction::waveFunction(string orbitalType, int nPart, int jas)
{
    nParticles = nPart;
    nDimensions = 3;

    slater = Slater(orbitalType, nParticles);

    if (jas == 0){
        jastrow = new NoJastrow(nParticles);
    } else {
        jastrow = new Jastrow(nParticles);
    }
}

void waveFunction::setAlpha(double alpha){
    slater.setAlpha(alpha);
}

void waveFunction::setBeta(double beta){
    jastrow->setBeta(beta);
}

void waveFunction::setR(double R){
    slater.setR(R);
}

int waveFunction::getNParticles(){
    return nParticles;
}

int waveFunction::getNDimensions(){
    return nDimensions;
}

void waveFunction::update(const mat &r){
    slater.update(r);
    gradientRatioSD = slater.getGradientRatio(r);
    gradientRatioJ = jastrow->getGradientRatio(r);
}

double waveFunction::getRatio(int particleNum, const mat &rNew, const mat &rOld){
    return slater.getRatio(particleNum, rNew)*jastrow->getRatio(particleNum, rNew, rOld);
}

mat waveFunction::getQuantumForceRatio(){
    return 2*(gradientRatioSD + gradientRatioJ);
}

double waveFunction::getLaplaceRatio(const mat &r){
    double value = 0;
    for (int i = 0; i < nParticles; i++){
        value += (gradientRatioSD(i,0)*gradientRatioJ(i,0) + gradientRatioSD(i,1)*gradientRatioJ(i,1)
                    + gradientRatioSD(i,2)*gradientRatioJ(i,2));
    }
    value *= 2;
    value += slater.getLaplaceRatio(r) + jastrow->getLaplaceRatio(r);
    return value;
}

double waveFunction::getAlphaDerivativeRatio(const mat &r){
    return slater.getAlphaDerivativeRatio(r);
}

double waveFunction::getBetaDerivativeRatio(const mat &r){
    return jastrow->getBetaDerivativeRatio(r);
}
