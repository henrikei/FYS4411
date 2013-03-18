#include "wfgeneral.h"

wfGeneral::wfGeneral(const int &nPart, const double &a, const double &b)
{
    nParticles = nPart;
    alpha = a;
    beta = b;
    nDimensions = 3;
    slater = Slater(nParticles, alpha);
    jastrow = Jastrow(nParticles, beta);
}

void wfGeneral::setAlpha(const double &a){
    alpha = a;
    slater.setAlpha(alpha);
}

void wfGeneral::update(const mat &r){
    slater.update(r);
}

double wfGeneral::getRatio(const int &particleNum, const mat &rNew, const mat &rOld){
    return slater.getRatio(particleNum, rNew)*jastrow.getRatio(particleNum, rNew, rOld);
}

mat wfGeneral::getQuantumForceRatio(const mat &r){
    quantumForceRatioSD = slater.getQuantumForceRatio(r);
    quantumForceRatioJ = jastrow.getQuantumForceRatio(r);
    return quantumForceRatioSD + quantumForceRatioJ;
}

double wfGeneral::getLaplaceRatio(const mat &r, const double &h){
    double value = 0;
    for (int i = 0; i < nParticles; i++){
        value += (quantumForceRatioSD(i,0)*quantumForceRatioJ(i,0) + quantumForceRatioSD(i,1)*quantumForceRatioJ(i,1)
                    + quantumForceRatioSD(i,2)*quantumForceRatioJ(i,2));
    }
    value *= 2;
    value += slater.getLaplaceRatio(r) + jastrow.getLaplaceRatio(r);
    return value;
}
