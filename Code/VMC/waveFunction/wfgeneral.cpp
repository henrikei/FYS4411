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
    gradientRatioSD = slater.getGradientRatio(r);
    gradientRatioJ = jastrow.getGradientRatio(r);
    return 2*(gradientRatioSD + gradientRatioJ);
}

double wfGeneral::getLaplaceRatio(const mat &r, const double &h){
    double value = 0;
    for (int i = 0; i < nParticles; i++){
        value += (gradientRatioSD(i,0)*gradientRatioJ(i,0) + gradientRatioSD(i,1)*gradientRatioJ(i,1)
                    + gradientRatioSD(i,2)*gradientRatioJ(i,2));
    }
    value *= 2;
    value += slater.getLaplaceRatio(r) + jastrow.getLaplaceRatio(r);
    return value;
}
