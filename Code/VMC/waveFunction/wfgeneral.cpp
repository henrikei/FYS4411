#include "wfgeneral.h"

wfGeneral::wfGeneral(const int &nPart, const double &alph)
{
    nParticles = nPart;
    alpha = alph;
    nDimensions = 3;
    slater = Slater(nParticles, alpha);
}

void wfGeneral::update(const mat &r){
    slater.update(r);
}

double wfGeneral::getRatio(const int &particleNum, const mat &r){
    return slater.getRatio(particleNum, r);
}

mat wfGeneral::getQuantumForceRatio(const mat &r){
    return slater.getQuantumForceRatio(r);
}

double wfGeneral::getLaplaceRatio(const mat &r, const double &h){
    return slater.getLaplaceRatio(r);
}
