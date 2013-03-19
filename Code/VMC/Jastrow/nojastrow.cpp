#include "nojastrow.h"

NoJastrow::NoJastrow(const int &nPart)
{
    nParticles = nPart;
    nDimensions = 3;
    gradientRatio = zeros(nParticles, nDimensions);
}

double NoJastrow::getRatio(const int &particleNum, const mat &rNew, const mat &rOld){
    (void)particleNum;
    (void)rNew;
    (void)rOld;
    return 1.0;
}

mat NoJastrow::getGradientRatio(const mat &r){
    (void)r;
    return gradientRatio;
}

double NoJastrow::getLaplaceRatio(const mat &r){
    (void)r;
    return 0.0;
}
