#include "berylliumsimplenum.h"

BerylliumSimpleNum::BerylliumSimpleNum()
{
    nParticles = 4;
    nDimensions = 3;
    quantumForce = zeros(nParticles, nDimensions);
}

double BerylliumSimpleNum::getValue(const mat &r){
    double value = 0;
    args.clear();
    for (int i = 0; i < nParticles; i++){
        double rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++){
            rSingleParticle += r(i,j)*r(i,j);
        }
        args.push_back(sqrt(rSingleParticle));
    }
    value = (psi1S(args.at(0))*psi2S(args.at(1)) - psi1S(args.at(1))*psi2S(args.at(0)))*
            (psi1S(args.at(2))*psi2S(args.at(3)) - psi1S(args.at(3))*psi2S(args.at(2)));
    return value;
}

double BerylliumSimpleNum::psi1S(const double &rSingleParticle){
    return exp(-alpha*rSingleParticle);
}

double BerylliumSimpleNum::psi2S(const double &rSingleParticle){
    return (1 - alpha*rSingleParticle/2)*exp(-alpha*rSingleParticle/2);
}
