#include "slater.h"

Slater::Slater()
{
}

Slater::Slater(const int &nPart, const double &alph)
{
    nParticles = nPart;
    nDimensions = 3;
    alpha = alph;

    slaterUp = zeros(nParticles/2, nParticles/2);
    slaterDown = zeros(nParticles/2, nParticles/2);
    invSlaterUp = zeros(nParticles/2, nParticles/2);
    invSlaterDown = zeros(nParticles/2, nParticles/2);

    orbitals = Orbitals(nParticles, alpha);
}

void Slater::update(const mat &R){
    for (int i = 0; i < nParticles/2; i++){
        for (int j = 0; j < nParticles/2; j++){
            slaterUp(i,j) = orbitals.getValue(i, j, R);
            slaterDown(i,j) = orbitals.getValue(nParticles/2 + i, j, R);
        }
    }
    invSlaterUp = inv(slaterUp);
    invSlaterDown = inv(slaterDown);
}


double Slater::getRatio(const int &particleNum, const mat &R){
    double value = 0;
    // Find which slater determinant the particle belongs to
    if (particleNum < nParticles/2){
        for (int j = 0; j < nParticles/2; j++){
            value += orbitals.getValue(particleNum, j, R)*invSlaterUp(j, particleNum);
        }
    } else {
        for (int j = 0; j < nParticles/2; j++){
            value += orbitals.getValue(particleNum, j, R)*invSlaterDown(j, particleNum - nParticles/2);
        }
    }
    return value;
}


mat Slater::getQuantumForceRatio(const mat &R){
    quantumForceRatio = zeros(nParticles, nDimensions);
    for (int i = 0; i < nParticles/2; i++){
        for (int j = 0; j < nParticles/2; j++){
            quantumForceRatio.row(i) += orbitals.getGradient(i, j, R)*invSlaterUp(j, i);
            quantumForceRatio.row(i + nParticles/2) += orbitals.getGradient(i + nParticles/2, j, R)*invSlaterDown(j, i);
        }
    }
    quantumForceRatio = 2*quantumForceRatio; //getRatio(particleNum, R);
    return quantumForceRatio;
}


double Slater::getLaplaceRatio(const mat &R){
    double value = 0;
    for (int i = 0; i < nParticles/2; i++){
        for (int j = 0; j < nParticles/2; j++){
            value += orbitals.getLaplacian(i, j, R)*invSlaterUp(j, i)
                    + orbitals.getLaplacian(i + nParticles/2, j, R)*invSlaterDown(j, i);
        }
    }
    return value;
}
