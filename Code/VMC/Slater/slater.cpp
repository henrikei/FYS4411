#include "slater.h"
#include "stdlib.h"

using namespace std;


Slater::Slater()
{
}

Slater::Slater(string orbitalType, int nPart)
{
    nParticles = nPart;
    nDimensions = 3;

    slaterUp = zeros(nParticles/2, nParticles/2);
    slaterDown = zeros(nParticles/2, nParticles/2);
    invSlaterUp = zeros(nParticles/2, nParticles/2);
    invSlaterDown = zeros(nParticles/2, nParticles/2);

    if (orbitalType == "Hydrogenic"){
        orbitals = new Hydrogenic(nParticles);
    } else if (orbitalType == "Diatomic"){
        orbitals = new Diatomic(nParticles);
    } else {
        cout << "Error: Orbital type not defined." << endl;
        exit(1);
    }
}

void Slater::setAlpha(double alpha){
    orbitals->setAlpha(alpha);
}

void Slater::setR(double R){
    orbitals->setR(R);
}

void Slater::update(const mat &R){
    for (int i = 0; i < nParticles/2; i++){
        for (int j = 0; j < nParticles/2; j++){
            slaterUp(i,j) = orbitals->getValue(i, j, R);
            slaterDown(i,j) = orbitals->getValue(nParticles/2 + i, j, R);
        }
    }
    try{
        invSlaterUp = inv(slaterUp);
        invSlaterDown = inv(slaterDown);
    }catch(std::runtime_error){
        invSlaterUp = zeros(nParticles/2, nParticles/2);
        invSlaterDown = zeros(nParticles/2, nParticles/2);
        cout << slaterUp << endl;
        cout << slaterDown << endl;
    }
}


double Slater::getRatio(const int &particleNum, const mat &R){
    double value = 0;
    // Find which slater determinant the particle belongs to
    if (particleNum < nParticles/2){
        for (int j = 0; j < nParticles/2; j++){
            value += orbitals->getValue(particleNum, j, R)*invSlaterUp(j, particleNum);
        }
    } else {
        for (int j = 0; j < nParticles/2; j++){
            value += orbitals->getValue(particleNum, j, R)*invSlaterDown(j, particleNum - nParticles/2);
        }
    }
    return value;
}


mat Slater::getGradientRatio(const mat &R){
    gradientRatio = zeros(nParticles, nDimensions);
    for (int i = 0; i < nParticles/2; i++){
        for (int j = 0; j < nParticles/2; j++){
            gradientRatio.row(i) += orbitals->getGradient(i, j, R)*invSlaterUp(j, i);
            gradientRatio.row(i + nParticles/2) += orbitals->getGradient(i + nParticles/2, j, R)*invSlaterDown(j, i);
        }
    }
    gradientRatio = gradientRatio;
    return gradientRatio;
}


double Slater::getLaplaceRatio(const mat &R){
    double value = 0;
    for (int i = 0; i < nParticles/2; i++){
        for (int j = 0; j < nParticles/2; j++){
            value += orbitals->getLaplacian(i, j, R)*invSlaterUp(j, i)
                    + orbitals->getLaplacian(i + nParticles/2, j, R)*invSlaterDown(j, i);
        }
    }
    return value;
}


double Slater::getAlphaDerivativeRatio(const mat &R){
    double value = 0;
    for (int i = 0; i < nParticles/2; i++){
        for (int j = 0; j < nParticles/2; j++){
            value += orbitals->getAlphaDerivative(i, j, R)*invSlaterUp(j, i)
                    + orbitals->getAlphaDerivative(i + nParticles/2, j, R)*invSlaterDown(j, i);
        }
    }
    return value;
}
