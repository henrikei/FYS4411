#include "orbitals.h"

Orbitals::Orbitals()
{
    nParticles = 2;
    nDimensions = 3;
    alpha = 2;
    value = 0;
    gradient = zeros(nDimensions);
    laplacian = 0;
}

double Orbitals::getValue(const int &quantumNum, const int &particleNum, double const &alpha, const mat &R){
    double rSingleParticle = 0;
    for (int i = 0; i < nDimensions; i++){
        rSingleParticle += R(particleNum, i)*R(particleNum, i);
    }
    rSingleParticle = sqrt(rSingleParticle);

    // calculating the value of the appropriate orbital
    if (quantumNum == 0){
        value = exp(-alpha*rSingleParticle);
    } else if (quantumNum == 1){
        value = (1 - alpha*rSingleParticle/2)*exp(-alpha*rSingleParticle/2);
    } else if (quantumNum == 2){
        value = alpha*R(particleNum,0)*exp(-alpha*rSingleParticle/2);
    } else if (quantumNum == 3){
        value = alpha*R(particleNum,1)*exp(-alpha*rSingleParticle/2);
    } else if (quantumNum == 4){
        value = alpha*R(particleNum,2)*exp(-alpha*rSingleParticle/2);
    }
    return value;
}

vec Orbitals::getGradient(const int &quantumNum, const int &particleNum, const double &alpha, const mat &R){
    double rSingleParticle = 0;
    for (int i = 0; i < nDimensions; i++){
        rSingleParticle += R(particleNum, i)*R(particleNum, i);
    }
    rSingleParticle = sqrt(rSingleParticle);

    // calculating the gradient of the appropriate orbital
    if (quantumNum == 0){
        gradient = -(alpha*R.row(particleNum)/rSingleParticle)*exp(-alpha*rSingleParticle);
    } else if (quantumNum == 1){
        gradient = (alpha*R.row(particleNum)/(4*rSingleParticle))(alpha*rSingleParticle - 4);
    } else if (quantumNum == 2){
        gradient(0) = -alpha*R(particleNum,0)*R(particleNum,0) + 2*rSingleParticle;
        gradient(1) = -alpha*R(particleNum,0)*R(particleNum,1);
        gradient(2) = -alpha*R(particleNum,0)*R(particleNum,2);
        gradient = gradient/(2*rSingleParticle);
    } else if (quantumNum == 3){
        gradient(0) = -alpha*R(particleNum,0)*R(particleNum,1);
        gradient(1) = -alpha*R(particleNum,1)*R(particleNum,1) + 2*rSingleParticle;
        gradient(2) = -alpha*R(particleNum,1)*R(particleNum,2);
        gradient = gradient/(2*rSingleParticle);
    } else if (quantumNum == 4){
        gradient(0) = -alpha*R(particleNum,0)*R(particleNum,2);
        gradient(1) = -alpha*R(particleNum,1)*R(particleNum,2);
        gradient(2) = -alpha*R(particleNum,2)*R(particleNum,2) + 2*rSingleParticle;
        gradient = gradient/(2*rSingleParticle);
    }
}

vec Orbitals::getLaplacian(const int &quantumNum, const int &particleNum, const double &alpha, const mat &R){
    double rSingleParticle = 0;
    for (int i = 0; i < nDimensions; i++){
        rSingleParticle += R(particleNum, i)*R(particleNum, i);
    }
    rSingleParticle = sqrt(rSingleParticle);

    // calculating the laplacian of the appropriate orbital
    if (quantumNum == 0){
        laplacian = alpha*(alpha*rSingleParticle - 2)/rSingleParticle;
    } else if (quantumNum == 1){
        laplacian = -alpha*(alpha*rSingleParticle - 8)*(alpha*rSingleParticle - 2)/(8*rSingleParticle);
    } else if (quantumNum == 2){
        laplacian = alpha*R(particleNum,0)*(alpha*rSingleParticle - 8)/(4*rSingleParticle);
    } else if (quantumNum == 3){
        laplacian = alpha*R(particleNum,1)*(alpha*rSingleParticle - 8)/(4*rSingleParticle);
    } else if (qantumNum == 4){
        laplacian = alpha*R(particleNum,2)*(alpha*rSingleParticle - 8)/(4*rSingleParticle);
    }
}
