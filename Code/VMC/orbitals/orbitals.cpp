#include "orbitals.h"


Orbitals::Orbitals()
{
}

Orbitals::Orbitals(const int &nPart, const double &a)
{
    nParticles = nPart;
    nDimensions = 3;
    alpha = a;
    gradient = zeros<rowvec>(nDimensions);
}

void Orbitals::setAlpha(const double &a){
    alpha = a;
}

double Orbitals::getValue(const int &particleNum, const int &quantumNum, const mat &R){
    double rSingleParticle = 0;
    double value = 0;
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
        value = R(particleNum,0)*exp(-alpha*rSingleParticle/2);
    } else if (quantumNum == 3){
        value = R(particleNum,1)*exp(-alpha*rSingleParticle/2);
    } else if (quantumNum == 4){
        value = R(particleNum,2)*exp(-alpha*rSingleParticle/2);
    }
    return value;
}

rowvec3 Orbitals::getGradient(const int &particleNum, const int &quantumNum, const mat &R){
    double rSingleParticle = 0;
    for (int i = 0; i < nDimensions; i++){
        rSingleParticle += R(particleNum, i)*R(particleNum, i);
    }
    rSingleParticle = sqrt(rSingleParticle);

    // calculating the gradient of the appropriate orbital
    if (quantumNum == 0){
        gradient = -(alpha*R.row(particleNum)/rSingleParticle)*exp(-alpha*rSingleParticle);
    } else if (quantumNum == 1){
        gradient = (alpha*R.row(particleNum)*(alpha*rSingleParticle - 4)/(4*rSingleParticle))*exp(-alpha*rSingleParticle/2);
    } else if (quantumNum == 2){
        gradient(0) = -alpha*R(particleNum,0)*R(particleNum,0) + 2*rSingleParticle;
        gradient(1) = -alpha*R(particleNum,0)*R(particleNum,1);
        gradient(2) = -alpha*R(particleNum,0)*R(particleNum,2);
        gradient = gradient*exp(-alpha*rSingleParticle/2)/(2*rSingleParticle);
    } else if (quantumNum == 3){
        gradient(0) = -alpha*R(particleNum,0)*R(particleNum,1);
        gradient(1) = -alpha*R(particleNum,1)*R(particleNum,1) + 2*rSingleParticle;
        gradient(2) = -alpha*R(particleNum,1)*R(particleNum,2);
        gradient = gradient*exp(-alpha*rSingleParticle/2)/(2*rSingleParticle);
    } else if (quantumNum == 4){
        gradient(0) = -alpha*R(particleNum,0)*R(particleNum,2);
        gradient(1) = -alpha*R(particleNum,1)*R(particleNum,2);
        gradient(2) = -alpha*R(particleNum,2)*R(particleNum,2) + 2*rSingleParticle;
        gradient = gradient*exp(-alpha*rSingleParticle/2)/(2*rSingleParticle);
    }
    return gradient;
}

double Orbitals::getLaplacian(const int &particleNum, const int &quantumNum, const mat &R){
    double rSingleParticle = 0;
    double laplacian = 0;
    for (int i = 0; i < nDimensions; i++){
        rSingleParticle += R(particleNum, i)*R(particleNum, i);
    }
    rSingleParticle = sqrt(rSingleParticle);

    // calculating the laplacian of the appropriate orbital
    if (quantumNum == 0){
        laplacian = alpha*(alpha*rSingleParticle - 2)*exp(-alpha*rSingleParticle)/rSingleParticle;
    } else if (quantumNum == 1){
        laplacian = -alpha*(alpha*rSingleParticle - 8)*(alpha*rSingleParticle - 2)*exp(-alpha*rSingleParticle/2)/(8*rSingleParticle);
    } else if (quantumNum == 2){
        laplacian = alpha*R(particleNum,0)*exp(-alpha*rSingleParticle/2)*(alpha*rSingleParticle - 8)/(4*rSingleParticle);
    } else if (quantumNum == 3){
        laplacian = alpha*R(particleNum,1)*exp(-alpha*rSingleParticle/2)*(alpha*rSingleParticle - 8)/(4*rSingleParticle);
    } else if (quantumNum == 4){
        laplacian = alpha*R(particleNum,2)*exp(-alpha*rSingleParticle/2)*(alpha*rSingleParticle - 8)/(4*rSingleParticle);
    }
    return laplacian;
}

double Orbitals::getAlphaDerivative(const int &particleNum, const int &quantumNum, const mat &R){
    double rSingleParticle = 0;
    double value = 0;
    for (int i = 0; i < nDimensions; i++){
        rSingleParticle += R(particleNum, i)*R(particleNum, i);
    }
    rSingleParticle = sqrt(rSingleParticle);

    // calculating the derivative of the appropriate orbital
    if (quantumNum == 0){
        value = -rSingleParticle*exp(-alpha*rSingleParticle);
    } else if (quantumNum == 1){
        value = -rSingleParticle*(4 - alpha*rSingleParticle)*exp(-alpha*rSingleParticle/2)/4;
    } else if (quantumNum == 2){
        value = -0.5*R(particleNum,0)*rSingleParticle*exp(-alpha*rSingleParticle/2);
    } else if (quantumNum == 3){
        value = -0.5*R(particleNum,1)*rSingleParticle*exp(-alpha*rSingleParticle/2);
    } else if (quantumNum == 4){
        value = -0.5*R(particleNum,2)*rSingleParticle*exp(-alpha*rSingleParticle/2);
    }
    return value;
}
