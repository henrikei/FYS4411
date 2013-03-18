#include "heliumsimpleanalytic.h"

heliumSimpleAnalytic::heliumSimpleAnalytic()
{
    nParticles = 2;
    nDimensions = 3;
    quantumForce = zeros(2,3);
}

double heliumSimpleAnalytic::getValue (const mat &r){
    double argument = 0;
    int nParticles = r.n_rows;
    int nDimensions = r.n_cols;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument += sqrt(rSingleParticle);
    }

    return exp(-argument*alpha);
}

mat heliumSimpleAnalytic::getQuantumForceRatio (const mat &r){
    double rSingleParticle = 0;
    for (int i = 0; i < nParticles; i++){
        rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++){
            rSingleParticle += r(i,j)*r(i,j);
        }
        rSingleParticle = sqrt(rSingleParticle);
        quantumForce.row(i) = r.row(i)/rSingleParticle;
    }
    return -2*alpha*quantumForce;
}

double heliumSimpleAnalytic::getLaplacian(const mat &r, const double &h){
    double laplacian = 0;
    int nParticles = r.n_rows;
    int nDimensions = r.n_cols;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        laplacian -= 2*alpha/sqrt(rSingleParticle);
    }
    laplacian += 2*alpha*alpha;
    laplacian = laplacian*getValue(r);
    return laplacian;
}
