#include "heliumjastrowanalytic.h"

heliumJastrowAnalytic::heliumJastrowAnalytic()
{
    nParticles = 2;
    nDimensions = 3;
    quantumForce = zeros(2,3);
}

double heliumJastrowAnalytic::getValue(const mat &r){
    int nParticles = r.n_rows;
    int nDimensions = r.n_cols;
    double argument1 = 0;
    for(int i = 0; i < nParticles; i++) {
        double rSingleParticle = 0;
        for(int j = 0; j < nDimensions; j++) {
            rSingleParticle += r(i,j) * r(i,j);
        }
        argument1 += sqrt(rSingleParticle)*alpha;
    }
    double argument2 = 0;
    for(int i = 0; i < nParticles; i++) {
        for(int j = i + 1; j < nParticles; j++) {
            double r12 = 0;
            for(int k = 0; k < nDimensions; k++) {
                r12 += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
            }
            r12 = sqrt(r12);
            argument2 += r12/(2*(1 + beta*r12));
        }
    }
    return exp(-argument1 + argument2);
}

mat heliumJastrowAnalytic::getQuantumForceRatio(const mat &r){
    double rSingleParticle = 0;
    double r12 = 0;
    for (int i = 0; i < nParticles; i++){
        rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++){
            rSingleParticle += r(i,j)*r(i,j);
        }
        rSingleParticle = sqrt(rSingleParticle);
        quantumForce.row(i) = -2*alpha*r.row(i)/rSingleParticle;
    }
    for (int i = 0; i < nParticles; i++){
        for (int j = i+1; j < nParticles; j++){
            for (int k = 0; k < nDimensions; k++){
                r12 += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
            }
        }
    }
    r12 = sqrt(r12);
    double denominator = (1 + beta*r12)*(1 + beta*r12)*r12;
    quantumForce.row(0) += (r.row(0) - r.row(1))/denominator;
    quantumForce.row(1) += (r.row(1) - r.row(0))/denominator;

    return quantumForce;
}

double heliumJastrowAnalytic::getLaplacian(const mat &r, const double &h){
    double r1 = 0;
    double r2 = 0;
    double r12 = 0;
    rowvec r1Vec = r.row(0);
    rowvec r2Vec = r.row(1);
    double laplacian = 0;

    for (int i = 0; i < 3; i++){
        r1 += r(0,i)*r(0,i);
        r2 += r(1,i)*r(1,i);
        r12 += (r(0,i) - r(1,i))*(r(0,i) - r(1,i));
    }

    r1 = sqrt(r1);
    r2 = sqrt(r2);
    r12 = sqrt(r12);

    laplacian = -2*alpha*(1/r1 + 1/r2) + 2*alpha*alpha - (alpha*(r1 + r2)*(1 - dot(r1Vec,r2Vec)/(r1*r2))/r12 - 1/(2*(1 + beta*r12)*(1 + beta*r12))
                                                          - 2/r12 + 2*beta/(1 + beta*r12))/((1 + beta*r12)*(1 + beta*r12));
    laplacian = laplacian*getValue(r);
    return laplacian;
}
