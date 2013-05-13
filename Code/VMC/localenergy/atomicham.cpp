#include "atomicham.h"

AtomicHam::AtomicHam()
{
}

double AtomicHam::getValue(const mat &r, waveFunction *wf, double charge){

    double kineticEnergy = -0.5*wf->getLaplaceRatio(r);

    double potentialEnergy = 0;
    double rSingleParticle = 0;
    double r12 = 0;
    int nParticles = r.n_rows;
    int nDimensions = r.n_cols;

    // Contribution from nucleus
    for (int i = 0; i < nParticles; i++){
        rSingleParticle = 0;
        for (int j = 0; j < nDimensions; j++){
            rSingleParticle += r(i,j)*r(i,j);
        }
        potentialEnergy += -charge/sqrt(rSingleParticle);
    }

    // Contribution from electron-electron interactions
    for (int i = 0; i < nParticles; i++){
        for (int j = i + 1; j < nParticles; j++){
            r12 = 0;
            for (int k = 0; k < nDimensions; k++){
                r12 += (r(i,k) - r(j,k))*(r(i,k) - r(j,k));
            }
            potentialEnergy += 1/sqrt(r12);
        }
    }
    return (kineticEnergy + potentialEnergy);
}

void AtomicHam::setR(double r){
    (void) r;
}
