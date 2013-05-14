#include "diatomicham.h"

DiatomicHam::DiatomicHam()
{
    rNuclei << 1.0 << 0.0 << 0.0;
}

double DiatomicHam::getValue(const mat &r, waveFunction *wf, double charge){

    double kineticEnergy = -0.5*wf->getLaplaceRatio(r);

    double potentialEnergy = 0;
    double rp1 = 0;     // electron-proton1 distance
    double rp2 = 0;     // electron-proton2 distance
    double r12 = 0;     // eletron-electron distance
    int nParticles = r.n_rows;
    int nDimensions = r.n_cols;

    // Contribution from electron-nuclei interactions
    for (int i = 0; i< nParticles; i++){
        rp1 = 0;
        rp2 = 0;
        for (int j = 0; j < nDimensions; j++){
            rp1 += (r(i,j) + rNuclei(j))*(r(i,j) + rNuclei(j));
            rp2 += (r(i,j) - rNuclei(j))*(r(i,j) - rNuclei(j));
        }
        potentialEnergy += -charge*(1/sqrt(rp1) + 1/sqrt(rp2));
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


    // Contribution from nucleus-nucleus interaction
    potentialEnergy += charge*charge/(2*rNuclei(0));

    return (kineticEnergy + potentialEnergy);

}

void DiatomicHam::setR(double r){
    rNuclei << r/2 << 0.0 << 0.0;
}
