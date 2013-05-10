#include "diatomicham.h"

DiatomicHam::DiatomicHam()
{
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
        for (int j = 0; j < nDimensions){
            rp1 += (r(i,j) + rNuclei(j))*(r(i,j) + rNuclei(j));
            rp2 += (r(i,j) - rNuclei(j))*(r(i,j) - rNuclei(j));
        }
        potentialEnergy += -charge*(1/sqrt(rp1) + 1/sqrt(rp2));
    }

}

void DiatomicHam::setNucleiDistance(double r){
    rNuclei << r/2 << 0.0 << 0.0;
}
