#include "heliumSimpleNum.h"
#include "wavefunction.h"

heliumSimpleNum::heliumSimpleNum()
{
    nParticles = 2;
    nDimensions = 3;
    quantumForce = zeros(2,3);
}

double heliumSimpleNum::getValue (const mat &r){
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
