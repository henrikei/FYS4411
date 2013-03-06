#include "heliumJastrowNum.h"
#include "wavefunction.h"

heliumJastrowNum::heliumJastrowNum()
{
    nParticles = 2;
    nDimensions = 3;
    quantumForce = zeros(2,3);
}

double heliumJastrowNum::getValue(const mat &r){
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
