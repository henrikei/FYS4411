#include "helium.h"
#include "wavefunction.h"

helium::helium() :
    alpha(138)
{
}

void helium::setAlpha (const double a){
    alpha = a;
}

double helium::getValue (const mat &r){
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

double helium::getLaplacian(const mat &r, const double &h){
    double laplacian = 0;
    mat rPlus = r;
    mat rMinus = r;
    int nParticles = r.n_rows;
    int nDimensions = r.n_cols;
    for (int i = 0; i < nParticles; i++){
        for (int j = 0; j < nDimensions; j++){
            rPlus(i, j) = r(i,j) + h;
            rMinus(i, j) = r(i,j) - h;

            laplacian += (getValue(rPlus) - 2*getValue(r) + getValue(rMinus));

            rPlus(i, j) = r(i, j);
            rMinus(i, j) = r(i, j);
        }
    }

    laplacian = laplacian/(h*h);
    return laplacian;
}
