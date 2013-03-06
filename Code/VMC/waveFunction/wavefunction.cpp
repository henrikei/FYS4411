#include "wavefunction.h"

waveFunction::waveFunction():
    alpha(128),
    beta(247),
    nParticles(-1),
    nDimensions(-1),
    quantumForce(zeros(2,2))
{
}

void waveFunction::setAlpha(const double &a){
    alpha = a;
}

void waveFunction::setBeta(const double &b){
    beta = b;
}

int waveFunction::getNParticles(){
    return nParticles;
}

int waveFunction::getNDimensions(){
    return nDimensions;
}

double waveFunction::getLaplacian(const mat &r, const double &h){
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
