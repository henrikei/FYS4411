#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"

using namespace arma;

class VMCSolver
{
public:
    VMCSolver();

    void setWaveFunction(waveFunction *w);
    void setLocalEnergy(const localEnergy &);
    void runMonteCarloIntegration();
    double getEnergy();
    double getVariance();

private:
    // double localEnergy(const mat &r);
    waveFunction *wf;
    localEnergy localE;

    int nDimensions;
    int charge;
    double stepLength;
    int nParticles;

    double h;
    double h2;

    long idum;

    double alpha;
    double beta;

    int nCycles;

    mat rOld;
    mat rNew;

    double energy;
    double variance;
};

#endif // VMCSOLVER_H
