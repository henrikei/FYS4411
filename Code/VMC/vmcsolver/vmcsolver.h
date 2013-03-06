#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"

class VMCSolver
{
public:
    VMCSolver();
    void setWaveFunction(waveFunction *w);
    void setLocalEnergy(const localEnergy &);
    double getEnergy();
    double getVariance();
    virtual void runMonteCarloIntegration()=0;
protected:
    waveFunction *wf;
    localEnergy localE;

    int nParticles;
    int nDimensions;

    int charge;

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
