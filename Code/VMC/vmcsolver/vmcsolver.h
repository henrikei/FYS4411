#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"

class VMCSolver
{
public:
    VMCSolver();
    void setCharge(const int &charg);
    void setWaveFunction(waveFunction *w);
    void setLocalEnergy(localEnergy *E);
    double getEnergy();
    double getVariance();
    double getdEdAlpha();
    double getdEdBeta();
    virtual void runMonteCarloIntegration()=0;
protected:
    waveFunction *wf;
    localEnergy *localE;

    int nParticles;
    int nDimensions;

    int charge;

    double h;
    double h2;

    long idum;

    double alpha;
    double beta;

    int nCycles;
    int thermalization;

    mat rOld;
    mat rNew;

    double energy;
    double variance;

    double dEdAlpha;
    double dEdBeta;
};

#endif // VMCSOLVER_H
