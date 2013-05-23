#ifndef VMCSOLVER_H
#define VMCSOLVER_H

#include <armadillo>
#include <fstream>
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"

using namespace std;
using namespace arma;


class VMCSolver
{
public:
    VMCSolver();
    void setCharge(const int &charg);
    void setWaveFunction(waveFunction *w);
    void setLocalEnergy(localEnergy *E);
    void setNumOfCycles(int n);
    double getEnergy();
    double getVariance();
    double getPotentialEnergy();
    double getdEdAlpha();
    double getdEdBeta();
    void calcEnergyGradients();
    void calcOneBodyDensity();
    virtual void runMonteCarloIntegration()=0;
    void init_MPI(int rank, int nprocs);
    int getMPIRank();

protected:
    waveFunction *wf;
    localEnergy *localE;

    int my_rank;
    int numprocs;

    int nParticles;
    int nDimensions;
    int charge;

    double h;
    double h2;

    int nCycles;
    int thermalization;

    int minimizer;
    int oneBody;

    mat rOld;
    mat rNew;

    double energy;
    double variance;
    double potentialEnergy;

    double dEdAlpha;
    double dEdBeta;
};

#endif // VMCSOLVER_H
