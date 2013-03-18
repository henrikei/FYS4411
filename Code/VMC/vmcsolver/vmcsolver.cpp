#include "vmcsolver.h"

VMCSolver::VMCSolver() :
    nParticles(-1),
    nDimensions(-1),
    charge(-1),
    h(0.001),
    idum(-time(NULL)),
    nCycles(100000)
{
}

void VMCSolver:: setCharge(const int &charg){
    charge = charg;
}

void VMCSolver::setWaveFunction(waveFunction *w){
    wf = w;
    nParticles = w->getNParticles();
    nDimensions = w->getNDimensions();
}

void VMCSolver::setLocalEnergy(const localEnergy &E){
    localE = E;
}

double VMCSolver::getEnergy(){
    return energy;
}

double VMCSolver::getVariance(){
    return variance;
}
