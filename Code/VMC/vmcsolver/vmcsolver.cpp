#include "vmcsolver.h"

VMCSolver::VMCSolver() :
    nParticles(-1),
    nDimensions(-1),
    charge(-1),
    h(0.01),
    idum(-time(NULL)),
    nCycles(1000000),
    thermalization(nCycles/10),
    minimizer(0),
    oneBody(0)
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

void VMCSolver::setLocalEnergy(localEnergy *E){
    localE = E;
}

double VMCSolver::getEnergy(){
    return energy;
}

double VMCSolver::getVariance(){
    return variance;
}

double VMCSolver::getdEdAlpha(){
    return dEdAlpha;
}

double VMCSolver::getdEdBeta(){
    return dEdBeta;
}

void VMCSolver::calcEnergyGradients(){
    minimizer = 1;
}

void VMCSolver::calcOneBodyDensity(){
    oneBody = 1;
}
