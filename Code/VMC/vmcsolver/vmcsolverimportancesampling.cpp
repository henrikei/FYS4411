#include "vmcsolverimportancesampling.h"
#include "lib.h"

VMCSolverImportanceSampling::VMCSolverImportanceSampling(const int &charg)
{
    timeStep = 0.001;
    charge = charg;
}

void VMCSolverImportanceSampling::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    double ratio2 = 0;

    mat quantumForceOld = zeros(nParticles, nDimensions);
    mat quantumForceNew = zeros(nParticles, nDimensions);

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE;

    // initial trial positions
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = randn()*sqrt(timeStep);
        }
    }
    rNew = rOld;

    wf->update(rOld);
    quantumForceOld = wf->getQuantumForceRatio(rOld);

    // Monte Carlo loop
    for(int cycle = 0; cycle < nCycles; cycle++){

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + randn()*sqrt(timeStep) + 0.5*quantumForceOld(i,j)*timeStep;
            }

            // Calculate slater-ratio and quantum force-ratio
            ratio2 = wf->getRatio(i, rNew, rOld);
            ratio2 *= ratio2;
            wf->update(rNew);
            quantumForceNew = wf->getQuantumForceRatio(rNew);

            // Check for step acceptance (if yes, update position and quantum force, if no, reset position)
            if(ran2(&idum) <= (getGreensFunctionRatio(rNew, rOld, quantumForceNew, quantumForceOld)*ratio2)){
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                }
                quantumForceOld = quantumForceNew;
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
                wf->update(rOld);
            }
            // update energies
            deltaE = localE.getValue(rNew, wf, h, charge);
            energySum += deltaE;
            energySquaredSum += deltaE*deltaE;
        }
    }
    energy = energySum/(nCycles * nParticles);
    variance = energySquaredSum/(nCycles * nParticles) - energy*energy;
}

// Calculates greens function ratio. y and x are position matrices. n is particle number.
double VMCSolverImportanceSampling::getGreensFunctionRatio(const mat &y, const mat &x, const mat &quantumForceNew, const mat &quantumForceOld){
    double argument1 = 0;
    double argument2 = 0;
    double argumentSum = 0;
    double greensFunctionRatio = 0;
    for (int i = 0; i < nParticles; i++){
        for (int j = 0; j < nDimensions; j++){
            argument1 += (y(i,j) - x(i,j) - 0.5*timeStep*quantumForceOld(i,j))*(y(i,j) - x(i,j) - 0.5*timeStep*quantumForceOld(i,j));
            argument2 += (x(i,j) - y(i,j) - 0.5*timeStep*quantumForceNew(i,j))*(x(i,j) - y(i,j) - 0.5*timeStep*quantumForceNew(i,j));
        }
    }
    argumentSum = (argument1 - argument2)/(2*timeStep);
    greensFunctionRatio = exp(argumentSum);
    return greensFunctionRatio;
}
