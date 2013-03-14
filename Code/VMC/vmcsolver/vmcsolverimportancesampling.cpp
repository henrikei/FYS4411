#include "vmcsolverimportancesampling.h"
#include "lib.h"

VMCSolverImportanceSampling::VMCSolverImportanceSampling() :
    timeStep(0.01)
{
}

void VMCSolverImportanceSampling::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

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

    // Monte Carlo loop
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function and quantum force
        waveFunctionOld = wf->getValue(rOld);
        quantumForceOld = wf->getQuantumForce(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + randn()*sqrt(timeStep) + 0.5*quantumForceOld(i,j)*timeStep;
            }

            // Recalculate the value of the wave function and quantum force
            waveFunctionNew = wf->getValue(rNew);
            quantumForceNew = wf->getQuantumForce(rNew);

            // Check for step acceptance (if yes, update position and quantum force, if no, reset position)
            if(ran2(&idum) <= (getGreensFunctionRatio(rNew, rOld, quantumForceNew, quantumForceOld, timeStep)*waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
                    quantumForceOld = quantumForceNew;
                }
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
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
double VMCSolverImportanceSampling::getGreensFunctionRatio(const mat &y, const mat &x, const mat &quantumForceNew, const mat &quantumForceOld, const double & timeStep){
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
    greensFunctionRatio = exp(argumentSum)/pow(2*M_PI*timeStep, 3*nParticles/2);
    return greensFunctionRatio;
}





