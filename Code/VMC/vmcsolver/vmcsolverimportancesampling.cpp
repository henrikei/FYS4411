#include "vmcsolverimportancesampling.h"
#include "lib.h"

VMCSolverImportanceSampling::VMCSolverImportanceSampling() :
    timeStep(0.02)
{
}

void VMCSolverImportanceSampling::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

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

    // Store the current value of the wave function
    waveFunctionOld = wf->getValue(rOld);


    // Monte Carlo loop
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = wf->getValue(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + randn()*sqrt(timeStep) + 0.5*wf->getQuantumForce(rOld)(i,j)*timeStep;
            }

            // Recalculate the value of the wave function
            waveFunctionNew = wf->getValue(rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (getGreensFunctionRatio(rNew, rOld, i, timeStep)*waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                    waveFunctionOld = waveFunctionNew;
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
double VMCSolverImportanceSampling::getGreensFunctionRatio(const mat &y, const mat &x, const int &n, const double & timeStep){
    double argumentNew = 0;
    double argumentOld = 0;
    double argumentSum = 0;
    double greensFunctionRatio = 0;
    for (int i = 0; i < nDimensions; i++){
        argumentNew += (y(n,i) - x(n,i) - 0.5*timeStep*wf->getQuantumForce(x)(n,i))*(y(n,i) - x(n,i) - 0.5*timeStep*wf->getQuantumForce(x)(n,i));
        argumentOld += (x(n,i) - y(n,i) - 0.5*timeStep*wf->getQuantumForce(y)(n,i))*(x(n,i) - y(n,i) - 0.5*timeStep*wf->getQuantumForce(y)(n,i));
    }
    argumentSum = (- argumentNew + argumentOld)/(2*timeStep);
    greensFunctionRatio = exp(argumentSum)/pow(2*M_PI*timeStep, 3*nParticles/2);
    return greensFunctionRatio;
}





