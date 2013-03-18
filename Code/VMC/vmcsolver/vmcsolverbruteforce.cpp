#include "vmcsolverbruteforce.h"
#include "lib.h"
#include "waveFunction/wavefunction.h"
#include "waveFunction/heliumSimpleNum.h"
#include "waveFunction/heliumJastrowNum.h"
#include "localenergy/localEnergy.h"
#include "Slater/slater.h"
#include "orbitals/orbitals.h"

#include <armadillo>
#include <iostream>
#include <sys/time.h>

using namespace arma;
using namespace std;

VMCSolverBruteForce::VMCSolverBruteForce() :
    stepLength(0),
    nDummyCycles(10000)
{
}

void VMCSolverBruteForce::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double waveFunctionOld = 0;
    double waveFunctionNew = 0;

    double energySum = 0;
    double energySquaredSum = 0;

    double deltaE;

    // initial trial positions
    double initialLength = 10;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = initialLength*(ran2(&idum) - 0.5);
        }
    }
    rNew = rOld;

    // Store the current value of the wave function
    waveFunctionOld = wf->getValue(rOld);


    // trial Monte Carlo loop to achieve acceptance rate of approximately 0.5
    stepLength = 0;
    double acceptRate = 0;
    while (abs(acceptRate - 0.5) > 0.05){
        int acceptCount = 0;
        stepLength += 0.01;
        for (int cycle = 0; cycle < nDummyCycles; cycle++){
            for (int i = 0; i < nParticles; i++){
                for (int j = 0; j < nDimensions; j++){
                    rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
                }
                waveFunctionNew = wf->getValue(rNew);
                if (ran2(&idum) < (waveFunctionNew*waveFunctionNew/(waveFunctionOld*waveFunctionOld))){
                    for (int j = 0; j < nDimensions; j++){
                        rOld(i,j) = rNew(i,j);
                    }
                    waveFunctionOld = waveFunctionNew;
                    acceptCount += 1;
                } else {
                    for (int j = 0; j < nDimensions; j++){
                        rNew(i,j) = rOld(i,j);
                    }
                }
            }
        }
        acceptRate = (double) acceptCount/(nDummyCycles*nParticles);
    }


    // actual Monte Carlo loop
    for(int cycle = 0; cycle < nCycles; cycle++) {

        // Store the current value of the wave function
        waveFunctionOld = wf->getValue(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            waveFunctionNew = wf->getValue(rNew);

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= (waveFunctionNew*waveFunctionNew) / (waveFunctionOld*waveFunctionOld)) {
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
