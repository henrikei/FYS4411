#include "vmcsolverbruteforce.h"
#include "lib.h"
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"
#include "Slater/slater.h"
#include "orbitals/orbitals.h"

#include <armadillo>
#include <iostream>
#include <sys/time.h>

using namespace arma;
using namespace std;

VMCSolverBruteForce::VMCSolverBruteForce(const int &charg)
{
    stepLength = 0;
    nDummyCycles = 1000;
    charge = charg;
}

void VMCSolverBruteForce::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);

    double ratio2 = 0;

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
    wf->update(rOld);

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
                ratio2 = wf->getRatio(i, rNew, rOld);
                ratio2 *= ratio2;
                if (ran2(&idum) < ratio2){
                    for (int j = 0; j < nDimensions; j++){
                        rOld(i,j) = rNew(i,j);
                    }
                    wf->update(rNew);
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
    for(int cycle = 0; cycle < (nCycles + thermalization); cycle++) {

        // Store the current value of the wave function
        wf->update(rOld);

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + stepLength*(ran2(&idum) - 0.5);
            }

            // Recalculate the value of the wave function
            ratio2 = wf->getRatio(i, rNew, rOld);
            ratio2 *= ratio2;

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(ran2(&idum) <= ratio2) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                }
                wf->update(rNew);
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }
            // update energies
            if (cycle >= thermalization){
                deltaE = localE.getValue(rNew, wf, charge);
                energySum += deltaE;
                energySquaredSum += deltaE*deltaE;
            }
        }
    }
    energy = energySum/(nCycles * nParticles);
    variance = energySquaredSum/(nCycles * nParticles) - energy*energy;
}
