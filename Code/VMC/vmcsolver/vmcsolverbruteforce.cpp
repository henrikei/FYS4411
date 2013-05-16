#include "vmcsolverbruteforce.h"
#include "lib.h"
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"
#include "Slater/slater.h"
#include "orbitals/orbitals.h"

#include <armadillo>
#include <iostream>
#include <sys/time.h>
#include <mpi.h>

using namespace arma;
using namespace std;

VMCSolverBruteForce::VMCSolverBruteForce()
{
    stepLength = 0;
    nDummyCycles = 1000;
}

void VMCSolverBruteForce::runMonteCarloIntegration()
{
    rOld = zeros<mat>(nParticles, nDimensions);
    rNew = zeros<mat>(nParticles, nDimensions);
    double ratio2 = 0;

    double energySum = 0;
    double energySquaredSum = 0;
    double deltaE;

    double dPsidAlpha = 0;
    double dPsidBeta = 0;
    double alphaTerm1 = 0;
    double alphaTerm2 = 0;
    double betaTerm1 = 0;
    double betaTerm2 = 0;

    int acceptCount = 0;

    ofstream outFileBlocking;
    ofstream outFileOneBody;
    stringstream name;

    name << "../Out/blocking" << my_rank << ".dat";
    outFileBlocking.open(name.str().c_str());


    // initial trial positions
    double initialLength = 10;
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = initialLength*(randu() - 0.5);
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
                    rNew(i,j) = rOld(i,j) + stepLength*(randu() - 0.5);
                }
                ratio2 = wf->getRatio(i, rNew, rOld);
                ratio2 *= ratio2;
                if (randu() < ratio2){
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
                rNew(i,j) = rOld(i,j) + stepLength*(randu() - 0.5);
            }

            // Recalculate the value of the wave function
            ratio2 = wf->getRatio(i, rNew, rOld);
            ratio2 *= ratio2;

            // Check for step acceptance (if yes, update position, if no, reset position)
            if(randu() <= ratio2) {
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                }
                wf->update(rNew);
                acceptCount += 1;
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
            }
            // update energies
            if (cycle >= thermalization){
                deltaE = localE->getValue(rNew, wf, charge);
                energySum += deltaE;
                energySquaredSum += deltaE*deltaE;

                // if minimizing, calculate energy gradients
                if(minimizer){
                    dPsidAlpha = wf->getAlphaDerivativeRatio(rNew);
                    alphaTerm1 += dPsidAlpha;
                    alphaTerm2 += deltaE*dPsidAlpha;

                    dPsidBeta = wf->getBetaDerivativeRatio(rNew);
                    betaTerm1 += dPsidBeta;
                    betaTerm2 += deltaE*dPsidBeta;
                }

                // if oneBody, write out the radial distances of all particles every 20th cycle
                if(oneBody && (cycle + i) % 20 == 0){
                    if(oneBody < 2){
                        stringstream filename;
                        filename << "../Out/onebody" << my_rank << ".dat";
                        outFileOneBody.open(filename.str().c_str());
                        oneBody = 2;
                    }
                    for (int j = 0; j < nParticles; j++){
                        for (int k = 0; k < nDimensions; k++){
                            outFileOneBody << rNew(j,k) << "  ";
                        }
                        outFileOneBody << endl;
                    }
                }
            }
        }
    }

    outFileBlocking.close();

    if(oneBody){
        outFileOneBody.close();
    }

    //ALLREDUCE
    double energySumMPI, energySquaredSumMPI;
    MPI_Allreduce(&energySum, &energySumMPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&energySquaredSum, &energySquaredSumMPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    energySum = energySumMPI/numprocs;
    energySquaredSum = energySumMPI/numprocs;


    energy = energySum/(nCycles * nParticles); //SKALER MED NPROCS OGSAA
    variance = energySquaredSum/(nCycles * nParticles) - energy*energy;

    if(minimizer ==1){

        double alphaTerm1MPI, alphaTerm2MPI;
        MPI_Allreduce(&alphaTerm1, &alphaTerm1MPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&alphaTerm2, &alphaTerm2MPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //maa vite numprocs
        alphaTerm1 = alphaTerm1MPI/numprocs;
        alphaTerm2 = alphaTerm2MPI/numprocs;

        alphaTerm1 = alphaTerm1/(nCycles*nParticles);
        alphaTerm2 = alphaTerm2/(nCycles*nParticles);

        betaTerm1 = betaTerm1/(nCycles*nParticles);
        betaTerm2 = betaTerm2/(nCycles*nParticles);

        dEdAlpha = 2*(alphaTerm2 - alphaTerm1*energy);
        dEdBeta = 2*(betaTerm2 - betaTerm1*energy);
    }

    //cout << "Acceptance ratio: " << (double)acceptCount/((double)((nCycles + thermalization)*nParticles)) << endl;
}
