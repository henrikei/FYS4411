#include "vmcsolverimportancesampling.h"
#include "lib.h"
#include <mpi.h>



VMCSolverImportanceSampling::VMCSolverImportanceSampling()
{
    timeStep = 0.01;
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
    double potentialEnergySum = 0;
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
    for(int i = 0; i < nParticles; i++) {
        for(int j = 0; j < nDimensions; j++) {
            rOld(i,j) = randn()*sqrt(timeStep);
        }
    }
    rNew = rOld;

    wf->update(rOld);
    quantumForceOld = wf->getQuantumForceRatio();



    // Monte Carlo loop
    for(int cycle = 0; cycle < (nCycles + thermalization); cycle++){

        // New position to test
        for(int i = 0; i < nParticles; i++) {
            for(int j = 0; j < nDimensions; j++) {
                rNew(i,j) = rOld(i,j) + randn()*sqrt(timeStep) + 0.5*quantumForceOld(i,j)*timeStep;
            }

            // Calculate slater-ratio and quantum force-ratio
            ratio2 = wf->getRatio(i, rNew, rOld);
            ratio2 *= ratio2;
            // In importance sampling wf needs to be updated before accepting/rejecting the
            // move because quantumforce is needed.
            wf->update(rNew);
            quantumForceNew = wf->getQuantumForceRatio();

            // Check for step acceptance (if yes, update position and quantum force, if no, reset position)
            if(randu() <= (getGreensFunctionRatio(rNew, rOld, quantumForceNew, quantumForceOld)*ratio2)){
                for(int j = 0; j < nDimensions; j++) {
                    rOld(i,j) = rNew(i,j);
                }
                quantumForceOld = quantumForceNew;
                acceptCount += 1;
            } else {
                for(int j = 0; j < nDimensions; j++) {
                    rNew(i,j) = rOld(i,j);
                }
                wf->update(rOld);
            }


            // calculate integrals
            if (cycle >= thermalization){
                localE->calculate(rNew, wf, charge);
                deltaE = localE->getTotal();
                energySum += deltaE;
                energySquaredSum += deltaE*deltaE;

                potentialEnergySum += localE->getPotential();

                // Write deltaE to file for later Blocking
                outFileBlocking << deltaE << endl;

                // if minimizing, calculate energy gradients
                if(minimizer){
                    dPsidAlpha = wf->getAlphaDerivativeRatio(rNew);
                    alphaTerm1 += dPsidAlpha;
                    alphaTerm2 += deltaE*dPsidAlpha;

                    dPsidBeta = wf->getBetaDerivativeRatio(rNew);
                    betaTerm1 += dPsidBeta;
                    betaTerm2 += deltaE*dPsidBeta;
                }

                // if oneBody, write out the positions of all particles after all particles
                // have been moved
                if(oneBody && (i+1) == nParticles){
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

    double energySumMPI, energySquaredSumMPI, potentialEnergySumMPI;
    MPI_Allreduce(&energySum, &energySumMPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&energySquaredSum, &energySquaredSumMPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&potentialEnergySum, &potentialEnergySumMPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    energySum = energySumMPI/numprocs;
    energySquaredSum = energySquaredSumMPI/numprocs;
    potentialEnergySum = potentialEnergySumMPI/numprocs;

    energy = energySum/(nCycles * nParticles);
    variance = energySquaredSum/(nCycles * nParticles) - energy*energy;
    potentialEnergy = potentialEnergySum/(nCycles * nParticles);

    if(minimizer ==1){

        double alphaTerm1MPI, alphaTerm2MPI, betaTerm1MPI, betaTerm2MPI;
        MPI_Allreduce(&alphaTerm1, &alphaTerm1MPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&alphaTerm2, &alphaTerm2MPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&betaTerm1, &betaTerm1MPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&betaTerm2, &betaTerm2MPI, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //maa vite numprocs
        alphaTerm1 = alphaTerm1MPI/numprocs;
        alphaTerm2 = alphaTerm2MPI/numprocs;
        betaTerm1 = betaTerm1MPI/numprocs;
        betaTerm2 = betaTerm2MPI/numprocs;

        alphaTerm1 = alphaTerm1/(nCycles*nParticles);
        alphaTerm2 = alphaTerm2/(nCycles*nParticles);
        betaTerm1 = betaTerm1/(nCycles*nParticles);
        betaTerm2 = betaTerm2/(nCycles*nParticles);

        dEdAlpha = 2*(alphaTerm2 - alphaTerm1*energy);
        dEdBeta = 2*(betaTerm2 - betaTerm1*energy);
    }

    // cout << "Acceptance ratio: " << (double)acceptCount/((double)((nCycles + thermalization)*nParticles)) << endl;
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

void VMCSolverImportanceSampling::setTimeStep(const double &dt){
    timeStep = dt;
}
