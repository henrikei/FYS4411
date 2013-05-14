#include "vmcsolverimportancesampling.h"
#include "lib.h"



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
    double deltaE;

    double dPsidAlpha = 0;
    double dPsidBeta = 0;
    double alphaTerm1 = 0;
    double alphaTerm2 = 0;
    double betaTerm1 = 0;
    double betaTerm2 = 0;

    int acceptCount = 0;

    ofstream outfile;

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
            wf->update(rNew);
            quantumForceNew = wf->getQuantumForceRatio();

            // Check for step acceptance (if yes, update position and quantum force, if no, reset position)
            if(ran2(&idum) <= (getGreensFunctionRatio(rNew, rOld, quantumForceNew, quantumForceOld)*ratio2)){
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

                // write out the radial distances of all particles every 20th cycle
                if(oneBody && (cycle + i) % 20 == 0){
                    if(oneBody < 2){
                        outfile.open("onebody.dat");
                        oneBody = 2;
                    }
                    for (int j = 0; j < nParticles; j++){
                        double radialDist = 0;
                        for (int k = 0; k < nDimensions; k++){
                            radialDist += rNew(j,k)*rNew(j,k);
                        }
                        outfile << sqrt(radialDist) << "  ";
                    }
                    outfile << endl;
                }
            }
        }
    }

    energy = energySum/(nCycles * nParticles);
    variance = energySquaredSum/(nCycles * nParticles) - energy*energy;

    if(oneBody){
        outfile.close();
    }

    if(minimizer ==1){
        alphaTerm1 = alphaTerm1/(nCycles*nParticles);
        alphaTerm2 = alphaTerm2/(nCycles*nParticles);

        betaTerm1 = betaTerm1/(nCycles*nParticles);
        betaTerm2 = betaTerm2/(nCycles*nParticles);

        dEdAlpha = 2*(alphaTerm2 - alphaTerm1*energy);
        dEdBeta = 2*(betaTerm2 - betaTerm1*energy);
    }

    //cout << "Acceptance ratio: " << (double)acceptCount/((double)((nCycles + thermalization)*nParticles)) << endl;
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
