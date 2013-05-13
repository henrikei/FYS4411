#include "vmcsolver/vmcsolver.h"
#include "vmcsolver/vmcsolverbruteforce.h"
#include "vmcsolver/vmcsolverimportancesampling.h"
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"
#include "localenergy/atomicham.h"
#include "localenergy/diatomicham.h"
#include "orbitals/orbitals.h"
#include "Jastrow/jastrow.h"
#include "minimizer/minimizer.h"

#include <iostream>
#include <sys/time.h>
#include <armadillo>

using namespace arma;
using namespace std;


int main()
{
    // Configuration
    int nParticles = 2;
    int charge = 1;
    double alpha = 1.285;
    double beta = 0.28;
    double R = 1.4;
    int jastrow = 1;
    int importanceSampling = 1;
    int minimize = 0;
    int oneBody = 0;
    string orbitalType = "Diatomic";
    string hamiltonianType = "DiatomicHam";



    // Initialization
    VMCSolver *solver;
    if (importanceSampling == 0){
        solver = new VMCSolverBruteForce();
    } else {
        solver = new VMCSolverImportanceSampling();
    }
    if (oneBody){
        solver->calcOneBodyDensity();
    }

    waveFunction *wf = new waveFunction(orbitalType, nParticles, jastrow);
    wf->setAlpha(alpha);
    wf->setBeta(beta);

    localEnergy *localE;
    if (hamiltonianType == "Atomic"){
        localE = new AtomicHam();
    } else if (hamiltonianType == "DiatomicHam"){
        localE = new DiatomicHam();
        wf->setR(R);
        localE->setR(R);
    } else {
        cout << "Error: Hamiltonian type not defined" << endl;
        exit(1);
    }

    solver->setCharge(charge);
    solver->setWaveFunction(wf);
    solver->setLocalEnergy(localE);


    // Run calculation
    if (minimize == 1){
        solver->calcEnergyGradients();
        Minimizer minimize;
        minimize.run(solver, wf, alpha, beta);
    }

    solver->runMonteCarloIntegration();
    cout << "Energy: " << solver->getEnergy() << endl << "Variance: " << solver->getVariance();

    return 0;
}
