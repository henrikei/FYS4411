#include "vmcsolver/vmcsolver.h"
#include "vmcsolver/vmcsolverbruteforce.h"
#include "vmcsolver/vmcsolverimportancesampling.h"
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"
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
    int nParticles = 4;
    int charge = 4;
    double alpha = 4;
    double beta = 1;
    int jastrow = 1;
    int importanceSampling = 1;

    srand(time(NULL));

    VMCSolver *solver;
    if (importanceSampling == 0){
        solver = new VMCSolverBruteForce();
    } else {
        solver = new VMCSolverImportanceSampling();
    }

    waveFunction *wf = new waveFunction(nParticles, alpha, beta, jastrow);
    localEnergy *localE = new localEnergy();
    solver->setCharge(charge);
    solver->setWaveFunction(wf);
    solver->setLocalEnergy(localE);

    Minimizer minimize;
    minimize.run(solver, wf, alpha, beta);

//    solver->runMonteCarloIntegration();
//    cout << "Energy: " << solver->getEnergy() << endl << "Variance: " << solver->getVariance();

    return 0;
}
