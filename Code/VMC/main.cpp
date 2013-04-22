#include "vmcsolver/vmcsolver.h"
#include "vmcsolver/vmcsolverbruteforce.h"
#include "vmcsolver/vmcsolverimportancesampling.h"
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"
#include "orbitals/orbitals.h"
#include "Jastrow/jastrow.h"

#include <iostream>
#include <sys/time.h>
#include <armadillo>

using namespace arma;
using namespace std;


int main()
{
    // Configuration
    int nParticles = 2;
    int charge = 2;
    double alpha = 2;
    double beta = 0.5;
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
    solver->runMonteCarloIntegration();
    cout << "Energy: " << solver->getEnergy() << endl << "Variance: " << solver->getVariance();


//    double alphamin = 1;
//    double alphamax = 3;
//    double betamin = 0;
//    double betamax = 1;
//    int n = 10;

//    ofstream ofile;
//    ofile.open("results.dat");
//    double alpha;
//    double beta;
//    double deltaAlpha = (double) (alphamax - alphamin)/n;
//    double deltaBeta = (double) (betamax - betamin)/n;
//    for (int i = 0; i < n + 1; i++){
//        for (int j = 0; j < n+1; j++){
//            alpha = alphamin + deltaAlpha*i;
//            beta = betamin + deltaBeta*j;
//            wf->setAlpha(alpha);
//            wf->setBeta(beta);
//            solver->runMonteCarloIntegration();
//            ofile << solver->getEnergy() << "  ";
//            cout << "j: " << j << endl;
//        }
//        ofile << endl;
//    }
//    ofile.close();
    return 0;
}
