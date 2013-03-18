#include "vmcsolver/vmcsolver.h"
#include "vmcsolver/vmcsolverbruteforce.h"
#include "vmcsolver/vmcsolverimportancesampling.h"
#include "waveFunction/wavefunction.h"
#include "waveFunction/heliumSimpleNum.h"
#include "waveFunction/heliumsimpleanalytic.h"
#include "waveFunction/heliumJastrowNum.h"
#include "waveFunction/heliumjastrowanalytic.h"
#include "waveFunction/berylliumsimplenum.h"
#include "waveFunction/wfgeneral.h"
#include "localenergy/localEnergy.h"
#include "orbitals/orbitals.h"
#include "Jastrow/jastrow.h"

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


int main()
{
//    Jastrow jastrow(4, 0.8);
//    mat rNew = zeros(4,3);
//    mat rOld = zeros(4,3);
//    for (int i = 0; i < 4; i++){
//        for (int j = 0; j < 3; j++){
//            rNew(i,j) = 0.3 +i*j - j + 2*i;
//            rOld(i,j) = 1.3 +0.2*i*j;
//        }
//    }
//    cout << "rNew :" << rNew << endl;
//    cout << "rOld :" << rOld << endl;
//    jastrow.getQuantumForceRatio(rNew);
//    cout << "quantumForceRatio: " << jastrow.getQuantumForceRatio(rNew) << endl;
//    cout << "laplaceRatio: " << jastrow.getLaplaceRatio(rNew) << endl;


    VMCSolver *solver = new VMCSolverImportanceSampling(4);
    waveFunction *wf = new wfGeneral(4, 3.925, 0.109);
    localEnergy localE;
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
