#include "vmcsolver/vmcsolver.h"
#include "vmcsolver/vmcsolverbruteforce.h"
#include "vmcsolver/vmcsolverimportancesampling.h"
#include "waveFunction/wavefunction.h"
#include "waveFunction/heliumSimpleNum.h"
#include "waveFunction/heliumsimpleanalytic.h"
#include "waveFunction/heliumJastrowNum.h"
#include "waveFunction/heliumjastrowanalytic.h"
#include "localenergy/localEnergy.h"

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


int main()
{   
    /*heliumJastrowNum *wfhelium = new heliumJastrowNum;
    localEnergy localE;
    wfhelium->setAlpha(1.844);
    wfhelium->setBeta(0.36);
    waveFunction *wf = wfhelium;
    VMCSolver *solver = new VMCSolver;
    solver->setWaveFunction(wf);
    solver->setLocalEnergy(localE);
    solver->runMonteCarloIntegration();
    cout << solver->getEnergy();*/
    VMCSolver *solver = new VMCSolverImportanceSampling();
    heliumJastrowAnalytic *wf = new heliumJastrowAnalytic;
    localEnergy localE;
    solver->setLocalEnergy(localE);
    solver->setWaveFunction(wf);
    double alphamin = 1;
    double alphamax = 3;
    double betamin = 0;
    double betamax = 1;
    int n = 40;

    ofstream ofile;
    ofile.open("results.dat");
    double alpha;
    double beta;
    double deltaAlpha = (double) (alphamax - alphamin)/n;
    double deltaBeta = (double) (betamax - betamin)/n;
    for (int i = 0; i < n + 1; i++){
        for (int j = 0; j < n+1; j++){
            alpha = alphamin + deltaAlpha*i;
            beta = betamin + deltaBeta*j;
            wf->setAlpha(alpha);
            wf->setBeta(beta);
            solver->runMonteCarloIntegration();
            ofile << solver->getEnergy() << "  ";
            cout << "j: " << j << endl;
        }
        ofile << endl;
    }
    ofile.close();
    return 0;
}
