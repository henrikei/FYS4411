#include "vmcsolver.h"
#include "waveFunction/wavefunction.h"
#include "waveFunction/helium.h"
#include "waveFunction/heliumwithjastrow.h"
#include "localenergy/localEnergy.h"

#include <iostream>
#include <armadillo>

using namespace arma;
using namespace std;


int main()
{   
    helium *wfhelium = new helium;
    localEnergy localE;
    wfhelium->setAlpha(1.8);
    waveFunction *wf = wfhelium;
    VMCSolver *solver = new VMCSolver;
    solver->setWaveFunction(wf);
    solver->setLocalEnergy(localE);
    solver->runMonteCarloIntegration();
    cout << solver->getEnergy();
    /*heliumWithJastrow *wf1 = new heliumWithJastrow;
    localEnergy localE;
    double alphamin = 1.5;
    double alphamax = 2;
    double betamin = 0;
    double betamax = 0.5;
    int n = 10;

    ofstream ofile;
    ofile.open("results.dat");
    double alpha;
    double beta;
    double deltaAlpha = (double) (alphamax - alphamin)/n;
    double deltaBeta = (double) (betamax - betamin)/n;
    for (int i = 0; i < n+1; i++){
        for (int j = 0; j < n+1; j++){
            alpha = alphamin + deltaAlpha*i;
            beta = betamin + deltaBeta*j;
            wf1->setAlpha(alpha);
            wf1->setBeta(beta);
            waveFunction *wf = wf1;
            VMCSolver *solver = new VMCSolver();
            solver->setWaveFunction(wf);
            solver->setLocalEnergy(localE);
            solver->runMonteCarloIntegration();
            ofile << alpha << "  " << beta << "  " << solver->getEnergy() << "  " << solver->getEnergy() << endl;
        }
    }
    ofile.close();*/
    return 0;
}
