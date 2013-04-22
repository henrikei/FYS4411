#include "minimizer.h"
#include "waveFunction/wavefunction.h"

Minimizer::Minimizer()
{
    hAlpha = 0.1;
    hBeta = 0.1;
    hAlphaMax = 0.2;
    hBetaMax = 0.2;
    maxIter = 40;
}

void Minimizer::run(VMCSolver *solver, waveFunction *wf, double alph, double bet){
    int counter = 0;
    double dAlphaOld = 0;
    double dAlphaNew = 0;
    double dBetaOld = 0;
    double dBetaNew = 0;

    alpha = alph;
    beta = bet;
    wf->setAlpha(alpha);
    wf->setBeta(beta);

    while (counter < maxIter){
        solver->runMonteCarloIntegration();

        dAlphaOld = dAlphaNew;
        dBetaOld = dBetaNew;

        dAlphaNew = solver->getdEdAlpha();
        dAlphaNew = dAlphaNew/fabs(dAlphaNew);
        dBetaNew = solver->getdEdBeta();
        dBetaNew = dBetaNew/fabs(dBetaNew);

        if(dAlphaOld/dAlphaNew < 0){
            hAlpha = hAlpha/2;
        } else if(hAlpha < hAlphaMax){
            hAlpha = 1.1*hAlpha;
        }

        if(dBetaOld/dBetaNew < 0){
            hBeta = hBeta/2;
        } else if(hBeta < hBetaMax){
            hBeta = 1.1*hBeta;
        }

        alpha -= hAlpha*dAlphaNew;
        beta -= hBeta*dBetaNew;
    }
}
