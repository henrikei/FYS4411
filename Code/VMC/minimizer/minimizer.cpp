#include "minimizer.h"
#include "waveFunction/wavefunction.h"

Minimizer::Minimizer()
{
    hAlpha = 0.1;
    hBeta = 0.1;
    hAlphaMax = 0.2;
    hBetaMax = 0.2;
    maxIter = 100;
    toler = 0.0001;
}

void Minimizer::run(VMCSolver *solver, waveFunction *wf, double alph, double bet){
    int counter = 0;
    double dAlphaOld = 0;
    double dAlphaNew = 0;
    double dBetaOld = 0;
    double dBetaNew = 0;

    alpha = alph;
    beta = bet;

    while (counter < maxIter){
        counter += 1;

        wf->setAlpha(alpha);
        wf->setBeta(beta);
        solver->runMonteCarloIntegration();

        dAlphaOld = dAlphaNew;
        dBetaOld = dBetaNew;

        dAlphaNew = solver->getdEdAlpha();
        dAlphaNew = dAlphaNew/fabs(dAlphaNew);
        dBetaNew = solver->getdEdBeta();
        dBetaNew = dBetaNew/fabs(dBetaNew);

        if(dAlphaOld/dAlphaNew < 0){
            hAlpha = 0.5*hAlpha;
        } else if(hAlpha < hAlphaMax){
            hAlpha = 1.25*hAlpha;
        }

        if(dBetaOld/dBetaNew < 0){
            hBeta = 0.5*hBeta;
        } else if(hBeta < hBetaMax){
            hBeta = 1.25*hBeta;
        }

        alpha -= hAlpha*dAlphaNew;
        beta -= hBeta*dBetaNew;

        cout << "Alpha = " << alpha << endl;
        cout << "Beta = " << beta << endl;

        if(hAlpha < toler && hBeta < toler){
            cout << "Reached tolerance" << endl;
            break;
        }
    }
}

double Minimizer::getAlpha(){
    return alpha;
}

double Minimizer::getBeta(){
    return beta;
}
