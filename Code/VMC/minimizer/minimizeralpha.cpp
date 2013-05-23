#include "minimizeralpha.h"

MinimizerAlpha::MinimizerAlpha()
{
    hAlpha = 0.01;
    hAlphaMax = 0.1;
    maxIter = 100;
    toler = 0.0001;
}

void MinimizerAlpha::run(VMCSolver *solver, waveFunction *wf, double alph){
    int counter = 0;
    double dAlphaOld = 0;
    double dAlphaNew = 0;

    alpha = alph;

    while (counter < maxIter){
        counter += 1;

        wf->setAlpha(alpha);
        solver->runMonteCarloIntegration();

        dAlphaOld = dAlphaNew;

        dAlphaNew = solver->getdEdAlpha();
        dAlphaNew = dAlphaNew/fabs(dAlphaNew);

        if(dAlphaOld/dAlphaNew < 0){
            hAlpha = 0.5*hAlpha;
        } else if(hAlpha < hAlphaMax){
            hAlpha = 1.25*hAlpha;
        }

        alpha -= hAlpha*dAlphaNew;

//        if (solver->getMPIRank() ==0 ){
//            cout << "Alpha = " << alpha << endl;
//        }
        if(hAlpha < toler){
            break;
        }
    }
}

double MinimizerAlpha::getAlpha(){
    return alpha;
}
