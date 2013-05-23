#ifndef MINIMIZERALPHABETA_H
#define MINIMIZERALPHABETA_H

#include <vmcsolver/vmcsolver.h>


class MinimizerAlphaBeta
{
public:
    MinimizerAlphaBeta();
    void run(VMCSolver *solver, waveFunction *wf, double alph, double bet);
    double getAlpha();
    double getBeta();
private:
    double alpha;
    double hAlpha;
    double hAlphaMax;
    double beta;
    double hBeta;
    double hBetaMax;
    double toler;
    int maxIter;
};

#endif // MINIMIZERALPHABETA_H
