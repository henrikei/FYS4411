#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "vmcsolver/vmcsolver.h"


class Minimizer
{
public:
    Minimizer();
    void run(VMCSolver *solver, waveFunction *wf, double alph, double bet);
    double getAlpha();
    double getBeta();
private:
    double alpha;
    double beta;
    double hAlpha;
    double hBeta;
    double hAlphaMax;
    double hBetaMax;
    double toler;
    int maxIter;
};

#endif // MINIMIZER_H
