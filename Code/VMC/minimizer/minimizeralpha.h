#ifndef MINIMIZERALPHA_H
#define MINIMIZERALPHA_H

#include <vmcsolver/vmcsolver.h>

class MinimizerAlpha
{
public:
    MinimizerAlpha();
    void run(VMCSolver *solver, waveFunction *wf, double alph);
    double getAlpha();
private:
    double alpha;
    double hAlpha;
    double hAlphaMax;
    double toler;
    int maxIter;
};

#endif // MINIMIZERALPHA_H
