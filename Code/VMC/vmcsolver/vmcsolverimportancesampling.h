#ifndef VMCSOLVERIMPORTANCESAMPLING_H
#define VMCSOLVERIMPORTANCESAMPLING_H

#include "vmcsolver.h"


class VMCSolverImportanceSampling : public VMCSolver
{
public:
    VMCSolverImportanceSampling(const int &charg);
    void runMonteCarloIntegration();
    double getGreensFunctionRatio(const mat &, const mat &, const mat &, const mat &);
private:
    double timeStep;
};

#endif // VMCSOLVERIMPORTANCESAMPLING_H
