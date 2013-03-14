#ifndef VMCSOLVERIMPORTANCESAMPLING_H
#define VMCSOLVERIMPORTANCESAMPLING_H

#include "vmcsolver.h"


class VMCSolverImportanceSampling : public VMCSolver
{
public:
    VMCSolverImportanceSampling();
    void runMonteCarloIntegration();
    double getGreensFunctionRatio(const mat &, const mat &, const mat &, const mat &, const double &);
private:
    double timeStep;
};

#endif // VMCSOLVERIMPORTANCESAMPLING_H
