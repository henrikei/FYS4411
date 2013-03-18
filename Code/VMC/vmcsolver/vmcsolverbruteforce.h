#ifndef VMCSOLVERBRUTEFORCE_H
#define VMCSOLVERBRUTEFORCE_H

#include "vmcsolver.h"

using namespace arma;

class VMCSolverBruteForce : public VMCSolver
{
public:
    VMCSolverBruteForce(const int &charg);
    void runMonteCarloIntegration();
private:
    double stepLength;
    int nDummyCycles;
};

#endif // VMCSOLVER_H
