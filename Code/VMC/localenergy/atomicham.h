#ifndef ATOMIC_H
#define ATOMIC_H

#include "localenergy/localEnergy.h"

class AtomicHam : public localEnergy
{
public:
    AtomicHam();
    void calculate(const mat &r, waveFunction *wf, double charge);
    void setR(double r);
};

#endif // ATOMIC_H
