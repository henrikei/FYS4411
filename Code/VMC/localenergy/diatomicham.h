#ifndef DIATOMICHAM_H
#define DIATOMICHAM_H

#include "localenergy/localEnergy.h"


class DiatomicHam : public localEnergy
{
public:
    DiatomicHam();
    double getValue(const mat &r, waveFunction *wf, double charge);
    void setNucleiDistance(double r);
private:
    rowvec3 rNuclei;
};

#endif // DIATOMIC_H
