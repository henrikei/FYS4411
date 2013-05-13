#ifndef HYDROGENIC_H
#define HYDROGENIC_H

#include "orbitals/orbitals.h"


class Hydrogenic : public Orbitals
{
public:
    Hydrogenic(int nPart);
    void setR(double R);
    double getValue(int particleNum, int quantumNum, const mat &R);
    rowvec3 getGradient(int particleNum, int quantumNum, const mat &R);
    double getLaplacian(int particleNum, int quantumNum, const mat &R);
    double getAlphaDerivative(int particleNum, int quantumNum, const mat &R);
};

#endif // HYDROGENIC_H
