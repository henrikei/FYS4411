#ifndef DIATOMIC_H
#define DIATOMIC_H

#include "orbitals/orbitals.h"


class Diatomic : public Orbitals
{
public:
    Diatomic(int nPart, double r, double a);
    double getValue(int particleNum, int quantumNum, const mat &R);
    rowvec3 getGradient(int particleNum, int quantumNum, const mat &R);
    double getLaplacian(int particleNum, int quantumNum, const mat &R);
private:
    rowvec3 rNuclei; //vector from origo to proton1
};

#endif // DIATOMIC_H
