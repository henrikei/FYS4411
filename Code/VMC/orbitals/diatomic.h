#ifndef DIATOMIC_H
#define DIATOMIC_H

#include "orbitals/orbitals.h"


class Diatomic : public Orbitals
{
public:
    Diatomic(int nPart);
    void setAlpha(double a);
    void setR(double R);
    double getValue(int particleNum, int quantumNum, const mat &R);
    rowvec3 getGradient(int particleNum, int quantumNum, const mat &R);
    double getLaplacian(int particleNum, int quantumNum, const mat &R);
    double getAlphaDerivative(int particleNum, int quantumNum, const mat &R);
private:
    mat rNuclei; //vector from origo to proton1
    Orbitals *hydrogenic;
};

#endif // DIATOMIC_H
