#ifndef SLATER_H
#define SLATER_H

#include <armadillo>
#include "orbitals/orbitals.h"


class Slater
{
public:
    Slater();
    void setNParticles(const int &n);
    void update(const mat &R);
    double getRatio(const int &particleNum, const mat &R);
    mat getQuantumForceRatio(const int &particleNum, const mat &R);
    double getLaplaceRatio(const mat &R);
private:
    int nParticles;
    int nDimensions;
    Orbitals orbitals;

    mat slaterUp;
    mat slaterDown;
    mat invSlaterUp;
    mat invSlaterDown;

    mat quantumForceRatio;
};

#endif // SLATER_H
