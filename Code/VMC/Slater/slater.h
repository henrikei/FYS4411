#ifndef SLATER_H
#define SLATER_H

#include <armadillo>
#include "orbitals/orbitals.h"
#include "orbitals/hydrogenic.h"
#include "orbitals/diatomic.h"

using namespace std;
using namespace arma;


class Slater
{
public:
    Slater();
    Slater(string orbitalType, int nPart);
    void setAlpha(double alpha);
    void setR(double R);
    void update(const mat &R);
    double getRatio(const int &particleNum, const mat &R);
    mat getGradientRatio(const mat &R);
    double getLaplaceRatio(const mat &R);
    double getAlphaDerivativeRatio(const mat &R);
private:
    int nParticles;
    int nDimensions;
    Orbitals *orbitals;

    mat slaterUp;
    mat slaterDown;
    mat invSlaterUp;
    mat invSlaterDown;

    mat gradientRatio;
};

#endif // SLATER_H
