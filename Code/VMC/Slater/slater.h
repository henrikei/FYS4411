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
    Slater(string orbitalType, int nPart, double alph);
    void setAlpha(const double &a);
    void update(const mat &R);
    double getRatio(const int &particleNum, const mat &R);
    mat getGradientRatio(const mat &R);
    double getLaplaceRatio(const mat &R);
    double getAlphaDerivativeRatio(const mat &R);
private:
    int nParticles;
    int nDimensions;
    double alpha;
    Orbitals *orbitals;

    mat slaterUp;
    mat slaterDown;
    mat invSlaterUp;
    mat invSlaterDown;

    mat gradientRatio;
};

#endif // SLATER_H
