#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include "Slater/slater.h"
#include "Jastrow/jastrow.h"
#include "Jastrow/nojastrow.h"

using namespace arma;

class waveFunction
{
public:
    waveFunction(const int &nPart, const double &a, const double &b, const int &jas);
    void setAlpha (const double &a);
    void setBeta (const double &a);
    int getNParticles();
    int getNDimensions();
    void update(const mat &r);
    double getRatio(const int &particleNum, const mat &rNew, const mat &rOld);
    mat getQuantumForceRatio();
    double getLaplaceRatio(const mat &r);
    double getAlphaDerivativeRatio(const mat &r);
    double getBetaDerivativeRatio(const mat &r);
protected:
    double alpha;
    double beta;
    int nParticles;
    int nDimensions;
    Slater slater;
    Jastrow *jastrow;
    mat gradientRatioSD;
    mat gradientRatioJ;
};

#endif // WAVEFUNCTION_H
