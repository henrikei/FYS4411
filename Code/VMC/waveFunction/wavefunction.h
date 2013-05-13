#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include "Slater/slater.h"
#include "Jastrow/jastrow.h"
#include "Jastrow/nojastrow.h"

using namespace arma;
using namespace std;


class waveFunction
{
public:
    waveFunction(string orbitalType,int nPart, int jas);
    void setAlpha (double alpha);
    void setBeta (double beta);
    void setR(double R);
    int getNParticles();
    int getNDimensions();
    void update(const mat &r);
    double getRatio(int particleNum, const mat &rNew, const mat &rOld);
    mat getQuantumForceRatio();
    double getLaplaceRatio(const mat &r);
    double getAlphaDerivativeRatio(const mat &r);
    double getBetaDerivativeRatio(const mat &r);
protected:
    int nParticles;
    int nDimensions;
    Slater slater;
    Jastrow *jastrow;
    mat gradientRatioSD;
    mat gradientRatioJ;
};

#endif // WAVEFUNCTION_H
