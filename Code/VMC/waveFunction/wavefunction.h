#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>
#include "Slater/slater.h"

using namespace arma;

class waveFunction
{
public:
    waveFunction();
    void setAlpha (const double &a);
    void setBeta (const double &a);
    int getNParticles();
    int getNDimensions();
    virtual void update(const mat &r);
    virtual double getValue(const mat &r);
    virtual double getRatio(const int &particleNum, const mat &r);
    virtual mat getQuantumForceRatio(const mat &r);
    virtual double getLaplacian(const mat &r, const double &h);
    virtual double getLaplaceRatio(const mat &r, const double &h);
protected:
    double alpha;
    double beta;
    int nParticles;
    int nDimensions;
    mat quantumForce;
};

#endif // WAVEFUNCTION_H
