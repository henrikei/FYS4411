#ifndef ORBITALS_H
#define ORBITALS_H

#include <armadillo>

using namespace arma;


class Orbitals
{
public:
    Orbitals();
    virtual void setAlpha(double a);
    virtual void setR(double R)=0;
    virtual double getValue(int particleNum, int quantumNum, const mat &R)=0;
    virtual rowvec3 getGradient(int particleNum, int quantumNum, const mat &R)=0;
    virtual double getLaplacian(int particleNum, int quantumNum, const mat &R)=0;
    virtual double getAlphaDerivative(int particleNum, int quantumNum, const mat &R)=0;
protected:
    int nParticles;
    int nDimensions;
    double alpha;
    rowvec3 gradient;
};

#endif // ORBITALS_H
