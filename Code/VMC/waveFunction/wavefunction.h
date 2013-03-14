#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <armadillo>

using namespace arma;

class waveFunction
{
public:
    waveFunction();
    void setAlpha (const double &);
    void setBeta (const double &);
    int getNParticles();
    int getNDimensions();
    virtual double getValue(const mat &r)=0;
    virtual mat getQuantumForce(const mat &r);
    virtual double getLaplacian(const mat &r, const double &);
protected:
    double alpha;
    double beta;
    int nParticles;
    int nDimensions;
    mat quantumForce;
};

#endif // WAVEFUNCTION_H
