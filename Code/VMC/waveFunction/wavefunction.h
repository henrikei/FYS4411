#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include<armadillo>

using namespace arma;

class waveFunction
{
public:
    waveFunction();
    virtual double getValue(const mat &r)=0;
    virtual double getLaplacian(const mat &r, const double &)=0;
};

#endif // WAVEFUNCTION_H
