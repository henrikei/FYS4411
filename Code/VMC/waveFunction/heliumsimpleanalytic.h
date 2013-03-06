#ifndef HELIUMSIMPLEANALYTIC_H
#define HELIUMSIMPLEANALYTIC_H

#include "wavefunction.h"


class heliumSimpleAnalytic: public waveFunction
{
public:
    heliumSimpleAnalytic();
    double getValue (const mat &r);
    mat getQuantumForce(const mat &r);
    double getLaplacian (const mat &, const double &);
};

#endif // HELIUMSIMPLEANALYTIC_H
