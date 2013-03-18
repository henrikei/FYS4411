#ifndef HELIUMJASTROWANALYTIC_H
#define HELIUMJASTROWANALYTIC_H

#include "waveFunction/wavefunction.h"

class heliumJastrowAnalytic :public waveFunction
{
public:
    heliumJastrowAnalytic();
    double getValue (const mat &);
    mat getQuantumForceRatio (const mat &r);
    double getLaplacian (const mat &, const double &);
};

#endif // HELIUMJASTROWANALYTIC_H
