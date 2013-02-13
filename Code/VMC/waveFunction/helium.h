#ifndef HELIUM_H
#define HELIUM_H

#include "wavefunction.h"

class helium: public waveFunction
{
private:
    double alpha;
public:
    helium();
    void setAlpha(const double);
    double getValue(const mat &);
    double getLaplacian(const mat &, const double &);
};

#endif // HELIUM_H
