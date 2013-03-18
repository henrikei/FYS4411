#ifndef HELIUM_H
#define HELIUM_H

#include "wavefunction.h"

class heliumSimpleNum: public waveFunction
{
public:
    heliumSimpleNum();
    double getValue(const mat &);
};

#endif // HELIUM_H
