#ifndef HELIUMWITHJASTROW_H
#define HELIUMWITHJASTROW_H

#include "wavefunction.h"

class heliumJastrowNum: public waveFunction
{
public:
    heliumJastrowNum();
    double getValue (const mat &);
};

#endif // HELIUMWITHJASTROW_H
