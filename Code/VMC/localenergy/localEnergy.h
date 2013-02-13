#ifndef LOCALENERGY_H
#define LOCALENERGY_H

#include<waveFunction/wavefunction.h>

class localEnergy
{
public:
    localEnergy();
    double getValue(const mat &, waveFunction*, const double &, const double &);
};

#endif // LOCALENERGY_H
