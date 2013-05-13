#ifndef LOCALENERGY_H
#define LOCALENERGY_H

#include<waveFunction/wavefunction.h>

class localEnergy
{
public:
    localEnergy();
    virtual double getValue(const mat &r, waveFunction*wf, double charge)=0;
    virtual void setR(double r)=0;
};

#endif // LOCALENERGY_H
