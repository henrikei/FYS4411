#ifndef LOCALENERGY_H
#define LOCALENERGY_H

#include<waveFunction/wavefunction.h>

class localEnergy
{
public:
    localEnergy();
    virtual void calculate(const mat &r, waveFunction*wf, double charge)=0;
    virtual void setR(double r)=0;
    double getKinetic();
    double getPotential();
    double getTotal();
protected:
    double potentialEnergy;
    double kineticEnergy;
    double totalEnergy;

};

#endif // LOCALENERGY_H
