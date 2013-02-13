#ifndef HELIUMWITHJASTROW_H
#define HELIUMWITHJASTROW_H

#include "wavefunction.h"

class heliumWithJastrow: public waveFunction
{
private:
    double alpha, beta;
public:
    heliumWithJastrow();
    void setAlpha (const double &);
    void setBeta (const double &);
    double getValue (const mat &);
};

#endif // HELIUMWITHJASTROW_H
