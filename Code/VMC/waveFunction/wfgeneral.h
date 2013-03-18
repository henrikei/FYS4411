#ifndef WFGENERAL_H
#define WFGENERAL_H

#include "wavefunction.h"


class wfGeneral : public waveFunction
{
public:
    wfGeneral(const int &nPart, const double &alph);
    void update(const mat &r);
    double getRatio(const int &particleNum, const mat &r);
    mat getQuantumForceRatio(const mat &r);
    double getLaplaceRatio(const mat &r, const double &h);
private:
    Slater slater;
};

#endif // WFGENERAL_H
