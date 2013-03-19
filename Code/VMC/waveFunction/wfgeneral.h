#ifndef WFGENERAL_H
#define WFGENERAL_H

#include "wavefunction.h"


class wfGeneral : public waveFunction
{
public:
    wfGeneral(const int &nPart, const double &a, const double &b);
    void setAlpha(const double &a);
    void update(const mat &r);
    double getRatio(const int &particleNum, const mat &rNew, const mat &rOld);
    mat getQuantumForceRatio(const mat &r);
    double getLaplaceRatio(const mat &r, const double &h);
private:
    Slater slater;
    Jastrow jastrow;
    mat gradientRatioSD;
    mat gradientRatioJ;
};

#endif // WFGENERAL_H
