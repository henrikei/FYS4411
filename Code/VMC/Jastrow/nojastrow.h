#ifndef NOJASTROW_H
#define NOJASTROW_H

#include "jastrow.h"


class NoJastrow : public Jastrow
{
public:
    NoJastrow(const int &nPart);
    double getRatio(const int &particleNum, const mat &rNew, const mat &rOld);
    mat getGradientRatio(const mat &r);
    double getLaplaceRatio(const mat &r);
};

#endif // NOJASTROW_H
