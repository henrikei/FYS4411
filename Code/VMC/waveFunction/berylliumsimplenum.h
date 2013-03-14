#ifndef BERYLLIUMSIMPLENUM_H
#define BERYLLIUMSIMPLENUM_H

#include <vector>
#include "wavefunction.h"

using namespace std;


class BerylliumSimpleNum : public waveFunction
{
public:
    BerylliumSimpleNum();
    double getValue(const mat &r);
private:
    vector<double> args;
    double psi1S(const double &rSingleParticle);
    double psi2S(const double &rSingleParticle);
};

#endif // BERYLLIUMSIMPLENUM_H
