#ifndef JASTROW_H
#define JASTROW_H

#include <armadillo>

using namespace arma;

class Jastrow
{
public:
    Jastrow();
    Jastrow(const int &nPart, const double &b);
    double getRatio(const int &particleNum, const mat &rNew, const mat &rOld);
    mat getQuantumForceRatio(const mat &r);
    double getLaplaceRatio(const mat &r);
private:
    double f(const double &r12, const int &particleNum1, const int &particleNum2);
    double dfdr(const double &r12, const int &particleNum1, const int &particleNum2);
    double d2fdr2(const double &r12, const int &particleNum1, const int &particleNum2);
    int nParticles;
    int nDimensions;
    double beta;

    rowvec3 r12Vec;
    mat quantumForceRatio;
};

#endif // JASTROW_H
