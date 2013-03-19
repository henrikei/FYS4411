#ifndef JASTROW_H
#define JASTROW_H

#include <armadillo>

using namespace arma;

class Jastrow
{
public:
    Jastrow();
    Jastrow(const int &nPart, const double &b);
    void setBeta(const double &b);
    virtual double getRatio(const int &particleNum, const mat &rNew, const mat &rOld);
    virtual mat getGradientRatio(const mat &r);
    virtual double getLaplaceRatio(const mat &r);
private:
    double f(const double &r12, const int &particleNum1, const int &particleNum2);
    double dfdr(const double &r12, const int &particleNum1, const int &particleNum2);
    double d2fdr2(const double &r12, const int &particleNum1, const int &particleNum2);
    double beta;

    rowvec3 r12Vec;
protected:
    mat gradientRatio;
    int nParticles;
    int nDimensions;
};

#endif // JASTROW_H
