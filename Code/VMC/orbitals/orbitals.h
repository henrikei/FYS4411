#ifndef ORBITALS_H
#define ORBITALS_H

#include <armadillo>

using namespace arma;


class Orbitals
{
public:
    Orbitals();
    void setNumParticles(const int &nPart);
    double getValue(const int &particleNum, const int &quantumNum, const mat &R);
    rowvec3 getGradient(const int &particleNum, const int &quantumNum, const mat &R);
    double getLaplacian(const int &particleNum, const int &quantumNum, const mat &R);
private:
    int nParticles;
    int nDimensions;
    double alpha;
    double value;
    rowvec3 gradient;
    double laplacian;
};

#endif // ORBITALS_H
