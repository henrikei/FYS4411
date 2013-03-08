#ifndef ORBITALS_H
#define ORBITALS_H

#include <armadillo>

using namespace arma;


class Orbitals
{
public:
    Orbitals();
    double getValue(const int &quantumNum, const int &particleNum, double const &alpha, const mat &R);
    vec getGradient(const int &quantumNum, const int &particleNum, double const &alpha, const mat &R);
    mat getLaplacian(const int &quantumNum, const int &particleNum, double const &alpha, const mat &R);
private:
    int nParticles;
    int nDimensions;
    double alpha;
    double value;
    vec gradient;
    double laplacian;
};

#endif // ORBITALS_H
