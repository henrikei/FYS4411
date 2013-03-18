#include "jastrow.h"

Jastrow::Jastrow()
{
}

Jastrow::Jastrow(const int &nPart, const double &b)
{
    nParticles = nPart;
    nDimensions = 3;
    beta = b;
}

double Jastrow::getRatio(const int &particleNum, const mat &rNew, const mat &rOld){
    double argument = 0;
    double fNew = 0;
    double fOld = 0;
    double r12 = 0;
    for (int i = 0; i < particleNum; i++){
        r12Vec = rNew.row(particleNum) - rNew.row(i);
        r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
        fNew = f(r12, particleNum, i);

        r12Vec = rOld.row(particleNum) - rOld.row(i);
        r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
        fOld = f(r12, particleNum, i);

        argument += fNew - fOld;
    }
    for (int i = particleNum + 1; i < nParticles; i++){
        r12Vec = rNew.row(particleNum) - rNew.row(i);
        r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
        fNew = f(r12, particleNum, i);

        r12Vec = rOld.row(particleNum) - rOld.row(i);
        r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
        fOld = f(r12, particleNum, i);

        argument += fNew - fOld;
    }
    return exp(argument);
}

mat Jastrow::getQuantumForceRatio(const mat &r){
    quantumForceRatio = zeros(nParticles, nDimensions);
    for (int k = 0; k < nParticles; k++){
        for (int i = 0; i < k; i++){
            r12Vec = r.row(k) - r.row(i);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            quantumForceRatio.row(k) += r12Vec*dfdr(r12, k, i)/r12;
        }
        for (int i = k + 1; i < nParticles; i++){
            r12Vec = r.row(i) - r.row(k);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            quantumForceRatio.row(k) -= r12Vec*dfdr(r12, k, i)/r12;
        }
    }
    return quantumForceRatio;
}

double Jastrow::getLaplaceRatio(const mat &r){
    double laplaceRatio = 0;
    for (int k = 0; k < nParticles; k++){
        for (int i = 0; i < k; i++){
            r12Vec = r.row(k) - r.row(i);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            laplaceRatio += (nDimensions - 1)*dfdr(r12, k, i)/r12 + d2fdr2(r12, k, i);
        }
        for (int i = k + 1; i < nParticles; i++){
            r12Vec = r.row(k) - r.row(i);
            double r12 = sqrt(r12Vec(0)*r12Vec(0) + r12Vec(1)*r12Vec(1) + r12Vec(2)*r12Vec(2));
            laplaceRatio += (nDimensions - 1)*dfdr(r12, k, i)/r12 + d2fdr2(r12, k, i);
        }
        laplaceRatio += quantumForceRatio(k,0)*quantumForceRatio(k,0) + quantumForceRatio(k,1)*quantumForceRatio(k,1)
                        + quantumForceRatio(k,2)*quantumForceRatio(k,2);
    }
    return laplaceRatio;
}

double Jastrow::f(const double &r12, const int &particleNum1, const int &particleNum2){
    double a = 0;
    if (((particleNum1 < nParticles/2) && (particleNum2 < nParticles/2)) || ((particleNum1 >= nParticles/2) && (particleNum2 >= nParticles/2))){
        a = 0.25;
    } else {
        a = 0.50;
    }
    return a*r12/(1 + beta*r12);
}

double Jastrow::dfdr(const double &r12, const int &particleNum1, const int &particleNum2){
    double a = 0;
    if (((particleNum1 < nParticles/2) && (particleNum2 < nParticles/2)) || ((particleNum1 >= nParticles/2) && (particleNum2 >= nParticles/2))){
        a = 0.25;
    } else {
        a = 0.50;
    }
    return a/((1 + beta*r12)*(1 + beta*r12));
}

double Jastrow::d2fdr2(const double &r12, const int &particleNum1, const int &particleNum2){
    double a = 0;
    if (((particleNum1 < nParticles/2) && (particleNum2 < nParticles/2)) || ((particleNum1 >= nParticles/2) && (particleNum2 >= nParticles/2))){
        a = 0.25;
    } else {
        a = 0.50;
    }
    return -2*a*beta/((1 + beta*r12)*(1 + beta*r12)*(1 + beta*r12));
}
