#include "localEnergy.h"
#include "waveFunction/wavefunction.h"

localEnergy::localEnergy()
{
}

double localEnergy::getKinetic(){
    return kineticEnergy;
}

double localEnergy::getPotential(){
    return potentialEnergy;
}

double localEnergy::getTotal(){
    return totalEnergy;
}
