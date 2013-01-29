#include "vmcsolver.h"

#include <iostream>

using namespace std;
/*
MINI::MINI(VMCSolver* vmc){
    this->vmcObject = vmc
    ...
}

void MINI::loopALPHA(){
    FOR ....
            vmcObject->setAlpha(...)
            vmcObject->solve()

}
*/

int main()
{
    VMCSolver *solver = new VMCSolver(1.844,0.36);
    solver->runMonteCarloIntegration();
    return 0;
}

