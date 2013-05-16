#include "vmcsolver/vmcsolver.h"
#include "vmcsolver/vmcsolverbruteforce.h"
#include "vmcsolver/vmcsolverimportancesampling.h"
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"
#include "localenergy/atomicham.h"
#include "localenergy/diatomicham.h"
#include "orbitals/orbitals.h"
#include "Jastrow/jastrow.h"
#include "minimizer/minimizer.h"

#include <iostream>
#include <sys/time.h>
#include <armadillo>
#include <mpi.h>

using namespace arma;
using namespace std;


int main()
{
    // Parallellization initialization
    int numprocs, my_rank;
    MPI_Init (NULL, NULL);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    srand(time(NULL) - my_rank);


    // Configuration
    int nParticles = 10;
    int charge = 10;
    double alpha = 10.267;
    double beta = 0.086;
    double R = 4;
    int jastrow = 1;
    int importanceSampling = 1;
    int nCycles = 1000000; //total nCycles
    int minimize = 0;
    int oneBody = 1;
    string orbitalType = "Hydrogenic";          // Hydrogenic or Diatomic
    string hamiltonianType = "AtomicHam";     // AtomicHam or DiatomicHam



    // Initialization
    VMCSolver *solver;
    if (importanceSampling == 0){
        solver = new VMCSolverBruteForce();
    } else {
        solver = new VMCSolverImportanceSampling();
    }

    solver->init_MPI(my_rank, numprocs);

    if (oneBody){
        solver->calcOneBodyDensity();
    }

    waveFunction *wf = new waveFunction(orbitalType, nParticles, jastrow);
    wf->setAlpha(alpha);
    wf->setBeta(beta);

    localEnergy *localE;
    if (hamiltonianType == "AtomicHam"){
        localE = new AtomicHam();
    } else if (hamiltonianType == "DiatomicHam"){
        localE = new DiatomicHam();
        wf->setR(R);
        localE->setR(R);
    } else {
        cout << "Error: Hamiltonian type not defined" << endl;
        exit(1);
    }

    solver->setCharge(charge);
    solver->setWaveFunction(wf);
    solver->setLocalEnergy(localE);
    solver->setNumOfCycles(nCycles);

    if (minimize == 1){
        solver->calcEnergyGradients();
        Minimizer minimizer;
        minimizer.run(solver, wf, alpha, beta);
        alpha = minimizer.getAlpha();
        beta = minimizer.getBeta();

        if (my_rank == 0) {
            cout << "alpha = " << alpha << ", beta = " << beta << endl;
        }
    }




    // Run calculation

//    ofstream ofile;
//    ofile.open("Hydrogen_Energy.dat");
//    int N = 100;
//    int nCycles1 = 10000;
//    int nCycles2 = 1000000;
//    double Rmin = 0.5;
//    double Rmax = 4.0;
//    alpha = 1.5;
//    beta = 0.5;
//    double deltaR = (Rmax - Rmin)/((double)N);
//    for (int i = 0; i < N; i++){
//        R = Rmin + i*deltaR;
//        wf->setR(R);
//        localE->setR(R);
//        if (minimize == 1){
//            solver->calcEnergyGradients();
//            solver->setNumOfCycles(nCycles1);
//            Minimizer minimizer;
//            minimizer.run(solver, wf, alpha, beta);
//            alpha = minimizer.getAlpha();
//            beta = minimizer.getBeta();
//            cout << "alpha = " << alpha << ", beta = " << beta;
//        }

//        solver->setNumOfCycles(nCycles2);
//        solver->runMonteCarloIntegration();
//        ofile << R << "  " << solver->getEnergy() << endl;

//        cout << ", R = " << R << ", E = " << solver->getEnergy() << endl;
//    }
//    ofile.close();
    solver->runMonteCarloIntegration();
    if (my_rank == 0){
        cout << "Energy: " << solver->getEnergy() << endl << "Variance (one proc) = " << solver->getVariance() << endl;
    }

    MPI_Finalize();
    return 0;
}
