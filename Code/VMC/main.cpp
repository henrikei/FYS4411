#include "vmcsolver/vmcsolver.h"
#include "vmcsolver/vmcsolverbruteforce.h"
#include "vmcsolver/vmcsolverimportancesampling.h"
#include "waveFunction/wavefunction.h"
#include "localenergy/localEnergy.h"
#include "localenergy/atomicham.h"
#include "localenergy/diatomicham.h"
#include "orbitals/orbitals.h"
#include "Jastrow/jastrow.h"
#include "minimizer/minimizeralpha.h"
#include "minimizer/minimizeralphabeta.h"

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
    int nParticles = 8;
    int charge = 4;
    double alpha = 3.766;
    double beta = 0.521;
    double R = 2.0;
    int jastrow = 1;
    int importanceSampling = 1;
    int nCycles = 100000; //total nCycles
    int minimize = 0;
    int oneBody = 1;
    string orbitalType = "Diatomic";          // Hydrogenic or Diatomic
    string hamiltonianType = "DiatomicHam";     // AtomicHam or DiatomicHam



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
        if (jastrow == 0){
            MinimizerAlpha minimizer;
            minimizer.run(solver, wf, alpha);
            alpha = minimizer.getAlpha();

            if (my_rank == 0) {
                cout << "alpha = " << alpha << endl;
            }
        } else {
            MinimizerAlphaBeta minimizer;
            minimizer.run(solver, wf, alpha, beta);
            alpha = minimizer.getAlpha();
            beta = minimizer.getBeta();

            if (my_rank == 0) {
                cout << "alpha = " << alpha << ", beta = " << beta << endl;
            }
        }
    }




//    // Run calculation
//    ofstream ofile;
//    ofile.open("../Out/Hydrogen_Energy.dat");
//    int N = 100;
//    int nCycles1 = 10000;
//    int nCycles2 = 10000;
//    double Rmin = 0.2;
//    double Rmax = 8.0;
//    alpha = 4.0;
//    beta = 0.5;
//    double deltaR = (Rmax - Rmin)/((double)N);
//    for (int i = 0; i < N; i++){
//        R = Rmin + i*deltaR;
//        wf->setR(R);
//        localE->setR(R);
//        if (minimize == 1){
//            solver->calcEnergyGradients();
//            solver->setNumOfCycles(nCycles1);
//            MinimizerAlphaBeta minimizer;
//            minimizer.run(solver, wf, alpha, beta);
//            alpha = minimizer.getAlpha();
//            beta = minimizer.getBeta();
//            cout << "alpha = " << alpha << ", beta = " << beta;
//        }

//        solver->setNumOfCycles(nCycles2);
//        solver->runMonteCarloIntegration();
//        if (my_rank == 0){
//            ofile << R << "  " << solver->getEnergy() << "  " << solver->getPotentialEnergy() << endl;
//            cout << ", R = " << R << ", E = " << solver->getEnergy() << ", EP = "<< solver->getPotentialEnergy() <<  endl;
//        }
//    }

//    ofile.close();
    solver->runMonteCarloIntegration();
    if (my_rank == 0){
        cout << "Energy: " << solver->getEnergy() << endl << "Variance (no blocking) = " << solver->getVariance() << endl;
    }

    MPI_Finalize();
    return 0;
}
