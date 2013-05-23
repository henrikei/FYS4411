TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

SOURCES += main.cpp \
    lib.cpp \
    waveFunction/wavefunction.cpp \
    localenergy/localEnergy.cpp \
    vmcsolver/vmcsolver.cpp \
    vmcsolver/vmcsolverbruteforce.cpp \
    vmcsolver/vmcsolverimportancesampling.cpp \
    orbitals/orbitals.cpp \
    Slater/slater.cpp \
    Jastrow/jastrow.cpp \
    Jastrow/nojastrow.cpp \
    orbitals/hydrogenic.cpp \
    orbitals/diatomic.cpp \
    localenergy/atomicham.cpp \
    localenergy/diatomicham.cpp \
    minimizer/minimizeralpha.cpp \
    minimizer/minimizeralphabeta.cpp

HEADERS += \
    lib.h \
    waveFunction/wavefunction.h \
    localenergy/localEnergy.h \
    vmcsolver/vmcsolver.h \
    vmcsolver/vmcsolverbruteforce.h \
    vmcsolver/vmcsolverimportancesampling.h \
    orbitals/orbitals.h \
    Slater/slater.h \
    Jastrow/jastrow.h \
    Jastrow/nojastrow.h \
    orbitals/hydrogenic.h \
    orbitals/diatomic.h \
    localenergy/diatomicham.h \
    localenergy/atomicham.h \
    minimizer/minimizeralpha.h \
    minimizer/minimizeralphabeta.h

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

release {
    DEFINES += ARMA_NO_DEBUG
}

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS = $$system(mpicc --showme:compile)
QMAKE_LFLAGS = $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS = $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
