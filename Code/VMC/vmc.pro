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
    Jastrow/nojastrow.cpp

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
    Jastrow/nojastrow.h

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

release {
    DEFINES += ARMA_NO_DEBUG
}
