TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

SOURCES += main.cpp \
    vmcsolver.cpp \
    lib.cpp \
    waveFunction/wavefunction.cpp \
    waveFunction/helium.cpp \
    waveFunction/heliumwithjastrow.cpp \
    localenergy/localEnergy.cpp

HEADERS += \
    vmcsolver.h \
    lib.h \
    waveFunction/wavefunction.h \
    localenergy/localEnergy.h \
    waveFunction/heliumWithJastrow.h \
    waveFunction/helium.h

