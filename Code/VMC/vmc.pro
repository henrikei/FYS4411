TEMPLATE = app
CONFIG += console
CONFIG -= qt

LIBS += -llapack -larmadillo

SOURCES += main.cpp \
    lib.cpp \
    waveFunction/wavefunction.cpp \
    localenergy/localEnergy.cpp \
    waveFunction/heliumSimpleNum.cpp \
    waveFunction/heliumJastrowNum.cpp \
    waveFunction/heliumsimpleanalytic.cpp \
    waveFunction/heliumjastrowanalytic.cpp \
    vmcsolver/vmcsolver.cpp \
    vmcsolver/vmcsolverbruteforce.cpp \
    vmcsolver/vmcsolverimportancesampling.cpp

HEADERS += \
    lib.h \
    waveFunction/wavefunction.h \
    localenergy/localEnergy.h \
    waveFunction/heliumJastrowNum.h \
    waveFunction/heliumSimpleNum.h \
    waveFunction/heliumsimpleanalytic.h \
    waveFunction/heliumjastrowanalytic.h \
    vmcsolver/vmcsolver.h \
    vmcsolver/vmcsolverbruteforce.h \
    vmcsolver/vmcsolverimportancesampling.h

QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3

release {
    DEFINES += ARMA_NO_DEBUG
}
