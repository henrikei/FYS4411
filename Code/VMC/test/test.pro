#-------------------------------------------------
#
# Project created by QtCreator 2013-02-13T18:39:17
#
#-------------------------------------------------

QT       += testlib

QT       -= gui

TARGET = tst_testvmc
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app
INCLUDEPATH = ../


SOURCES += tst_testvmc.cpp \
    ../waveFunction/wavefunction.cpp \
    ../waveFunction/helium.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

HEADERS += \
    ../waveFunction/wavefunction.h \
    ../waveFunction/helium.h
