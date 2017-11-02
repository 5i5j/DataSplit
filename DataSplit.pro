TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    datapartition.cpp \
    pseudoinverse.cpp

HEADERS += \
    datapartition.h \
    pseudoinverse.h
