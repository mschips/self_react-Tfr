DISTFILES += \
    ../bcinflow09rand0.par \
    ../leda-cellrep12.par \
    ../parameter_files/bcinflow02.par \
    ../parameter_files/bcinflow09.par \
    ../parameter_files/bcinflow09rand0.par \
    ../parameter_files/bcinflow09rand0abtest.par \
    ../parameter_files/bcinflow09rand0multiag001.par \
    ../parameter_files/bcinflow12.par \
    ../parameter_files/bcinflow12rand0.par \
    ../parameter_files/bcinflow13.par \
    ../parameter_files/bcinflow14.par \
    ../parameter_files/bcinflow15.par \
    ../parameter_files/bcinflow16.par \
    ../parameter_files/leda-cellrep12.par \
    ../parameter_files/leda-cr12_bcinflow.par \
    ../parameter_files/tmp.par

HEADERS += \
    affinityspace.h \
    cell.h \
    cellman.h \
    cellthis.h \
    dynarray.h \
    grid.h \
    gridpoint.h \
    kinetics.h \
    ode.h \
    odelist.h \
    random.h \
    setparam.h \
    signals.h \
    space.h \
    ss.h \
    track.h \
    brainbow.h \
    antibody.h \
    sequencespace.h

SOURCES += \
    cell.cpp \
    cellman.cpp \
    cellthis.cpp \
    grid.cpp \
    gridpoint.cpp \
    hyphasma.cpp \
    kinetics.cpp \
    ode.cpp \
    odelist.cpp \
    random.cpp \
    setparam.cpp \
    signals.cpp \
    space.cpp \
    ss.cpp \
    track.cpp \
    brainbow.cpp \
    antibody.cpp \
    sequencespace.cpp

QMAKE_CXXFLAGS_WARN_ON += -Wno-unused-parameter

