PRJROOT := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
# Path for objects
export BUILDPATH := $(PRJROOT)build
# Path for resulting binaries
export BINPATH := $(PRJROOT)bin
# Name for the resulting binary
export BINARY := hyphasma
# Extension for files to be automatically included
export SRCEXT := cpp
# Path to main source files
export SRCPATH := $(PRJROOT)src
# Path to tests
export TESTPATH := $(PRJROOT)tests
# Path to tools
export TOOLPATH := $(PRJROOT)tools
# Path to parameter files
export PARPATH := $(PRJROOT)parameter_files

### Compiler + Linker Options #################################################
# Includes
INCLUDES += -I$(SRCPATH)
LIBS +=
LDFLAGS +=
# C++ flags
CPPFLAGS += 
CXXFLAGS += -g -Wall -O2 -std=c++11
CXXFLAGS_D = -g3 -Wall -fno-inline -O0
###############################################################################

### Shell stuff ###############################################################
# Quiet by default. Comment next line for verbose output
CMDPRE = @
MKDIR = $(CMDPRE)mkdir -pv
RM = $(CMDPRE)rm -rfv
PRINT = $(CMDPRE)printf
CXX := $(CMDPRE)$(CXX)
AR := $(CMDPRE)$(AR)
###############################################################################

# Find source files
SOURCES = $(shell find $(SRCPATH) -name '*.$(SRCEXT)' -printf '%T@\t%p\n' \
	| sort -k 1nr | cut -f2-)
OBJECTS = $(SOURCES:$(SRCPATH)/%.$(SRCEXT)=$(BUILDPATH)/%.o)
OBJECTS_D = $(OBJECTS:%.o=%_D.o)

