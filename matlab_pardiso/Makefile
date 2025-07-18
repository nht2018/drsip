# I use gcc/9.1.0 to compile on linux server


OS_NAME := $(shell uname -s)
ifeq ($(OS_NAME),Linux)
$(info *** Linux operating system detected. ***)
# OS specific variables
MEXSUFFIX   = mexa64
MEX         = mex
CXX         = g++
else 
$(warning The mex version of Panua-Pardiso is currently available for Linux.)	
endif

# Recommended location for the shared library is the home directory
PARDISOHOME = ${HOME}/pardiso/panua-pardiso-20230908-linux/lib
$(info *** Looking for the shared Pardiso library at "$(PARDISOHOME)". ***)


# Do not edit anything past here (unless you know what you're doing).
MEXFLAGS   = -cxx CXX=$(CXX) CC=$(CXX) LD=$(CXX)            \
             -L$(PARDISOHOME) -lpardiso -lmwlapack -lmwblas \
             -lgfortran -lpthread -lm -R2017b

TARGETS = pardisoinit.$(MEXSUFFIX) pardisoreorder.$(MEXSUFFIX) \
          pardisofactor.$(MEXSUFFIX) pardisosolve.$(MEXSUFFIX) \
          pardisofree.$(MEXSUFFIX)

all: $(TARGETS)

%.$(MEXSUFFIX): %.cpp common.cpp matlabmatrix.cpp sparsematrix.cpp \
                pardisoinfo.cpp
	$(MEX) $(MEXFLAGS) -output $@ $^

clean:
	rm -f $(TARGETS)

