CXX=c++
CXXFLAGS=-g3 -Wall -Wno-unknown-pragmas -std=c++23 -I$(INC) -I$(REBOUNDINC)
SRC=src
OBJ=obj
INC=include
TESTDIR=test

LIBDIR=lib
LIB=$(LIBDIR)/libmantle.a

REBOUNDDIR=rebound
REBOUNDINC=$(REBOUNDDIR)/src

SET_OPENMP = 0

ifeq ($(SET_OPENMP), 1)
	OPT=-DOPENMP
	OPT+=-fopenmp
	export OPENMP=1
endif

all: rebound $(LIB)

release: CXXFLAGS=-O3 -DNDEBUG -Werror -std=c++23 -pedantic -Wall -Wno-unknown-pragmas -I$(INC) -I$(REBOUNDINC)
release: clean
release: rebound $(LIB)

$(LIBDIR)/libmantle.a: $(OBJ)/mantle.o $(OBJ)/mantle_config.o $(OBJ)/mantle_particle.o
	ar rvs $(LIBDIR)/libmantle.a $(OBJ)/mantle.o $(OBJ)/mantle_config.o $(OBJ)/mantle_particle.o

$(OBJ)/mantle.o $(OBJ)/mantle_config.o $(OBJ)/mantle_particle.o: $(SRC)/mantle.cpp $(SRC)/mantle_config.cpp $(SRC)/mantle_particle.cpp
	$(CXX) $(CXXFLAGS) -c $(SRC)/mantle_particle.cpp -o $(OBJ)/mantle_particle.o
	$(CXX) $(CXXFLAGS) -c $(OPT) $(SRC)/mantle_config.cpp -o $(OBJ)/mantle_config.o
	$(CXX) $(CXXFLAGS) -c $(OPT) $(SRC)/mantle.cpp -o $(OBJ)/mantle.o

rebound:
	echo "Building rebound"
	$(MAKE) -C $(REBOUNDDIR)

test:
	echo "Building test"
	$(MAKE) -C $(TESTDIR) SET_OPENMP=$(SET_OPENMP) CXX=$(CXX)

clean:
	$(RM) -r $(LIBDIR)/* $(OBJ)/*

.PHONY: all rebound test clean

