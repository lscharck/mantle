CXX=c++
CXXFLAGS=-g3 -Wall -std=c++23 -I$(INC) -I$(REBOUNDINC)
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
	CXXFLAGS+=-fopenmp
endif

all: rebound $(LIB)

release: CXXFLAGS=-O2 -DNDEBUG -Werror -std=c++23 -pedantic -Wall -I$(INC) -I$(REBOUNDINC)
release: clean
release: rebound $(LIB)

$(LIBDIR)/libmantle.a: $(OBJ)/mantle.o $(OBJ)/mantle_particle.o $(OBJ)/mantle_particle_restricted.o
	ar rvs $(LIBDIR)/libmantle.a $(OBJ)/mantle.o $(OBJ)/mantle_particle.o $(OBJ)/mantle_particle_restricted.o

$(OBJ)/mantle.o $(OBJ)/mantle_particle.o $(OBJ)/mantle_particle_restricted.o: $(SRC)/mantle.cpp $(SRC)/mantle_particle.cpp $(SRC)/mantle_particle_restricted.cpp
	$(CXX) $(CXXFLAGS) -c $(SRC)/mantle_particle.cpp -o $(OBJ)/mantle_particle.o
	$(CXX) $(CXXFLAGS) -c $(SRC)/mantle_particle_restricted.cpp -o $(OBJ)/mantle_particle_restricted.o
	$(CXX) $(CXXFLAGS) -c $(OPT) $(SRC)/mantle.cpp -o $(OBJ)/mantle.o

rebound:
	echo "Building rebound"
	$(MAKE) -C $(REBOUNDDIR)

test:
	echo "Building test"
	$(MAKE) -C $(TESTDIR)

clean:
	$(RM) -r $(LIBDIR)/* $(OBJ)/*

.PHONY: all rebound test clean

