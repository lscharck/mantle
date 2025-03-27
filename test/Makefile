CXX=c++
CXXFLAGS=-g3 -Wall -std=c++23 -I$(INC)
INC=../include

LIBDIR=../lib
LIB=$(LIBDIR)/libmantle.a

REBOUND=../rebound/src

SET_OPENMP = 0
ifeq ($(SET_OPENMP), 1)
	OPT=-fopenmp
endif

all: clean test run_test

test: test.o
	$(CXX) -L$(REBOUND) $(OPT) -Wl,-rpath,../rebound/src/ -lrebound test.o -L$(LIBDIR) -lmantle -o test

test.o: test.cpp
	$(CXX) $(CXXFLAGS) -c test.cpp -o test.o

run_test:
	./test
	diff data.txt real.txt

clean:
	$(RM) -r test test.o data.txt

.PHONY: all run_test clean

