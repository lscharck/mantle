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

test: test.o smoke_test.o meta_test.o
	$(CXX) -L$(REBOUND) $(OPT) -Wl,-rpath,../rebound/src/ -lrebound test.o smoke_test.o meta_test.o -L$(LIBDIR) -lmantle -o test

test.o smoke_test.o meta_test.o: test.cpp smoke_test.cpp meta_test.cpp
	$(CXX) $(CXXFLAGS) -c meta_test.cpp -o meta_test.o
	$(CXX) $(CXXFLAGS) -c smoke_test.cpp -o smoke_test.o
	$(CXX) $(CXXFLAGS) -c test.cpp -o test.o

run_test:
	./test
	diff data.txt real.txt
	diff ode_data.txt ode_real.txt

clean:
	$(RM) -r test test.o smoke_test.o meta_test.o data.txt ode_data.txt

.PHONY: all run_test clean

