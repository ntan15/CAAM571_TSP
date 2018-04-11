# Copyright 2016, Gurobi Optimization, Inc.

GCC = '/bin/g++'
CPP = $(GCC)
PLATFORM = linux64
GUROBIDIR = /clear/apps/gurobi-7.0.1
INC      = $(GUROBIDIR)/include/
CLIB     = -Wl,-rpath,$(GUROBIDIR)/lib/ -lgurobi70
CPPLIB   = -L$(GUROBIDIR)/lib/ -lgurobi_c++
JSRC     = ../java
CLASSDIR = -classpath $(GUROBIDIR)/lib/gurobi.jar:.
JFLAG    = -d . $(CLASSDIR)
CFLAGS   = -Wall -std=c++14 -O2 -Wno-sign-compare


OBJECTS=$(SOURCE:.cc=.0)
SOURCE= tsp.cpp

tsp: $(OBJECTS)
	scl enable devtoolset-4 'g++ -m64 -g -o tsp $(OBJECTS) -I$(INC) $(CLIB) $(CPPLIB) -lpthread -lm $(CFLAGS)'


