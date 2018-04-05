#!/bin/sh
# COMPILER = g++

# COMPILATION_FLAGS = -c -frounding-math -O3 -I. -I./Headers -std=c++11
# LINKING_FLAGS = -static -lCGAL -lboost_filesystem -lboost_system -lboost_program_options -lboost_thread -lpthread

# trackLoop: trackLoop.o SimplicialComplex.o
# 	$(COMPILER) trackLoop.o -o trackLoop $(LINKING_FLAGS)

# trackLoop.o: trackLoop.cpp SimplicialComplex.cpp
# 	$(COMPILER) trackLoop.cpp -o trackLoop.o $(COMPILATION_FLAGS)

# clean:
# 	\rm -f *.o trackLoop

LIBS="-I/usr/local/include -I/usr/local/include/eigen3 -static -Igmpq -lCGAL -Igmp -lgmp -lann -lboost_system -lboost_filesystem -lboost_program_options -lboost_thread -lpthread"

g++ -O3 -frounding-math -I. -I./Headers -o trackLoop --std=c++11  *.cpp  $LIBS
#  -O3

# g++  -O3 --std=c++11  *.cpp Wrappers/*.cpp SimPers/*.cpp Graphs/*.cpp GIComplex/*.cpp  $LIBS -o sibaco-dec3