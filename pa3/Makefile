CC=g++
CXXFLAGS=-std=c++17 -static -O2 -Wall -D_GLIBCXX_ISE_CXX11_ABI=1  # for release
# CXXFLAGS=-std=c++17 -g -static -Wall -D_GLIBCXX_ISE_CXX11_ABI=1  # for debug
LDFLAGS=-Llib -lDetailPlace -lGlobalPlace -lLegalizer -lPlacement -lParser -lPlaceCommon
SOURCES=src/ObjectiveFunction.cpp src/Optimizer.cpp src/GlobalPlacer.cpp src/main.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=place

all: $(SOURCES) bin/$(EXECUTABLE)
	
bin/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(CXXFLAGS) $(LDFLAGS) -o $@

clean:
	rm -rf *.o bin/$(EXECUTABLE)

TEST_DIR=benchmark
LOG_DIR=log

define TEST
	./bin/$(EXECUTABLE) -aux $(TEST_DIR)/$1
endef

t1: all
	$(call TEST,ibm01/ibm01-cu85.aux)
.PHONY: t1

t2: all
	$(call TEST,ibm05/ibm05.aux)
.PHONY: t2