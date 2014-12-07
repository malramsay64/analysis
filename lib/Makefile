include ../system

BIN_PATH=../bin/$(SYS_NAME)

CXX = g++
CXXFLAGS := -O2 -std=c++11 #-std=c++0x
LDFLAGS = -lm -lpthread

ifeq ($(SYS_NAME), unix)
	CXXFLAGS := $(CXXFLAGS) -pthread -Wl,--no-as-needed
endif

simple=input.o particle.o vect.o frame.o functions.o main.o 
MODULES= input.o functions.o output.o neighbours.o frame.o parallel.o particle.o vect.o mol_functions.o main.o 

VPATH=.:$(BIN_PATH)

all: analysis 
	
analysis: $(MODULES)
	$(CXX) -o $(BIN_PATH)/analysis $(addprefix $(BIN_PATH)/, $(MODULES)) $(CXXFLAGS) $(LDFLAGS)
	ln -sf $(BIN_PATH)/analysis ../test/analysis

%.o : %.cpp | $(BIN_PATH)
	$(CXX) $(CXXFLAGS) -c $< -o $(BIN_PATH)/$@

$(BIN_PATH): 
	-mkdir -p $(BIN_PATH)

.PHONY: clean
clean:
	-rm $(BIN_PATH)/*
