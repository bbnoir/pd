default: all

CXX = g++
CXXFLAGS = -std=c++11 -fopenmp
EXE = bin/fp
SRCS_DIR = src
SRCS = $(wildcard $(SRCS_DIR)/*.cpp)
OBJS_DIR = build
OBJS = $(patsubst $(SRCS_DIR)/%.cpp, $(OBJS_DIR)/%.o, $(SRCS))

TEST_DIR = input_pa2
LOG_DIR = log
OUTPUT_DIR = out

DIRS = $(SRCS_DIR) $(OBJS_DIR) $(TEST_DIR) $(LOG_DIR) $(OUTPUT_DIR) bin

all: CXXFLAGS += -DDEBUG -O3 -DLOG 
all: $(DIRS) $(EXE)

debug: CXXFLAGS += -g -O3 -DDEBUG -DLOG
debug: $(DIRS) $(EXE)

release: CXXFLAGS += -O3 -DNDEBUG
release: $(DIRS) $(EXE)

$(EXE): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJS_DIR)/main.o: $(SRCS_DIR)/main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJS_DIR)/%.o: $(SRCS_DIR)/%.cpp $(SRCS_DIR)/%.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(DIRS):
	mkdir -p $@

clean:
	rm -f $(EXE) $(OBJS_DIR)/*.o output.rpt output.rpt.HPWL
.PHONY: clean

ALPHA = 0.5

define TEST
	./$(EXE) $(ALPHA) $(TEST_DIR)/$1.block $(TEST_DIR)/$1.nets $(OUTPUT_DIR)/$1.rpt
	bash evaluator/evaluator.sh $(TEST_DIR)/$1.block $(TEST_DIR)/$1.nets $(OUTPUT_DIR)/$1.rpt $(ALPHA)
endef

t1: default
	$(call TEST,ami33)
.PHONY: t1
