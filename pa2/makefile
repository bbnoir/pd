default: release

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
	./$(EXE) $(ALPHA) $(TEST_DIR)/$1.block $(TEST_DIR)/$1.nets $(OUTPUT_DIR)/$1.rpt | tee $(LOG_DIR)/$1.log
	bash evaluator/evaluator.sh $(TEST_DIR)/$1.block $(TEST_DIR)/$1.nets $(OUTPUT_DIR)/$1.rpt $(ALPHA) | tee -a $(LOG_DIR)/$1.log
endef

t1: default
	$(call TEST,ami33)
.PHONY: t1

t2: default
	$(call TEST,ami49)
.PHONY: t2

t3: default
	$(call TEST,apte)
.PHONY: t3

t4: default
	$(call TEST,hp)
.PHONY: t4

t5: default
	$(call TEST,xerox)
.PHONY: t5

tall: default t1 t2 t3 t4 t5
.PHONY: tall

define VISUALIZE
	python3 scripts/draw.py -b $(TEST_DIR)/$1.block -r $(OUTPUT_DIR)/$1.rpt -o viz/$1.png -t $1
endef

v1:
	$(call VISUALIZE,ami33)
.PHONY: v1

v2:
	$(call VISUALIZE,ami49)
.PHONY: v2

v3:
	$(call VISUALIZE,apte)
.PHONY: v3

v4:
	$(call VISUALIZE,hp)
.PHONY: v4

v5:
	$(call VISUALIZE,xerox)
.PHONY: v5

vall: v1 v2 v3 v4 v5
.PHONY: vall