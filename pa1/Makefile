CC=g++
LDFLAGS=-std=c++11 -O3 -lm -g -fopenmp
SOURCES=src/partitioner.cpp src/main.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=fm
INCLUDES=src/cell.h src/net.h src/partitioner.h

all: $(SOURCES) bin/$(EXECUTABLE)

bin/$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o:  %.c  ${INCLUDES}
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o bin/$(EXECUTABLE)

re: clean all
.PHONY: re

TEST_DIR = input_pa1
LOG_DIR = log
OUTPUT_DIR = out

CUSTOM_DIRS = $(LOG_DIR) $(OUTPUT_DIR)

$(CUSTOM_DIRS):
	mkdir -p $@

dirs: $(CUSTOM_DIRS)

clean_custom:
	rm -rf $(CUSTOM_DIRS)

define TEST
	$(eval RUNTIME := $(shell { time -p ./bin/$(EXECUTABLE) $(TEST_DIR)/$1.dat $(OUTPUT_DIR)/$1.dat > $(LOG_DIR)/$1.log; } 2>&1 | grep real | awk '{print $$2}'))
	@cat $(LOG_DIR)/$1.log
	@echo "Runtime: $(RUNTIME) seconds"
	./evaluator/evaluator.sh $(TEST_DIR)/$1.dat $(OUTPUT_DIR)/$1.dat $(RUNTIME) | tee $(LOG_DIR)/$1_eval.log
endef

t1: all dirs
	$(call TEST,input_0)
.PHONY: t1

t2: all dirs
	$(call TEST,input_1)
.PHONY: t2

t3: all dirs
	$(call TEST,input_2)
.PHONY: t3

t4: all dirs
	$(call TEST,input_3)
.PHONY: t4

t5: all dirs
	$(call TEST,input_4)
.PHONY: t5

t6: all dirs
	$(call TEST,input_5)
.PHONY: t6

tall: clean_custom re t1 t2 t3 t4 t5 t6
.PHONY: tall