EXP_DIR := ..
PROJ_DIR := ${EXP_DIR}/PolyBenchC-4.2.1

SRC_DIR := $(PROJ_DIR)/datamining/correlation
UTIL_DIR := ${EXP_DIR}/Util
UTILS_DIR := $(PROJ_DIR)/utilities

BIN_DIR := $(PROJ_DIR)/bin
OBJ_DIR := $(PROJ_DIR)/obj

CC := gcc
CFLAGS := -O3 -g -Wall
LIB := -lm
INC := -I$(UTIL_DIR) -I$(UTILS_DIR)

SIZES := small medium large extra
SIZE_FLAGS.small := -DSMALL_DATASET
SIZE_FLAGS.medium := -DMEDIUM_DATASET
SIZE_FLAGS.large := -DLAREGE_DATASET
SIZE_FLAGS.extra := -DEXTRALARGE_DATASET

TIME_INS := -DPOLYBENCH_TIME

UTIL_SRCS := $(UTIL_DIR)/util.c $(UTILS_DIR)/polybench.c

TARGETS := $(foreach size,$(SIZES),correlation.$(size) correlation_block.$(size))

all: $(addprefix $(BIN_DIR)/, $(TARGETS))

$(BIN_DIR)/correlation.%.out: $(SRC_DIR)/correlation.c $(UTIL_SRCS)
	$(CC) $(CFLAGS) $(INC) $(TIME_INS) $(SIZE_FLAGS.$*) -o $@ $^ $(LIB)

$(BIN_DIR)/correlation_block.%.out: $(SRC_DIR)/correlation_block.c $(UTIL_SRCS)
	$(CC) $(CFLAGS) $(INC) $(TIME_INS) $(SIZE_FLAGS.$*) -o $@ $^ $(LIB)

# Provide compatibility if you still want "correlation.small" style naming
$(BIN_DIR)/correlation.%: $(BIN_DIR)/correlation.%.out
	cp $< $@

$(BIN_DIR)/correlation_block.%: $(BIN_DIR)/correlation_block.%.out
	cp $< $@

.PHONY: clean
clean:
	rm -f $(BIN_DIR)/* $(OBJ_DIR)/*


