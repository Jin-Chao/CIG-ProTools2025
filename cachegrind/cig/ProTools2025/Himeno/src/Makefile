# Directories
EXP_DIR := ../..
PROJ_DIR := $(EXP_DIR)/Himeno
SRC_DIR := $(PROJ_DIR)/src
BIN_DIR := $(PROJ_DIR)/bin
OBJ_DIR := $(PROJ_DIR)/obj
UTIL_DIR := $(EXP_DIR)/Util
INC := -I$(UTIL_DIR)

# Compiler and flags
CC := gcc
RM := rm -f

SRC_FILES := himenoBMTxpa himenoBMTxpa_pad himenoBMTxpa_aos
OPT_LEVELS := O0 O3

CFLAGS.O0 := -fopenmp -g -O0
CFLAGS.O3 := -fopenmp -g -O3

UTIL_OBJ := $(OBJ_DIR)/util.o
UTIL_SRC := $(UTIL_DIR)/util.c

# Generate all binary targets
BINARIES := $(foreach opt,$(OPT_LEVELS), \
  $(foreach src,$(SRC_FILES), \
    $(BIN_DIR)/$(patsubst himenoBMTxpa%,bmt%,$(src)).$(opt)))

# Build all
all: $(BINARIES)

# Build shared util.o with -O3 (assumed stable across builds)
$(UTIL_OBJ): $(UTIL_SRC)
	$(CC) $(CFLAGS.O3) $(INC) -c -o $@ $<

# Rule for generating object files per optimization level
$(OBJ_DIR)/%.O0.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS.O0) $(INC) -c -o $@ $<

$(OBJ_DIR)/%.O3.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS.O3) $(INC) -c -o $@ $<

# Rule for building binaries from per-opt object files
define build_rule
$(BIN_DIR)/$(patsubst himenoBMTxpa%,bmt%,$(src)).$(opt): $(OBJ_DIR)/$(src).$(opt).o $(UTIL_OBJ)
	$$(CC) $$(CFLAGS.$(opt)) $$(INC) -o $$@ $$^
endef

# Evaluate build rules
$(foreach opt,$(OPT_LEVELS), \
  $(foreach src,$(SRC_FILES), \
    $(eval $(call build_rule,src=$(src),opt=$(opt)))))

# Clean
clean:
	$(RM) $(BIN_DIR)/* $(OBJ_DIR)/*.o

