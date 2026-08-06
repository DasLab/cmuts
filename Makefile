# Author: Hamish M. Blair <hmblair@stanford.edu>

# Conditional so CC from the environment wins over the system cc.
CC       ?= cc
WARNINGS := -Wall -Wextra -Wpedantic

# pkg-config where it works, plain -l otherwise; environment modules often
# ship no .pc file. Overridable on the command line.
HTS_CFLAGS  := $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LIBS    := $(shell pkg-config --libs   htslib 2>/dev/null || echo -lhts)
HDF5_CFLAGS := $(shell pkg-config --cflags hdf5   2>/dev/null)
HDF5_LIBS   := $(shell pkg-config --libs   hdf5   2>/dev/null || echo -lhdf5)

CFLAGS   := -std=c11 -O2 $(WARNINGS) -pthread -Iinclude $(HTS_CFLAGS) $(HDF5_CFLAGS) -MMD -MP
LDLIBS   := $(HTS_LIBS) $(HDF5_LIBS) -pthread

BUILD := build
BIN   := $(BUILD)/cmuts
SRC   := $(wildcard src/*.c)
OBJ   := $(SRC:src/%.c=$(BUILD)/%.o)
DEP   := $(OBJ:.o=.d)

all: $(BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDLIBS)

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)

.PHONY: all clean

-include $(DEP)
