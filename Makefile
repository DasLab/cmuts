# Author: Hamish M. Blair <hmblair@stanford.edu>

# Assigned conditionally so that CC from the environment wins: on systems where
# the compiler comes from an environment module, plain cc is the untouched
# system one, which may be far older than the module's.
CC       ?= cc
WARNINGS := -Wall -Wextra -Wpedantic

# htslib and HDF5 are located with pkg-config where it can describe them, and
# otherwise assumed to be on the compiler's default search path, which is how
# environment modules present them. Neither route is universal: HPC HDF5
# builds often ship no .pc file at all, and htslib's declares private
# dependencies that older pkg-config insists on resolving even for a shared
# link. Any of these four may be overridden on the command line.
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
