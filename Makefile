# Author: Hamish M. Blair <hmblair@stanford.edu>

CC       := cc
WARNINGS := -Wall -Wextra -Wpedantic
# Homebrew ships HDF5 built against MPI, whose headers pull in mpi.h. Where
# that is the build in use, open-mpi's flags are needed too, even though
# nothing here is parallel at the MPI level. A serial HDF5 needs no such thing,
# so the package is added only when it is present.
PACKAGES := htslib hdf5
PACKAGES += $(shell pkg-config --exists ompi && echo ompi)

CFLAGS   := -std=c11 -O2 $(WARNINGS) -pthread -Iinclude $(shell pkg-config --cflags $(PACKAGES)) -MMD -MP
LDLIBS   := $(shell pkg-config --libs $(PACKAGES)) -pthread

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
