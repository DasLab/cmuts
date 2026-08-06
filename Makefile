# Author: Hamish M. Blair <hmblair@stanford.edu>

CC       := cc
WARNINGS := -Wall -Wextra -Wpedantic
PACKAGES := htslib hdf5

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
