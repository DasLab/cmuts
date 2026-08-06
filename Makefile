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

CFLAGS    := -std=c11 -O2 $(WARNINGS) -pthread -Iinclude $(HTS_CFLAGS) $(HDF5_CFLAGS) -MMD -MP
CMUTS_LIBS := $(HTS_LIBS) $(HDF5_LIBS) -pthread
GEN_LIBS   := $(HTS_LIBS)

# Defaults under the home directory so that installing needs no privileges and
# works on machines where /usr/local is not writable. DESTDIR is prepended only
# at install time, so that a package can be staged somewhere other than where
# it will finally live.
PREFIX  ?= $(HOME)/.local
BINDIR  ?= $(PREFIX)/bin
INSTALL ?= install

NAME     := cmuts
GEN_NAME := cmuts-gen

BUILD    := build
BIN      := $(BUILD)/$(NAME)
GEN_BIN  := $(BUILD)/$(GEN_NAME)

SRC      := $(wildcard src/*.c)
GEN_SRC  := $(wildcard tools/*.c)
OBJ      := $(SRC:src/%.c=$(BUILD)/src/%.o)
GEN_OBJ  := $(GEN_SRC:tools/%.c=$(BUILD)/tools/%.o)
DEP      := $(OBJ:.o=.d) $(GEN_OBJ:.o=.d)

# The generator shares the command line parser and nothing else: it has its own
# arguments, but no reason to read them differently.
GEN_SHARED := $(BUILD)/src/cli.o

all: $(BIN) $(GEN_BIN)

$(BIN): $(OBJ)
	$(CC) $(OBJ) -o $@ $(CMUTS_LIBS)

$(GEN_BIN): $(GEN_OBJ) $(GEN_SHARED)
	$(CC) $(GEN_OBJ) $(GEN_SHARED) -o $@ $(GEN_LIBS)

$(BUILD)/src/%.o: src/%.c | $(BUILD)/src
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD)/tools/%.o: tools/%.c | $(BUILD)/tools
	$(CC) $(CFLAGS) -Itools -c $< -o $@

$(BUILD)/src $(BUILD)/tools:
	mkdir -p $@

install: $(BIN) $(GEN_BIN)
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 755 $(BIN) $(DESTDIR)$(BINDIR)/$(NAME)
	$(INSTALL) -m 755 $(GEN_BIN) $(DESTDIR)$(BINDIR)/$(GEN_NAME)

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(NAME) $(DESTDIR)$(BINDIR)/$(GEN_NAME)

check: $(BIN) $(GEN_BIN)
	tests/run.sh

clean:
	rm -rf $(BUILD)

.PHONY: all check clean install uninstall

-include $(DEP)
