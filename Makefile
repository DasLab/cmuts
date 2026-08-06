# Author: Hamish M. Blair <hmblair@stanford.edu>

# Conditional so CC from the environment wins over the system cc.
CC       ?= cc

OPT      ?= -O2

# Clean as it stands. Left out: -Wconversion and -Wsign-conversion (three
# sites, two of them htslib's KSEQ_INIT expansion in fasta.c) and
# -Wswitch-enum (one switch in cli.c relying on its default).
WARNINGS := -Wall -Wextra -Wpedantic -Wshadow -Wmissing-prototypes \
            -Wstrict-prototypes -Wpointer-arith -Wwrite-strings -Wundef \
            -Wvla -Wdouble-promotion -Wcast-qual -Wformat=2 \
            -Wredundant-decls -Wnull-dereference

# pkg-config where it works, plain -l otherwise; environment modules often
# ship no .pc file. Overridable on the command line.
# -isystem rather than -I, so the warnings and clang-tidy apply here and not to
# htslib and HDF5.
HTS_CFLAGS  := $(subst -I,-isystem ,$(shell pkg-config --cflags htslib 2>/dev/null))
HTS_LIBS    := $(shell pkg-config --libs   htslib 2>/dev/null || echo -lhts)
HDF5_CFLAGS := $(subst -I,-isystem ,$(shell pkg-config --cflags hdf5   2>/dev/null))
HDF5_LIBS   := $(shell pkg-config --libs   hdf5   2>/dev/null || echo -lhdf5)

# Where there is no .pc file the headers arrive through CPATH instead, which
# the compiler treats as -I and so warns about. Restated as -isystem, which
# wins for a directory named both ways.
CPATH_CFLAGS := $(addprefix -isystem ,$(subst :, ,$(CPATH)))

CFLAGS    := -std=c11 $(OPT) $(WARNINGS) -pthread -Iinclude $(HTS_CFLAGS) $(HDF5_CFLAGS) $(CPATH_CFLAGS) -MMD -MP
CMUTS_LIBS := $(HTS_LIBS) $(HDF5_LIBS) -pthread -lm
GEN_LIBS   := $(HTS_LIBS)

# A sanitizer has to reach both the compiler and the linker.
ifdef SANITIZE
SANFLAGS   := -fsanitize=$(SANITIZE) -fno-omit-frame-pointer -g
CFLAGS     += $(SANFLAGS)
CMUTS_LIBS += $(SANFLAGS)
GEN_LIBS   += $(SANFLAGS)
endif

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

# cmuts alone. The generator writes test fixtures and benchmark inputs, which
# is work done from the build tree; nothing looks for it on PATH.
install: $(BIN)
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 755 $(BIN) $(DESTDIR)$(BINDIR)/$(NAME)

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(NAME)

# Prefers a virtual environment holding the test dependencies, which the
# programs themselves do not need:
#
#     uv venv .venv && uv pip install --python .venv/bin/python --group dev
PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)

check: $(BIN) $(GEN_BIN)
	$(PYTHON) -m pytest

# Separate object trees, so that an object from an ordinary build cannot link
# into a binary that checks nothing. Thread and address cannot be combined, and
# neither finds leaks on macOS, which has no LeakSanitizer.
tsan:
	@$(MAKE) BUILD=$(BUILD)/tsan OPT=-O1 SANITIZE=thread all

asan:
	@$(MAKE) BUILD=$(BUILD)/asan OPT=-O1 SANITIZE=address,undefined all

# clang-tidy reads the compile database, not this file. A clang that is not the
# system one also needs the SDK spelled out.
CLANG_TIDY ?= $(firstword $(wildcard /opt/homebrew/opt/llvm/bin/clang-tidy) clang-tidy)
SDK        := $(if $(filter Darwin,$(shell uname)),-isysroot $(shell xcrun --show-sdk-path),)
TIDY_FLAGS := -std=c11 $(SDK) -Iinclude $(HTS_CFLAGS) $(HDF5_CFLAGS)

compile_commands.json: Makefile $(SRC)
	@{ printf '[\n'; \
	   first=1; \
	   for f in $(SRC); do \
	     [ $$first -eq 1 ] || printf ',\n'; \
	     first=0; \
	     printf '  { "directory": "%s", "file": "%s", "command": "%s" }' \
	       "$(CURDIR)" "$$f" "$(CC) $(TIDY_FLAGS) -c $$f"; \
	   done; \
	   printf '\n]\n'; } > $@

lint: compile_commands.json
	@$(CLANG_TIDY) -p . --quiet $(SRC)

clean:
	rm -rf $(BUILD) compile_commands.json

.PHONY: all asan check clean install lint tsan uninstall

-include $(DEP)
