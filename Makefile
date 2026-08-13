# Author: Hamish M. Blair <hmblair@stanford.edu>

# Conditional so CC from the environment wins over the system cc.
CC       ?= cc

OPT      ?= -O2

# One directory per variant, none of them inside another, so that an object
# from an ordinary build cannot link into a binary that checks nothing and each
# survives the others being built.
BUILD_ROOT := build
VARIANT    := release

# A sanitizer is chosen by name, which stands for a variant of its own and the
# flags that fill it: `make SAN=asan` builds one and `make check SAN=asan` runs
# the tests over it. Thread and address cannot be combined, and neither finds
# leaks on macOS, which has no LeakSanitizer.
SANITIZE_asan := address,undefined
SANITIZE_tsan := thread

ifdef SAN
ifeq ($(SANITIZE_$(SAN)),)
$(error unknown sanitizer: $(SAN))
endif
VARIANT  := $(SAN)
OPT      := -O1
endif

BUILD    := $(BUILD_ROOT)/$(VARIANT)

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

# A sanitizer has to reach both the compiler and the linker; what it adds to
# the link is attached to each program's libraries, below. Halting on the first
# diagnostic is what makes one visible under the tests, which keep the output of
# a program that exits cleanly; undefined behavior is otherwise reported and
# stepped over.
ifdef SAN
SANFLAGS   := -fsanitize=$(SANITIZE_$(SAN)) -fno-sanitize-recover=all \
              -fno-omit-frame-pointer -g
CFLAGS     += $(SANFLAGS)
endif

# Defaults under the home directory so that installing needs no privileges and
# works on machines where /usr/local is not writable. DESTDIR is prepended only
# at install time, so that a package can be staged somewhere other than where
# it will finally live.
PREFIX  ?= $(HOME)/.local
BINDIR  ?= $(PREFIX)/bin
INSTALL ?= install

NAME     := cmuts

# One directory under apps/ per program, named for the binary it builds. A
# program is its own sources, its own private headers, and whichever members of
# the library it refers to; adding one means adding a directory and a word here.
PROGRAMS := cmuts-hmm cmuts-gen cmuts-sub

LIB      := $(BUILD)/lib$(NAME).a

# Sources under src/ are the library and nothing else: no entry point lives
# there, so every program can link the whole of it.
LIB_SRC  := $(wildcard src/*.c)
LIB_OBJ  := $(LIB_SRC:src/%.c=$(BUILD)/src/%.o)

app_sources = $(wildcard apps/$(1)/*.c)
app_objects = $(patsubst apps/$(1)/%.c,$(BUILD)/apps/$(1)/%.o,$(call app_sources,$(1)))

APP_SRC  := $(foreach p,$(PROGRAMS),$(call app_sources,$(p)))
APP_OBJ  := $(foreach p,$(PROGRAMS),$(call app_objects,$(p)))
BINS     := $(addprefix $(BUILD)/,$(PROGRAMS))
SRC      := $(LIB_SRC) $(APP_SRC)
DEP      := $(LIB_OBJ:.o=.d) $(APP_OBJ:.o=.d)

# What each program needs beyond the library. The archive is searched, not
# swallowed, so a program links only the members it refers to and only the
# libraries those in turn require: the generator reaches cli.o alone and so
# needs nothing of HDF5, and the subtraction reads and writes output files
# without ever opening an alignment.
LIBS_cmuts-hmm := $(HTS_LIBS) $(HDF5_LIBS) -pthread -lm
LIBS_cmuts-gen := $(HTS_LIBS)
LIBS_cmuts-sub := $(HDF5_LIBS) -lm

ifdef SAN
$(foreach p,$(PROGRAMS),$(eval LIBS_$(p) += $(SANFLAGS)))
endif

all: $(BINS)

# Objects first and the archive last, which is the order a linker resolves in.
define program_rule
$(BUILD)/$(1): $(call app_objects,$(1)) $$(LIB)
	$$(CC) $$^ -o $$@ $$(LIBS_$(1))
endef

$(foreach p,$(PROGRAMS),$(eval $(call program_rule,$(p))))

$(LIB): $(LIB_OBJ)
	$(AR) rcs $@ $^

# A program's own directory is on its include path, so a header beside its
# source is private to it and one in include/ is shared.
$(BUILD)/apps/%.o: apps/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I$(<D) -c $< -o $@

$(BUILD)/src/%.o: src/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

install: $(BINS)
	$(INSTALL) -d $(DESTDIR)$(BINDIR)
	$(INSTALL) -m 755 $^ $(DESTDIR)$(BINDIR)

uninstall:
	rm -f $(addprefix $(DESTDIR)$(BINDIR)/,$(PROGRAMS))

# Prefers a virtual environment holding the test dependencies, which the
# programs themselves do not need:
#
#     uv venv .venv && uv pip install --python .venv/bin/python --group dev
PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)

# The build directory goes first on PATH, which is how the tests reach these
# programs and not an installed copy of them. They refuse anything found
# outside this repository, so the tests cannot be run without it.
check: $(BINS)
	PATH=$(CURDIR)/$(BUILD):$$PATH $(PYTHON) -m pytest

# Every page under docs/ naming a program in one of its markers is written
# from that program, so this builds them all first: a table can only describe
# what was built here, never a version of it that is no longer around.
docs: $(BINS)
	@$(PYTHON) scripts/document.py $(BUILD) docs

# The blocks are rewritten first, so the site cannot render a page describing a
# program as it was. --strict fails on a broken link or a page left out.
#
#     uv pip install --python .venv/bin/python --group docs
site: docs
	@$(PYTHON) -m mkdocs build --strict

# A clang that is not the system one needs the SDK spelled out.
CLANG_TIDY ?= $(firstword $(wildcard /opt/homebrew/opt/llvm/bin/clang-tidy) clang-tidy)
SDK        := $(if $(filter Darwin,$(shell uname)),-isysroot $(shell xcrun --show-sdk-path),)
DB_FLAGS   := $(CFLAGS) $(SDK)

compile_commands.json: Makefile $(SRC)
	@{ printf '[\n'; \
	   first=1; \
	   for f in $(SRC); do \
	     [ $$first -eq 1 ] || printf ',\n'; \
	     first=0; \
	     printf '  { "directory": "%s", "file": "%s", "command": "%s" }' \
	       "$(CURDIR)" "$$f" "$(CC) $(DB_FLAGS) -I$$(dirname $$f) -c $$f"; \
	   done; \
	   printf '\n]\n'; } > $@

lint: compile_commands.json
	@$(CLANG_TIDY) -p . --quiet $(SRC)

# Clear every variant, not just the one this invocation names.
clean:
	rm -rf $(BUILD_ROOT) compile_commands.json

.PHONY: all check clean docs install lint site uninstall

-include $(DEP)
