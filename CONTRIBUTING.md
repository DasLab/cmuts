# Environment and Builds

## Setting up

Beyond the build dependencies in the README, the tests need `samtools` on `PATH` and a Python environment holding the test dependencies:

```sh
uv venv .venv
uv pip install --python .venv/bin/python --group dev
```

Additionally, rendering the documentation requires the `docs` group:

```sh
uv pip install --python .venv/bin/python --group docs
```

## Tests

```sh
make check
```

This builds the programs and runs pytest over `tests/` against the binaries in the build directory. Pass `PYTHON` to `make` to use a binary other than `.venv/bin/python`.

The tests will reject a binary outside the repository, so directly invoking them requires prepending it to `PATH`:

```sh
PATH=$PWD/build/release:$PATH .venv/bin/python -m pytest tests/test_filtering.py -x
```

## Linter

```sh
make lint
```

This writes `compile_commands.json` and runs clang-tidy over every source under `src/` and `apps/` with the checks in `.clang-tidy`. Set `CLANG_TIDY` to change the linter binary.

## Sanitizers

`SAN` names a sanitizer to build under, and applies to whichever target follows it:

```sh
make SAN=asan         # address and undefined behavior
make check SAN=tsan   # thread, over the tests
```

The tests capture the standard error of every program they run, and a failure reports only the exit status, so `make check` sends each report to a file of its own under the build directory of the variant:

```
build/tsan/logs/tsan.<pid>
```

A value already in `ASAN_OPTIONS` or `TSAN_OPTIONS` is read after the one `make` sets and so overrides it.

## Documentation

Documentation on program arguments and the output HDF5 structure is automatically generated from the binaries.

```sh
make docs   # rewrite the generated blocks under docs/
make site   # additionally render the site into site/
```

If a program's arguments or output format changes, you must run `make docs` and include its changes in the commit.

To preview the site locally:

```sh
.venv/bin/sphinx-autobuild docs site
```

## Static Binaries

`scripts/static-deps.sh` builds htslib, HDF5 and their compression libraries as static archives under one prefix, with sources pinned by version and checksum:

```sh
scripts/static-deps.sh deps-static
make STATIC_PREFIX=deps-static
```

On macOS system libraries stay dynamic; a musl compiler (`scripts/static-deps.sh deps-static musl-gcc`, or any Alpine gcc) makes the whole link static. `make check STATIC_PREFIX=deps-static` runs the tests against the result.

Pushing a `v*` tag runs the `Release` workflow, which builds these binaries for Linux (x86_64, aarch64, via Alpine containers) and macOS (arm64), runs the test suite against each, and attaches `scripts/package.sh` tarballs to a draft GitHub release.

# Contributing

## Writing Tests

Tests which run a program over an alignment take the `data` fixture, which is parametrized over the alignments specified in `tests/datasets.py` and over every format htslib reads, yielding one `Dataset` per combination. When adding a new test, ensure it takes this fixture and that it tests a contract which holds regardless of the input alignment and of the format that alignment arrives in. The `catalogue` fixture yields a lookup by name and format, for the few tests needing a second dataset.

A test requiring specific formats pins `fmt` instead of taking all three. A whole module pins it by overriding the fixture, whereas a single test pins it with `indirect`.

To ensure tests are not vacuous, the `falsifiable` fixture is used to signal whether the test is falsifiable on the current dataset. Place the call between the setup and the assertions of the test. A test in which no dataset declares itself falsifiable fails the run, signalling the need for a new dataset in `tests/datasets.py`.

## Continuous Integration

Two workflows run on every push to `main` and on every pull request against it.

`CI` runs `make check` under gcc and clang on Linux and under clang on macOS, and `make lint` on Linux. The sanitizers are not among them, taking longer than the rest of the run put together; run those locally.

`Documentation` renders the site, and fails if the generated docs are not current, rather than automatically updating them for you. A push to `main` deploys the rendered site as well.
