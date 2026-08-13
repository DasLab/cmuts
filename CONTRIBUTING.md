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

`SAN` names a sanitizer and integrates with the tests:

```sh
make SAN=asan         # address and undefined behavior
make check SAN=tsan   # thread
```

## Documentation

Documentation on program arguments and the output HDF5 structure is automatically generated from the binaries.

```sh
make docs   # rewrite the generated blocks under docs/
make site   # additionally render the site into site/
```

To preview the site locally, rebuilding as the sources change:

```sh
.venv/bin/sphinx-autobuild docs site
```

# Contributing

## Writing Tests

Tests which run a program over an alignment are parametrized over the alignments specified in `tests/datasets.py`. When adding a new test, ensure it follows this parametrization and that it tests a contract which holds regardless of the input alignment.

Tests which build their own inputs are the exception, and the generator produces neither kind: `test_subtraction.py` writes the output files it reads, and `test_marginalization.py` writes alignments with specific CIGARs.

To ensure tests are not vacuous, the `checked` fixture is used to signal that the test is falsifiable on the current dataset. A test in which no dataset passes a truthy object to `checked` fails the run, signalling the need for a new dataset in `tests/datasets.py`.
