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

Tests which run a program over an alignment take the `data` fixture, which is parametrized over the alignments specified in `tests/datasets.py` and yields one `Dataset` per test. When adding a new test, ensure it takes this fixture and that it tests a contract which holds regardless of the input alignment. The `catalogue` fixture yields a lookup by name, for the few tests needing a second dataset.

Tests which build their own inputs are the exception, and the generator produces neither kind: `test_subtraction.py` writes the output files it reads, and `test_marginalization.py` writes alignments with specific CIGARs.

To ensure tests are not vacuous, the `falsifiable` fixture is used to signal whether the test is falsifiable on the current dataset. Place the call between the setup and the assertions of the test. A test in which no dataset declares itself falsifiable fails the run, signalling the need for a new dataset in `tests/datasets.py`.

## Test Support Modules

Each module under `tests/` owns one concern, and a test file should reach for the one that matches what it needs:

| Module | Holds |
| --- | --- |
| `programs.py` | Locating and running the binaries and samtools. Every subprocess starts here. |
| `alignments.py` | The `Dataset` type, generating one with `cmuts-gen`, and the transforms over one. |
| `oracle.py` | The expected result, derived from the alignment and the FASTA by way of samtools. |
| `outputs.py` | The HDF5 layout: the field declarations, writing a file, and reading one back. |
| `subtraction.py` | What `cmuts-sub` should compute from its inputs. |
| `datasets.py` | The catalogue of alignments every test runs over. |
| `filters.py` | The read-filtering criteria a run is given. |

`outputs.py` declares each field once in `FIELDS`, and derives the groups tests reason about (`PADDED_FIELDS`, `COUNT_FIELDS`, `PER_REFERENCE_FIELDS`, and so on) from that declaration. Adding a field to the output means adding it there, not to a list in a test.
