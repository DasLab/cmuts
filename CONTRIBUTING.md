# Contributing

## Setting up

Beyond the build dependencies in the README, the tests need `samtools` on `PATH` and a Python environment holding the test dependencies:

```sh
uv venv .venv
uv pip install --python .venv/bin/python --group dev
```

`make` uses `.venv/bin/python` where it exists and `python3` otherwise; set `PYTHON` to point somewhere else.

## Tests

```sh
make check
```

This builds the programs and runs pytest over `tests/` against the binaries in the build directory. The tests will reject a binary outside the repository, so directly invoking them requires prepending it to `PATH`:

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
