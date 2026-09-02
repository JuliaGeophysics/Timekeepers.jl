# Contributing to Timekeepers.jl

Thank you for your interest in contributing! This guide covers how to
report issues, suggest improvements, and submit code.

## Reporting bugs

Open a [GitHub issue](../../issues) with:

- A short description of the problem.
- Steps to reproduce it (input files, commands, Julia version).
- The full error message or unexpected output.

For a reader or writer bug, the most useful thing you can attach is a short
excerpt of the offending file — a few records and the header block are usually
enough to reproduce a parsing failure.

## Suggesting features

Open a GitHub issue labelled **enhancement** describing the use case and
expected behaviour. New instrument formats are welcome; see below for what a
format module needs to provide.

## Submitting code

1. Fork the repository and create a branch from `main`.
2. Install the project: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
3. Make your changes.
4. Run the tests: `julia --project=. test/runtests.jl`
5. Open a pull request against `main`.

### Code style

- Follow standard Julia conventions (4-space indent, lowercase functions).
- Add docstrings for new public functions; the docs build runs with
  `checkdocs = :exports`, so an undocumented export fails CI.
- Each source file opens with a short comment saying what it is for — keep that
  up to date when you change a file's role.
- Keep commits focused; one logical change per commit.

### Adding an instrument format

A format module should provide the same three entry points as the existing
ones, so it slots into `read_timekeeper` / `write_timekeeper` without special
cases:

- `read_<format>(path; kwargs...) -> TimekeeperRun`
- `load_<format>(path; kwargs...) -> TimeArray`
- `write_<format>(path, run_or_timearray) -> String`

Then extend `_detect_format` and `_detect_output_format` in
`src/TimekeeperIO.jl`. Carry any auxiliary columns through the reader in
metadata and write them back in their original slots, so a read/write cycle is
lossless — the existing tests check exactly this for every format.

### Tests

Add or update tests in `test/` for any new functionality. All tests must
pass before a PR will be merged. The suite writes synthetic files and reads
them back, so no external data is needed — follow that pattern rather than
depending on a recording that is not in the repository.

### Documentation

Documentation lives in `docs/src/` and is built with Documenter:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The result is in `docs/build/`. New exported functions should be added to the
appropriate `@docs` block in `docs/src/api.md`.

## Code of Conduct

Contributors are expected to be respectful and constructive. Harassment
of any kind will not be tolerated.

## License

By contributing you agree that your contributions will be licensed under
the [MIT License](LICENSE).
