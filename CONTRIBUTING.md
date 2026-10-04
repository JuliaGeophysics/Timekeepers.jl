# Contributing to Timekeepers.jl

Thank you for your interest in this project. This guide tells you how to
report a problem, suggest an improvement and send code.

## Reporting bugs

Open a [GitHub issue](../../issues) that gives:

- a short description of the problem;
- the steps that cause the problem again (input files, commands, Julia
  version);
- the full error message or the unexpected output.

For a problem with a reader or a writer, attach a short part of the file that
causes it. A few records and the header block are usually sufficient to cause
the parse failure again.

## Suggesting features

Open a GitHub issue with the label **enhancement**. Describe the use and the
behavior that you expect. We accept new instrument formats. The subsequent
sections tell you what a format module must supply.

## Submitting code

1. Fork the repository and make a branch from `main`.
2. Install the project: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
3. Make your changes.
4. Do the tests: `julia --project=. -e 'using Pkg; Pkg.test()'`
5. Open a pull request against `main`.

### Code style

- Use the standard Julia conventions (an indent of 4 spaces, function names
  in lowercase).
- Add a docstring for each new public function. The docs build uses
  `checkdocs = :exports`. Thus, an export without a docstring causes a CI
  failure.
- Each source file starts with a short comment that tells its purpose. When
  you change the purpose of a file, update this comment.
- Keep each commit to one logical change.

### Adding an instrument format

A format module must supply the same three entry points as the current
modules. Then `read_timekeeper` and `write_timekeeper` can use it without
special cases:

- `read_<format>(path; kwargs...) -> TimekeeperRun`
- `load_<format>(path; kwargs...) -> TimeArray`
- `write_<format>(path, run_or_timearray) -> String`

Then extend `_detect_format` and `_detect_output_format` in
`src/TimekeeperIO.jl`. The reader must keep the auxiliary columns in the
metadata. The writer must put them back in their original positions. Thus, a
cycle of read and write does not lose data. The current tests do a check of
this for each format.

### Tests

Add or update tests in `test/` for each new function. All the tests must pass
before we merge a pull request. The tests write synthetic files and read them
again. Thus, they need no external data. Use the same method. Do not use a
recording that is not in the repository.

### Documentation

The documentation is in `docs/src/`. Documenter and DocumenterVitepress make
the site. At the first use, DocumenterVitepress gets the Node tools that it
needs:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The build writes the site in `docs/build/`. Add each new exported function to
the correct `@docs` block in `docs/src/api.md`.

## Code of Conduct

Contributors must be respectful and constructive. We do not accept
harassment of any type.

## License

When you contribute, you agree that your contributions have the
[MIT License](LICENSE).
