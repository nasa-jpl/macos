## CLI syntax conventions

### Command names and minimum match

Command names are shown with a **capitalized prefix and lowercase
tail**: `PERturb`, `SUMmarize`, `PROpagate`.  The capitalized part is
the *minimum-match abbreviation* — the dispatcher accepts any prefix
at least that long, case-insensitively.  So `per`, `pert`, and
`perturb` all invoke PERturb.  When two commands share a prefix the
longer spelling wins disambiguation; when in doubt type the full name.

### The dialog model

Most commands take no arguments on the command line.  Instead the
engine prompts for each input in turn, showing the expected form and
often a default in brackets:

```
MACOS>perturb
Enter element number [1]: 3
Enter 6-vector perturbation:
 (dx dy dz rx ry rz): 0 0 1e-6 0 0 0
```

- Hitting return accepts the bracketed default.
- Vectors are entered space-separated on one line.
- Fortran-style numbers are accepted: `1e-6`, `1d-6`, `0.25`, `1d22`.
- `Ctrl-C` at the `MACOS>` prompt exits cleanly.

### Sessions, journals, prescriptions

A typical first session:

```
$ macos
Enter model size: 256
MACOS>old
Enter file name: Cassegrain
MACOS>summarize
MACOS>spot
 ...
MACOS>quit
```

- The **model size** prompt (before the first `MACOS>`) sizes the
  runtime arrays; 256 is a good default for the bundled examples.
- `OLD` loads an existing `.in` prescription (the `.in` extension is
  implied); `NEW` builds one interactively; `SAVe` writes the current
  system back out.
- `JOUrnal` starts recording your commands to a `.jou` file;
  `EXEcute` replays one.  The examples directory pairs each
  `<name>.in` with a `<name>.jou` that drives a complete analysis.
- `VALidate` checks a prescription file for structural errors.

### Reading Part I entries

Each entry gives the command name (with minimum match), prerequisite
tags, a one-line effect, and — where expanded — the dialog, behavior
notes, and related commands.  Every entry applies equally to SMACOS
(next chapter): the command string is the same; only the argument
delivery differs.
