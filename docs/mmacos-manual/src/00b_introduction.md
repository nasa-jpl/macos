## INTRODUCTION

This manual is written for an optical engineer who works with an AI assistant. The assistant types;
you decide. Every example in the book is a request you could give the assistant today, followed by
what a competent assistant does with it and what came out when it did. Nothing here is a mock-up:
the listings, figures and numbers were produced by running the library, and the build stops if a
number moves (the front matter, *How this manual is made*).

### The three voices

Each worked example is told in three voices, always in the same order.

- ***You ask*** is the request, in the words a user would type. It carries the numbers that matter
  (aperture, f-number, field, wavelength, a contrast floor) and nothing about how to get there.
  Treat it as a template for your own request: change the numbers, keep the shape.
- ***The assistant reads*** is the list of files the assistant opens before it acts. They are paths
  in the two repositories (`macos`, `MACOS_resources`), so the list is portable: hand it to any
  assistant and it starts from the same ground. The list is also the manual's answer to "where is
  that written down?"
- ***It runs*** is what the assistant ran, unchanged, and what it printed and drew. The gray
  listings are the console output of that run; the figures are the files it wrote; a number in the
  prose is pinned to the run that produced it.

### How to use the manual

- **To learn the library, read Part I in order.** Its eight sections follow one 4 m telescope from
  a sketch to an error budget, and each section consumes the one before, so the thread shows how
  the pieces fit: the design layer, the instruments, segmentation, sensitivities, metrology, the
  control loop, the budget, and the coronagraph.
- **To get something done, go to Part II.** It is the template catalog, one page per template in
  the order of `templates/00_INDEX.md`. Find the page nearest to what you need, paste its closing
  request to your assistant with your numbers in it, and compare what comes back with the page.
- **To check what the assistant did, use the pinned numbers.** A number printed in this manual is
  what the engine measured when the page was built. If your run prints a different one on the
  same inputs, something changed, and the assistant should be able to say what.
- **To find the rules, go to Part III.** It points to the references that live beside the code,
  the design rules and agent guides, the library's own record, and the engine's manual, and it has
  the one chapter written for the assistant rather than for you.

### Getting started, step by step

If you have not worked with an AI assistant before, this is the mechanical part. It takes an
afternoon the first time and a minute after that.

1. **Get the code.** Clone the two repositories side by side, so that `macos/` and
   `MACOS_resources/` share a parent directory; the build scripts and the library find each other
   by that relative path. Your administrator has the repository addresses and the branch to use.
2. **Build the engine, then the interface.** Follow *Installation* in the front matter: build the
   gfortran engine in `macos/`, then `make` in `MACOS_resources/mmacos/`, then run the fast test
   suite once. A green suite is your proof that the installation works before the assistant
   touches anything.
3. **Start MATLAB and run the setup script** (`mmacos_setup.m`, once per session). Type
   `macos.init(128)` and then `help macos.trace`; if both answer, MATLAB sees the engine.
4. **Open a terminal in the parent directory** of the two repositories and start your assistant
   there. With Claude Code, that is the command `claude`; in VS Code, open the parent directory as
   the workspace and open the assistant's panel; other assistants have their own command. The
   directory you start in matters: the assistant reads the `CLAUDE.md` files it finds there
   (`macos/CLAUDE.md`, `macos/macos_f90/CLAUDE.md`, `MACOS_resources/mmacos/CLAUDE.md`), and those
   files carry the working rules that this manual's examples relied on.
5. **Point the assistant at this manual.** Its sources are
   `macos/docs/mmacos-manual/src/*.md`, one file per section, and the built copy is
   `macos/docs/mmacos-manual/build/mmacosMan.docx`. A first request can be as plain as: *Read
   docs/mmacos-manual/src/00b_introduction.md and src/01_telescope.md, then tell me what you would
   do with the request in Section 1.*
6. **Make your first request from a page in this book.** Copy a *You ask* paragraph, put your own
   numbers in it, and send it. Ask the assistant to list the files it will read before it acts; if
   that list looks like the page's *The assistant reads*, you are on the same footing the manual
   was written from.
7. **Check the result against the page.** Compare what the assistant printed and drew with the
   page's *It runs*. The same template on the same inputs prints the same pinned numbers; a
   difference is something to ask about, not something to accept.
8. **Keep MATLAB open.** The assistant runs MATLAB for you in batch mode, but the figures, the
   `.mat` workspaces and the `.in` prescriptions it writes are ordinary files beside the template;
   open them yourself whenever you want to look for yourself.

### What to expect from the assistant

A good assistant reads before it writes, runs a template before it writes a script, and reports the
numbers it measured rather than the numbers it expected. When it changes an optical rule it reads
the reference first; when a check fails it says so, with the output. The requests in this manual
were answered that way, and the reading lists under *The assistant reads* are the shortest route
to getting the same from yours.
