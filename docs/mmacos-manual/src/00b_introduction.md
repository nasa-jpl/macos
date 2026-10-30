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

### What to expect from the assistant

A good assistant reads before it writes, runs a template before it writes a script, and reports the
numbers it measured rather than the numbers it expected. When it changes an optical rule it reads
the reference first; when a check fails it says so, with the output. The requests in this manual
were answered that way, and the reading lists under *The assistant reads* are the shortest route
to getting the same from yours.
