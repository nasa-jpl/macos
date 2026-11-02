::: {custom-style="Title"}
mmacos Manual
:::

```{=openxml}
<w:p><w:pPr><w:spacing w:before="6600" w:after="0"/></w:pPr></w:p>
```

::: {custom-style="Subtitle"}
MATLAB and MACOS with an AI Assistant
:::

::: {custom-style="Subtitle"}
Draft --- October 2026
:::

::: {custom-style="Subtitle"}
Jet Propulsion Laboratory, California Institute of Technology
:::

\newpage

#### Contact

For information on MACOS contact:

Scott Basinger <Scott.A.Basinger@jpl.nasa.gov>

Jet Propulsion Laboratory 306-438, 4800 Oak Grove Drive,
Pasadena CA 91109

#### Authorship

David C. Redding, Scott Basinger, John Lou, Luis Marchen, Laura Needels, Carl Nissly, Norbert Sigrist and Sam Sirlin

#### Acknowledgment

The MACOS program was developed at the Jet Propulsion Laboratory, California Institute of Technology.

MACOS software and documentation are copyright 1992-2026 by the California Institute of Technology, Jet Propulsion Laboratory.

U.S. Government sponsorship of parts of MACOS under NASA Contract NAS7-918 is acknowledged.

#### Installation

*mmacos* is the MATLAB interface to the MACOS engine. It lives in the `MACOS_resources` repository beside the engine's `macos` repository; build the gfortran engine first (the MACOS Manual, Installation: `source ./makems.sh release gfortran`), then the interface (`mmacos/README.md`, Build):

    cd MACOS_resources/mmacos
    make                                  # the mex, linked against ../macos/build_release_gfortran
    ./run_mmacos_tests.sh fast            # the fast suite

and in MATLAB, once per session:

    run('<path-to>/MACOS_resources/mmacos/mmacos_setup.m');

#### How this manual is made

Every listing, figure and quoted number in this manual is generated: each chapter has a driver (`docs/mmacos-manual/tools/transcripts/<chapter>.m`) that runs the template the chapter narrates, unchanged, in a scratch directory, and `make transcripts` writes what it printed and drew into the chapter. A quoted number is pinned (five significant figures on a signal, a bound on a numerical zero), so a change in the engine or a template that moves a number stops the build until the number is re-pinned on purpose. Each chapter speaks in three voices: *You ask* (the request a user gives the assistant), *The assistant reads* (the files it reads before it acts, as paths you can hand to any assistant) and *It runs* (the template, its output, its figures, its numbers).

\newpage

::: {custom-style="TOC Heading"}
Contents
:::

```{=openxml}
<w:p><w:r><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> TOC \o "1-3" \h \z \u </w:instrText></w:r><w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>(contents -- update the field)</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r></w:p>
```
