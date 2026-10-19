::: {custom-style="Title"}
MACOS Manual
:::

::: {custom-style="Subtitle"}
Modeling and Analysis for Controlled Optical Systems
:::

::: {custom-style="Subtitle"}
Version 4.2 (beta) --- October 2026
:::

::: {custom-style="Subtitle"}
Jet Propulsion Laboratory, California Institute of Technology
:::

\newpage

::: {custom-style="TOC Heading"}
Contents
:::

```{=openxml}
<w:p><w:r><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> TOC \o "1-3" \h \z \u </w:instrText></w:r><w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>(contents -- update the field)</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r></w:p>
```

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

MACOS is distributed as two git repositories on GitHub: the engine and this manual in `macos`, and the MATLAB and Python bindings, the templates and the design tools in `MACOS_resources`. Clone them side by side (the build scripts expect `MACOS_resources` next to `macos`):

    git clone git@github.com:nasa-jpl/macos.git
    git clone git@github.com:nasa-jpl/MACOS_resources.git
    cd macos
    source ./makems.sh release            # the engine and the macos command-line program (ifx)
    source ./makems.sh release gfortran   # the same with gfortran (the build the MATLAB interface links)

`makems.sh` builds the engine library and the interactive program into `build_release[_gfortran]/bin/macos`; `makeall.sh` adds the SMACOS driver and the GMI MATLAB interface. The requirements (CMake, a Fortran compiler -- Intel ifx or gfortran -- and gcc), the Linux, macOS and Windows notes and the debug builds are in `BUILD.md` at the top of the `macos` repository; the branch model is in `BRANCHES.md`. The MATLAB interface (`MACOS_resources/mmacos`, built by its `run_mmacos_tests.sh`) and the Python package (`MACOS_resources/pymacos`) each carry their own README. Access to the repositories is by invitation; write to the contact above.
