"""
Phase: rewrite the obsolete '9.x Common Blocks and Include Files'
sub-section of SECTION 9 (SMACOS).

The original 3.2 manual described data sharing between SMACOS and a
calling program via Fortran 77 COMMON blocks (smacosvars.cmn,
elt.cmn, src.cmn, param.cmn).  Since v3.2 those COMMONs were replaced
with Fortran 90 modules; the canonical pattern lives in
MACOS_resources/GMI/GMI.F (Generic MACOS Interface).

This patch:
  * Renames the H4 heading to 'Sharing Data via Fortran Modules'.
  * Deletes the original sub-section body (smacosvars.cmn, elt.cmn,
    src.cmn dumps and surrounding prose) up to the next heading.
  * Inserts new BodyText / CodeBlock paragraphs that explain:
      - which modules SMACOS exposes,
      - how to import only what you need (`use ..., only : ...`),
      - the macos_init_all(modelsize) initialization requirement,
      - a complete minimal example modeled on GMI.F.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_CODE, P_LIST
from docx_helpers import (iter_paragraphs, get_style, get_text, wtag,
                          find_para_parent, replace_para_text, make_para)


HEAD = "Sharing Data via Fortran Modules"

INTRO = (
    "SMACOS shares data with the calling program through a small set "
    "of Fortran 90 modules.  This replaces the COMMON-block / .cmn "
    "include-file mechanism used in MACOS 2.x; calling programs only "
    "need a list of `use` statements -- no manual array sizing, no "
    "shared-memory layout to keep in sync, no preprocessor model-size "
    "switches.  The canonical reference implementation is "
    "MACOS_resources/GMI/GMI.F (the MATLAB MEX interface), which "
    "exercises every module described here.")

MODULES_HEAD = "The SMACOS modules"

MODULES_INTRO = (
    "A user program that drives SMACOS typically opens with this "
    "block:")

USE_BLOCK = (
    "        use param_mod        ! compile-time array maxima (mElt, mRay, ...)\n"
    "        use src_mod          ! source / chief ray / wavelength\n"
    "        use elt_mod          ! every per-element array\n"
    "        use macos_mod        ! global flags and model state\n"
    "        use math_mod         ! shared math helpers (DMPROD, DOT, ...)\n"
    "        use sourcsub_mod     ! source-grid setup (SetSourceRayGrid, etc.)\n"
    "        use smacos_mod       ! the SMACOS dispatcher itself\n"
    "        ! optional, for advanced uses:\n"
    "        use traceutil_mod, only : CRIncidDir, opdRayMask\n"
    "        use design_optim_mod ! optimization (AVAR/MVAR/CALib)\n"
    "        implicit none")

MODULES_NOTE = (
    "Most callers only need the first seven; the others are pulled in "
    "as needed.  The `only :` form is preferred when you reference "
    "just a handful of names from a large module, both for compile "
    "speed and to make the dependency explicit.")

WHAT_HEAD = "What's in each module"

WHAT_BULLETS = [
    "param_mod -- compile-time array bounds: mElt (max elements), "
    "mRay (max rays), mPix, mDP, mGridMat, etc.  Every other module "
    "uses these to declare allocatable arrays.",
    "src_mod -- source description: ChfRayPos(3), ChfRayDir(3), "
    "Wavelen, Flux, Aperture, GridType, SegCoord(3,mElt), and the "
    "wavefront ray pos / dir grids RayPos(3,mRay) / RayDir(3,mRay).",
    "elt_mod -- per-element data: KrElt(mElt), KcElt(mElt), "
    "psiElt(3,mElt), VptElt, RptElt, MonCoef(mMonCoef,mElt), "
    "FFCoef(mFFCoef,mElt), pFF / xFF / yFF / zFF, GridMat, plus "
    "EltID(mElt) and SrfType(mElt) flags.",
    "macos_mod -- session state: ifPlot, ifGrid, ifBuild, "
    "macos_realloc, plus the macos_init_all(modelsize) entry "
    "point that allocates / re-allocates everything.",
    "smacos_mod -- the SMACOS subroutine plus its dispatch helpers.  "
    "`call SMACOS('OLD', CARG, DARG, IARG, LARG, RARG, OPDMat, "
    "RaySpot, RMSWFE, PixArray)` is the call signature.",
    "math_mod -- helpers shared with macos_f90 internals: DMPROD "
    "(matrix product), DOT, DXPROD, DUNITIZE, DORTHOGANALIZE, etc.",
]

INIT_HEAD = "Initializing the runtime"

INIT_INTRO = (
    "Module-scope arrays are ALLOCATABLE; they have zero size until "
    "the runtime is initialized for a chosen model size.  Call "
    "macos_init_all once before the first SMACOS call:")

INIT_CODE = (
    "        integer, parameter :: modelsize = 256\n"
    "        call macos_init_all(modelsize)   ! sizes mElt, mRay, mPix, ...\n"
    "        call dopt_init                   ! sizes optimization arrays\n"
    "        macos_init = .true.")

INIT_NOTE = (
    "Valid model sizes are 128, 256, 512, 1024, 2048, 4096, and 8192.  "
    "Re-call macos_init_all with a different size to grow or shrink "
    "the runtime; existing arrays are deallocated and replaced.  In "
    "interactive macos this is the MREset command; in SMACOS-driven "
    "code it's a single subroutine call.")

EXAMPLE_HEAD = "A minimal SMACOS driver"

EXAMPLE_INTRO = (
    "The following is a complete program that loads a prescription, "
    "traces, and prints the OPD RMS.  GMI.F is essentially a longer "
    "version of this same pattern wrapped behind a MATLAB MEX "
    "gateway.")

EXAMPLE_CODE = (
    "      program smacos_driver\n"
    "        use param_mod\n"
    "        use src_mod,    only : ChfRayPos, ChfRayDir, Wavelen\n"
    "        use elt_mod,    only : RptElt, KrElt, EltName, nElt\n"
    "        use macos_mod,  only : macos_init\n"
    "        use sourcsub_mod, only : SetSourceRayGrid\n"
    "        use smacos_mod, only : SMACOS\n"
    "        implicit none\n"
    "\n"
    "        ! SMACOS argument-list buffers\n"
    "        character(len=132) :: command\n"
    "        character(len=32)  :: CARG(9)\n"
    "        real*8             :: DARG(9), RMSWFE\n"
    "        integer            :: IARG(9)\n"
    "        logical            :: LARG\n"
    "        real*4             :: RARG(9)\n"
    "        real*8, allocatable :: OPDMat(:,:), RaySpot(:,:)\n"
    "        real*4, allocatable :: PixArray(:,:)\n"
    "\n"
    "        integer :: modelsize\n"
    "\n"
    "        ! 1. initialize the runtime\n"
    "        modelsize = 256\n"
    "        call macos_init_all(modelsize)\n"
    "        macos_init = .true.\n"
    "        allocate(OPDMat(mpts,mpts), RaySpot(mRay,2),\n"
    "     &           PixArray(mPix,mPix))\n"
    "\n"
    "        ! 2. load a prescription\n"
    "        command  = 'OLD'\n"
    "        CARG(1)  = 'mySystem.in'\n"
    "        call SMACOS(command, CARG, DARG, IARG, LARG, RARG,\n"
    "     &              OPDMat, RaySpot, RMSWFE, PixArray)\n"
    "\n"
    "        ! 3. set up source ray grid (uses src_mod state)\n"
    "        call SetSourceRayGrid(npts, Dicr, Djcr, dxSource, .false.,\n"
    "     &                        .false.)\n"
    "\n"
    "        ! 4. build linear model + trace\n"
    "        command = 'BUILD'\n"
    "        IARG(1) = nElt    ! terminal element\n"
    "        call SMACOS(command, CARG, DARG, IARG, LARG, RARG,\n"
    "     &              OPDMat, RaySpot, RMSWFE, PixArray)\n"
    "\n"
    "        command = 'OPD'\n"
    "        IARG(1) = nElt\n"
    "        call SMACOS(command, CARG, DARG, IARG, LARG, RARG,\n"
    "     &              OPDMat, RaySpot, RMSWFE, PixArray)\n"
    "\n"
    "        ! 5. read results back from module-shared state\n"
    "        write(*,*) 'OPD RMS =', RMSWFE\n"
    "        write(*,*) 'KrElt(1) =', KrElt(1),\n"
    "     &             '  RptElt(:,1) =', RptElt(:,1)\n"
    "      end program smacos_driver")

PASSING_HEAD = "Reading and writing module data"

PASSING_BODY = (
    "Once `use elt_mod` is in scope, every element-array is directly "
    "addressable.  Reading is a simple subscript: `KrElt(3)`, "
    "`RptElt(:,iElt)`, `MonCoef(:,iElt)`.  Writing follows the same "
    "pattern -- but if a write changes anything that affects the ray "
    "trace (frame geometry, surface coefficients, etc.), the linear "
    "model becomes stale and you must re-run BUILD before exercising "
    "it.  GMI.F's ApplyPerturbationToOpticalSystem is a worked "
    "example: it mutates ZernCoef / MonZernCoef / TElt directly via "
    "the module, then issues `command = 'PERTURB'` for the dispatcher "
    "to recompute frames.")

NOMINAL_HEAD = "Nominal-state save/restore"

NOMINAL_BODY = (
    "Programs that apply repeated perturbations (sensitivity loops, "
    "tolerancing) need to restore the optical state between "
    "iterations.  GMI.F shows the pattern: snapshot every per-element "
    "array that PERTURB modifies into a local *Nom array on first "
    "entry, then copy back at the start of every later call.  For "
    "FreeForm surfaces this includes pFF/xFF/yFF/zFF and "
    "pData/xData/yData/zData; SMACOS does not snapshot these "
    "automatically.  Do not store snapshots in a Fortran 77 COMMON "
    "block; allocate them in the calling module instead.")

LEGACY_HEAD = "Legacy COMMON-block support"

LEGACY_BODY = (
    "The old smacosvars.cmn / elt.cmn / src.cmn / param.cmn include "
    "files were retired in v3.2.  Calling programs that still depend "
    "on them must be rewritten to `use` the modules above.  In "
    "almost every case this is a one-for-one mapping: `KrElt(i)` "
    "from elt.cmn becomes `use elt_mod, only : KrElt` and the same "
    "subscript expression.  See git history of macos_f90/elt_mod.F "
    "for the full alias map.")


def find_heading(root, level: str, text_substr: str):
    """Return first <w:p> with style=level whose text contains substr."""
    for p in iter_paragraphs(root):
        if get_style(p) == level and text_substr in get_text(p):
            return p
    return None


def transform(root):
    h = find_heading(root, "Heading4", "Common Blocks and Include Files")
    if h is None:
        # Already migrated; idempotent.
        return {"already migrated": 1}

    # Rename the heading.
    replace_para_text(h, HEAD, style="Heading4")

    # Delete every paragraph after h up to the next non-empty Heading.
    parent = find_para_parent(root, h)
    children = list(parent)
    idx = children.index(h) + 1
    deleted = 0
    while idx < len(children):
        c = children[idx]
        if c.tag == wtag("p"):
            cs = get_style(c)
            ct = get_text(c).strip()
            if cs in ("Heading2", "Heading3", "Heading4") and ct:
                break
        parent.remove(c)
        deleted += 1
        children = list(parent)

    # Build the replacement paragraph stream.
    new_paras = [
        make_para(P_BODY, INTRO),
        make_para("Heading4", MODULES_HEAD),
        make_para(P_BODY, MODULES_INTRO),
        make_para(P_CODE, USE_BLOCK),
        make_para(P_BODY, MODULES_NOTE),
        make_para("Heading4", WHAT_HEAD),
    ]
    new_paras.extend(make_para(P_LIST, b) for b in WHAT_BULLETS)
    new_paras.extend([
        make_para("Heading4", INIT_HEAD),
        make_para(P_BODY, INIT_INTRO),
        make_para(P_CODE, INIT_CODE),
        make_para(P_BODY, INIT_NOTE),
        make_para("Heading4", EXAMPLE_HEAD),
        make_para(P_BODY, EXAMPLE_INTRO),
        make_para(P_CODE, EXAMPLE_CODE),
        make_para("Heading4", PASSING_HEAD),
        make_para(P_BODY, PASSING_BODY),
        make_para("Heading4", NOMINAL_HEAD),
        make_para(P_BODY, NOMINAL_BODY),
        make_para("Heading4", LEGACY_HEAD),
        make_para(P_BODY, LEGACY_BODY),
    ])

    # Insert immediately after h.
    children = list(parent)
    h_idx = children.index(h)
    for offset, np in enumerate(new_paras, start=1):
        parent.insert(h_idx + offset, np)

    return {
        "old paragraphs deleted": deleted,
        "new paragraphs inserted": len(new_paras),
    }


if __name__ == "__main__":
    apply(__file__, transform)
