C***********************************************************************
C  MODULE: plot_mode_mod
C
C  Session-only display-convention state for Giza/PGPlot raster plots.
C
C  ifAstroMode = .TRUE.  (default)  : astronomy / negative convention —
C                                     large values render dark, small
C                                     values render bright. Matches the
C                                     legacy PGPlot default MACOS users
C                                     are accustomed to.
C              = .FALSE.              : conventional / positive — large
C                                     values render bright. Matches the
C                                     default Giza GRAY ramp.
C
C  Toggled at runtime via the IMGMODE command (POS | NEG | SHOW).
C  Applies to raster displays (OPD, INT, PIX, etc. via GRAY/PGIMAG).
C  Non-raster plots (SLICE, WIRE, CONTOUR, COLUMN, ROW, SPOT) always
C  render black ink on white background regardless of this flag.
C
C  Annotation strings are set by the command handler before the draw
C  call; the draw routine emits them at the bottom of the viewport
C  via pgmtxt('B', ...) and clears them on exit so they don't bleed
C  into the next plot.
C***********************************************************************

	MODULE plot_mode_mod
	  IMPLICIT NONE

	  LOGICAL :: ifAstroMode = .TRUE.

C  Bottom-of-plot annotation: one or two lines drawn by the raster plot
C  routine after PGLABEL.  Caller sets these, draw routine clears them.
	  LOGICAL            :: annotOn    = .FALSE.
	  CHARACTER(LEN=128) :: annotLine1 = ' '
	  CHARACTER(LEN=128) :: annotLine2 = ' '

C  Color-wedge label printed below the wedge by GRAY.  Caller sets it
C  per-command (e.g., 'OPD (m)', 'Intensity', 'Pixel value'); cleared
C  back to the generic default after the draw.
	  CHARACTER(LEN=64) :: wedgeLabel = 'pixel value'

	CONTAINS

C  Clear annotation + wedge label so the next plot command starts fresh.
	  SUBROUTINE ClearAnnot()
	    annotOn    = .FALSE.
	    annotLine1 = ' '
	    annotLine2 = ' '
	    wedgeLabel = 'pixel value'
	  END SUBROUTINE ClearAnnot

	END MODULE plot_mode_mod
