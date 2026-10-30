## APPENDIX A: Examples and Tutorial

*Demonstrations*

This appendix contains prescriptions and macros that illustrate various aspects of MACOS usage. Electronic versions of these files are in `docs/macos-manual/examples/` of the MACOS repository, and each listing below is that file as it stands there. The interested user can execute these examples to gain experience using MACOS.

### Cassegrain Telescope Examples

#### Cassegrain.jou

<!-- transcript: listing_Cassegrain_jou -->

#### Cassegrain.in

<!-- transcript: listing_Cassegrain_in -->

#### CassWithExitPupil.jou

<!-- transcript: listing_CassWithExitPupil_jou -->

#### CassWithExitPupil.in

<!-- transcript: listing_CassWithExitPupil_in -->

### Diffractive Elements Examples

#### GratingExample.jou

<!-- transcript: listing_GratingExample_jou -->

#### GratingExample.in

<!-- transcript: listing_GratingExample_in -->

#### HOEExample.jou

<!-- transcript: listing_HOEExample_jou -->

#### HOEExample.in

<!-- transcript: listing_HOEExample_in -->

### Non-Sequential Surface Examples

#### CornerCube.jou

<!-- transcript: listing_CornerCube_jou -->

#### CornerCube.in

<!-- transcript: listing_CornerCube_in -->

#### Luneberg.jou

<!-- transcript: listing_Luneberg_jou -->

#### Luneberg.in

<!-- transcript: listing_Luneberg_in -->

### Adaptive Optics Example

#### AOExample.jou

<!-- transcript: listing_AOExample_jou -->

#### AOExample.in

<!-- transcript: listing_AOExample_in -->

#### AOmirror.test

<!-- transcript: listing_AOmirror_test -->

### Segmented Telescope Example

#### SegDemo.jou

<!-- transcript: listing_SegDemo_jou -->

#### SegDemo.in

<!-- transcript: listing_SegDemo_in -->

### Near-Field Propagation Example

#### CoroExample.jou

<!-- transcript: listing_CoroExample_jou -->

#### CoroExample.in

<!-- transcript: listing_CoroExample_in -->

### Image Simulation Example

#### ImageDemo.jou

<!-- transcript: listing_ImageDemo_jou -->

#### ImageDemo.in

<!-- transcript: listing_ImageDemo_in -->

### S-MACOS Example

#### smacos_dvr.F

The SMACOS test driver built by `makesd.sh` (`macos_f90/smacos_dvr.F`): it opens the graphics
device, loads `ZGD_test_files/SegDemo3.in` with `OLD`, and calls `OPD` and `SPOT` through the
`SMACOS` subroutine (Section 9).

<!-- transcript: listing_smacos_dvr_F -->
