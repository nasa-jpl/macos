## APPENDIX A: Examples and Tutorial

*Demonstrations*

This appendix contains prescriptions and macros that illustrate various aspects of MACOS usage. Electronic versions of these files have been provided with the MACOS distribution tapes. The interested user can execute these examples to gain experience using MACOS.

### Cassegrain Telescope Examples

#### Cassegrain.jou

    % This macro generates example figures
    % for Section 2 of the MACOS User Manual
    old Cassegrain modify ngridpts=51 q
    spot 1 Tout % Fig 3a
    draw;;;; % Fig 3b
    spot 5 Tout % Fig 4a
    perturb 2 global 0,1e-4,0 0,0,0
    spot 5 Tout % Fig 4b

#### Cassegrain.in

    % Cassegrain telescope example
    ChfRayDir= 0d0 0d0 1d0 ChfRayPos= 0d0 0d0 -5d0 zSource= 1.d+21
    IndRef= 1d0
    Extinc= 0d0 Wavelen= 1.d-06
    Flux= 1d0 GridType= Circular Aperture= 4d0 Obscratn= 0d0 nGridpts= 15

xGrid= -1d0 0d0 0d0 yGrid= 0d0 -1d0 0d0 nElt= 5

iElt= 1

EltName= SecMirObs Element= Obscuring Surface= Conic

KrElt= -1.d+20

[Modeling and Analysis for Controlled Optical Systems 167](#modeling-and-analysis-for-controlled-optical-systems)

    KcElt= 0d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 -4.2d+00 RptElt= 0d0 0d0 -4.2d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 3
    ObsType= Circle
    ObsVec= 5.d-01 0d0 0d0 ObsType= Rectangle
    ObsVec= -0.003125 0.003125 -3 3
    ObsType= Rectangle
    ObsVec= -3 3 -0.003125 0.003125
    xObs= 1d0 0d0 0d0 ApType= Circular ApVec= 2d0 0d0 0d0
    zElt= 1.d+20 PropType= Geometric nECoord= -6
    iElt= 2
    EltName= Primary Element= Reflector Surface= Conic
    KrElt= -1.08d+01
    KcElt= -1d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 0d0 RptElt= 0d0 0d0 0d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= Circular
    ApVec= 2.1d0 0d0 0d0 zElt= 5.4d+00
    PropType= Geometric nECoord= -6
    iElt= 3
    EltName= Secondary Element= Reflector Surface= Conic
    KrElt= -3.526787501D+00 KcElt= -2.670556022D+00
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 -4.061145902D+00
    RptElt= 0d0 0d0 -5.4d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1.338854098D+00
    PropType= Geometric nECoord= -6
    iElt= 4
    EltName= Ref_surf Element= Reference Surface= Conic
    KrElt= -5.560145902D+00
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -4.060145902D+00 RptElt= 0d0 0d0 -4.060145902D+00

+--------------+--------------------------------+---------------------+
| IndRef=      | 1d0 1d22                       |                     |
| Extinc=      |                                |                     |
| nObs=        | 0                              |                     |
| ApType=      |                                |                     |
| zElt=        | None 5.560145902D+00           |                     |
| PropType=    |                                |                     |
| nECoord=     | FarField                       |                     |
|              |                                |                     |
|              | -6                             |                     |
+==============+================================+=====================+
| iElt=        | 5                              |                     |
| EltName=     |                                |                     |
| Element=     | Detector FocalPlane Flat       |                     |
| Surface=     |                                |                     |
| KrElt=       | -1d22                          |                     |
| KcElt=       |                                |                     |
| psiElt=      | 0d0                            |                     |
| VptElt=      |                                |                     |
| RptElt=      | 0d0 0d0 1d0 0d0 0d0 1.5d+00    |                     |
| IndRef=      |                                |                     |
| Extinc=      | 0d0 0d0 1.5d+00                |                     |
| nObs=        |                                |                     |
| ApType=      | 1d0 1d22                       |                     |
| zElt=        |                                |                     |
| PropType=    | 0                              |                     |
| nECoord=     |                                |                     |
|              | None 0d0                       |                     |
|              |                                |                     |
|              | Geometric                      |                     |
|              |                                |                     |
|              | -6                             |                     |
+--------------+--------------------------------+---------------------+
| nOutCord=    | 5                              | 0d0 0d0 0d0         |
|              |                                |                     |
| Tout=        | 1d0 0d0 0d0 0d0                |                     |
+--------------+--------------------------------+---------------------+
|              | 0d0 1d0 0d0 0d0 0d0 0d0 0d0    | 0d0 0d0 0d0 0d0 0d0 |
|              | 1d0 0d0 0d0 0d0 0d0            | 0d0 1d0 0d0 0d0     |
|              |                                |                     |
|              | 0d0 0d0 0d0 0d0                | 0d0 0d0 1d0         |
+--------------+--------------------------------+---------------------+

#### CassWithExitPupil.jou

    % This macro generates more example figures
    % for Section 2 of the MACOS User Manual
    old CassWithExitPupil modify ngridpts=128 q
    perturb 3 global 0,0,0 0,0,1e-4
    stop obj 0,0,0
    fex 5 yes
    opd 5 Tout % Fig 6a
    perturb 2 global 0,1e-4,0 0,0,0
    opd 5 Tout % Fig 6b
    old CassWithExitPupil modify ngridpts=64 q stop obj 0,0,0
    fex 5 yes
    intensity 6 % Fig 10a
    perturb 2 global 0,1e-4,0 0,0,0
    fex 5 yes
    intensity 6 % Fig 10b

#### CassWithExitPupil.in

    % Cassegrain telescope example
    % This version has an exit pupil
    ChfRayDir= 0d0 0d0 1d0 ChfRayPos= 0d0 0d0 -5d0 zSource= 1.d+21
    IndRef= 1d0
    Extinc= 0d0 Wavelen= 1.d-06
    Flux= 1d0 GridType= Circular Aperture= 4d0 Obscratn= 0d0 nGridpts= 15
    xGrid= -1d0 0d0 0d0 yGrid= 0d0 -1d0 0d0 nElt= 6
    iElt= 1
    EltName= SecMirObs Element= Obscuring Surface= Conic
    KrElt= -1.d+20 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 -4.2d+00 RptElt= 0d0 0d0 -4.2d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 3
    ObsType= Circle
    ObsVec= 5.d-01 0d0 0d0 ObsType= Rectangle
    ObsVec= -0.003125 0.003125 -3 3
    ObsType= Rectangle
    ObsVec= -3 3 -0.003125 0.003125
    xObs= 1d0 0d0 0d0 ApType= Circular ApVec= 2d0 0d0 0d0
    zElt= 1.d+20 PropType= Geometric nECoord= -6
    iElt= 2
    EltName= Primary Element= Reflector Surface= Conic
    KrElt= -1.08d+01
    KcElt= -1d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 0d0 RptElt= 0d0 0d0 0d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= Circular
    ApVec= 2.1d0 0d0 0d0 zElt= 5.4d+00
    PropType= Geometric nECoord= -6
    iElt= 3
    EltName= Secondary Element= Reflector Surface= Conic
    KrElt= -3.526787501D+00 KcElt= -2.670556022D+00
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 -4.061145902D+00
    RptElt= 0d0 0d0 -5.4d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1.338854098D+00
    PropType= Geometric nECoord= -6
    iElt= 4
    EltName= Return1 Element= Return Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 1.5d+00 RptElt= 0d0 0d0 1.5d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 0d0 PropType= Geometric nECoord= -6
    iElt= 5
    EltName= ExitPupil Element= Return Surface= Conic
    KrElt= -5.560145902D+00
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -4.060145902D+00 RptElt= 0d0 0d0 -4.060145902D+00
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 5.560145902D+00
    PropType= FarField nECoord= -6
    iElt= 6
    EltName= Detector Element= FocalPlane Surface= Flat
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 1.5d+00 RptElt= 0d0 0d0 1.5d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 0d0 PropType= Geometric nECoord= -6
    nOutCord= 5
    Tout= 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0

### Diffractive Elements Examples

#### GratingExample.jou

    % This macro generates example figures
    % for Section 4 of the MACOS User Manual old gratingExample
    mod OrderHOE(1)=0 q draw;;;;
    mod OrderHOE(1)=1 q
    draw;;;; % Fig ??
    mod OrderHOE(1)=2 q draw;;;;
    mod OrderHOE(1)=-1 q draw;;;;
    mod OrderHOE(1)=-2 q
    draw;;;; % Fig ??

#### GratingExample.in

    ChfRayDir= 0d0 0d0 1d0 ChfRayPos= 0d0 0d0 -6d0 zSource= 1d22
    IndRef= 1d0
    Extinc= 0d0 Wavelen= 6.d-07
    Flux= 1d0 GridType= Circular Aperture= 1d0 Obscratn= 0d0 nGridpts= 65
    xGrid= 1d0 0d0 0d0 yGrid= 0d0 1d0 0d0 nElt= 2
    iElt= 1
    EltName= Grating Element= Grating Surface= Flat
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 0d0 RptElt= 0d0 0d0 0d0

+------------------------------------------------+--------------+------+
| IndRef= 1d0                                    |              |      |
|                                                |              |      |
| Extinc= 1d22                                   |              |      |
|                                                |              |      |
| h1HOE= 1d0 0d0 0d0                             |              |      |
|                                                |              |      |
| OrderHOE= 1d0                                  |              |      |
|                                                |              |      |
| RuleWidth= 18.97366596d-7 nObs= 0              |              |      |
|                                                |              |      |
| ApType= None                                   |              |      |
|                                                |              |      |
| zElt= 1d22 PropType= Geometric nECoord= -6     |              |      |
+================================================+==============+======+
| iElt= 2                                        |              |      |
|                                                |              |      |
| EltName= FP Element= FocalPlane Surface= Flat  |              |      |
|                                                |              |      |
| KrElt= 1d22                                    |              |      |
|                                                |              |      |
| KcElt= 0d0                                     |              |      |
|                                                |              |      |
| psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 -3d0       |              |      |
| RptElt= 0d0 0d0 0d0 IndRef= 1d0                |              |      |
|                                                |              |      |
| Extinc= 0d0                                    |              |      |
|                                                |              |      |
| nObs= 0                                        |              |      |
|                                                |              |      |
| ApType= None                                   |              |      |
|                                                |              |      |
| zElt= 1d22 PropType= Geometric nECoord= -6     |              |      |
+------------------------------------------------+--------------+------+
| nOutCord= 5                                    | 0d0 0d0      | 0d0  |
|                                                |              |      |
| Tout= 1d0 0d0 0d0 0d0                          |              |      |
+------------------------------------------------+--------------+------+
| 0d0 1d0 0d0 0d0                                | 0d0 0d0      | 0d0  |
+------------------------------------------------+--------------+------+
| 0d0 0d0 0d0 1d0                                | 0d0 0d0      | 0d0  |
+------------------------------------------------+--------------+------+
| 0d0 0d0 0d0 0d0                                | 1d0 0d0      | 0d0  |
+------------------------------------------------+--------------+------+
| 0d0 0d0 0d0 0d0                                | 0d0 0d0      | 1d0  |
+------------------------------------------------+--------------+------+
| **A.2.3 HOEExample.jou**                       |              |      |
+------------------------------------------------+--------------+------+

    % This macro generates example figures
    % for Section 4 of the MACOS User Manual
    loa HOEExample mod ngridpts=15 q map
    ray 76
    stop elt 1 0,0 fex;;;
    int;;
    draw;;;; % Fig ??
    mod Element(1)=reflector q ray 76
    fex;;; int;;

+--------------------------------+-------------------------------------+
| draw;;;;                       | \% Fig ??                           |
+================================+=====================================+
| **A.2.4 HOEExample.in**        |                                     |
+--------------------------------+-------------------------------------+
| ChfRayDir= 0d0 0d0 1d0         |                                     |
| ChfRayPos= 0d0 0d0 -4d0        |                                     |
| zSource= 1d22                  |                                     |
|                                |                                     |
| IndRef= 1d0                    |                                     |
|                                |                                     |
| Extinc= 0d0 Wavelen= 6.d-07    |                                     |
|                                |                                     |
| Flux= 1d0 GridType= Circular   |                                     |
| Aperture= 5.d-01 Obscratn= 0d0 |                                     |
| nGridpts= 64                   |                                     |
|                                |                                     |
| xGrid= 1d0 0d0 0d0 yGrid= 0d0  |                                     |
| 1d0 0d0 nElt= 2                |                                     |
+--------------------------------+-------------------------------------+
| iElt= 1                        |                                     |
|                                |                                     |
| EltName= HOEonMirror Element=  |                                     |
| HOE Surface= Conic             |                                     |
|                                |                                     |
| KrElt= -6d0 KcElt= -1d0        |                                     |
|                                |                                     |
| psiElt= 0d0 0d0 -1d0 VptElt=   |                                     |
| 0d0 0d0 0d0 RptElt= 0d0 0d0    |                                     |
| 0d0 IndRef= 1d0                |                                     |
|                                |                                     |
| Extinc= 1d22                   |                                     |
|                                |                                     |
| h1HOE= 0d0 0d0 -1d22 h2HOE=    |                                     |
| 1d0 0d0 -3d0                   |                                     |
|                                |                                     |
| OrderHOE= 1d0 WaveHOE= 6.d-07  |                                     |
|                                |                                     |
| nObs= 0                        |                                     |
|                                |                                     |
| ApType= None                   |                                     |
|                                |                                     |
| zElt= 1d22 PropType= Geometric |                                     |
| nECoord= -6                    |                                     |
+--------------------------------+-------------------------------------+
| iElt= 2                        |                                     |
|                                |                                     |
| EltName= FP Element=           |                                     |
| FocalPlane Surface= Flat       |                                     |
|                                |                                     |
| KrElt= -1d22 KcElt= 0d0        |                                     |
|                                |                                     |
| psiElt= 0d0 0d0 1d0 VptElt=    |                                     |
| 0d0 0d0 -3d0 RptElt= 0d0 0d0   |                                     |
| 0d0 IndRef= 1d0                |                                     |
|                                |                                     |
| Extinc= 0d0                    |                                     |
|                                |                                     |
| nObs= 0                        |                                     |
|                                |                                     |
| ApType= None                   |                                     |
|                                |                                     |
| zElt= 1d22 PropType= Geometric |                                     |
| nECoord= -6                    |                                     |
+--------------------------------+-------------------------------------+
| nOutCord= 5                    | 0d0 0d0 0d0                         |
|                                |                                     |
| Tout= 1d0 0d0 0d0 0d0          |                                     |
+--------------------------------+-------------------------------------+
| 0d0 1d0 0d0 0d0                | 0d0 0d0 0d0                         |
+--------------------------------+-------------------------------------+
| 0d0 0d0 0d0 1d0                | 0d0 0d0 0d0                         |
+--------------------------------+-------------------------------------+
| 0d0 0d0 0d0 0d0                | 1d0 0d0 0d0                         |
+--------------------------------+-------------------------------------+
| 0d0 0d0 0d0 0d0                | 0d0 0d0 1d0                         |
+--------------------------------+-------------------------------------+

### Non-Sequential Surface Examples

#### CornerCube.jou

    old CornerCube draw;;;xz
    mod ChfRayPos(1)=0.762 ChfRayPos(2)=0.762 Aperture=0.125 nGridpts=11
    q draw;;;yz

#### CornerCube.in

    ChfRayDir= 0d0 0d0 -1d0 ChfRayPos= 1.d-05 1.d-05 1.d+01 zSource= 1.d22
    IndRef= 1d0
    Extinc= 0d0 Wavelen= 1.059d-04
    Flux= 1d0 GridType= Circular Aperture= 2.54d+00 Obscratn= 0d0 nGridpts= 65
    xGrid= 0d0 1d0 0d0 yGrid= -1d0 0d0 0d0 nElt= 5
    iElt= 1
    EltName= Ref Srf Element= Reference Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 3.5d+00 RptElt= 0d0 0d0 0d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 0d0 PropType= NFP1surf nECoord= 0
    iElt= 2
    EltName= CC1 Element= NSReflector Surface= Conic
    KrElt= 1d22
    KcElt= 0d0
    psiElt= 8.164994677D-01 0d0 5.773461867D-01
    VptElt= 0d0 0d0 0d0 RptElt= 0d0 0d0 0d0
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 3
    EltName= CC2 Element= NSReflector Surface= Conic
    KrElt= 1d22
    KcElt= 0d0
    psiElt= -4.082497338D-01 7.071092812D-01 5.773461867D-01
    VptElt= 0d0 0d0 0d0 RptElt= 0d0 0d0 0d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 4
    EltName= CC3 Element= NSReflector Surface= Conic
    KrElt= 1d22
    KcElt= 0d0
    psiElt= -4.082497338D-01 -7.071092812D-01 5.773461867D-01
    VptElt= 0d0 0d0 0d0 RptElt= 0d0 0d0 0d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 5
    EltName= FP Element= FocalPlane Surface= Flat
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 3.5d+00 RptElt= 0d0 0d0 0d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1.d+03 PropType= Geometric nECoord= 0
    nOutCord= 5
    Tout= 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0
    0d0 0d0 0d0 0d0 0d0 0d0 1d0

#### Luneberg.jou

    % Luneberg.jou: a MACOS example macro. This macro
    % demonstrates use of non-sequential refracting
    % elements, polarization, and the gain command.
    old Luneberg draw 0 16 XZ
    beam cos 8.7500000000D-01 1.5
    mod nGridPts=128 q int 16
    mod nGridPts=64 q gain 17
    col
    gain 17 129
    pol 1,0 0,0
    mod nGridPts=128 q int 16
    mod nGridPts=64 q gain 17
    col
    gain 17 129

#### Luneberg.in

    % Luneberg.in: a MACOS example .in-file. This
    % macro demonstrates use of non-sequential
    % refractors, polarization, and the gain command.
    ChfRayDir= 0d0 0d0 1d0
    ChfRayPos= 0d0 0d0 -2.121332605D-01
    zSource= -1.d-02 IndRef= 1d0
    Extinc= 0d0 Wavelen= 2.15d-02
    Flux= 1d0 GridType= Circular Aperture= 1.75d+00 Obscratn= 0d0 nGridpts= 33
    xGrid= 1d0 0d0 0d0 yGrid= 0d0 1d0 0d0 nElt= 17
    iElt= 1
    EltName= Shell1 Element= NSRefractor Surface= Conic
    KrElt= -2.0353d-01
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 -2.0353d-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.056692495D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 2
    EltName= Shell2 Element= NSRefractor Surface= Conic
    KrElt= -1.887578424D-01
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -1.887578424D-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.084805420D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 3
    EltName= Shell3 Element= NSRefractor Surface= Conic
    KrElt= -1.872644104D-01
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -1.872644104D-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.147803142D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 4
    EltName= Shell4 Element= NSRefractor Surface= Conic
    KrElt= -1.771655689D-01
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -1.771655689D-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.182543332D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 5
    EltName= Shell5 Element= NSRefractor Surface= Conic
    KrElt= -1.529604025D-01
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -1.529604025D-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.235415610D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 6
    EltName= Shell6 Element= NSRefractor Surface= Conic
    KrElt= -1.387929631D-01
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -1.387929631D-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.265593240D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 7
    EltName= Shell7 Element= NSRefractor Surface= Conic
    KrElt= -1.297999186D-01
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -1.297999186D-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.283513736D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 8
    EltName= Shell8 Element= NSRefractor Surface= Conic
    KrElt= -1.087287483D-01
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -1.087287483D-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.312911218D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 9
    EltName= Shell9 Element= NSRefractor Surface= Conic
    KrElt= -9.780754704D-02
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -9.780754704D-02
    RptElt= 0d0 0d0 0d0 IndRef= 1.340947910D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 10
    EltName= Shell10 Element= NSRefractor Surface= Conic
    KrElt= -8.777974770D-02
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -8.777974770D-02
    RptElt= 0d0 0d0 0d0 IndRef= 1.338212289D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 11
    EltName= Shell11 Element= NSRefractor Surface= Conic
    KrElt= -7.943573990D-02
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -7.943573990D-02
    RptElt= 0d0 0d0 0d0 IndRef= 1.371777076D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 12
    EltName= Shell12 Element= NSRefractor Surface= Conic
    KrElt= -6.183854602D-02
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -6.183854602D-02
    RptElt= 0d0 0d0 0d0 IndRef= 1.386006296D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 13 EltName= Shell13 Element= NSRefractor Surface= Conic
    KrElt= -4.663931773D-02
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -4.663931773D-02
    RptElt= 0d0 0d0 0d0 IndRef= 1.402911090D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 14
    EltName= Shell14 Element= NSRefractor Surface= Conic
    KrElt= -2.956340242D-02
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -2.956340242D-02
    RptElt= 0d0 0d0 0d0 IndRef= 1.412d+00 Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 15
    EltName= Aperture Element= Return Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 2.221359167D-01
    RptElt= 0d0 0d0 0d0 IndRef= 1.0488d+00
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1.d+01 PropType= Geometric nECoord= -6
    iElt= 16
    EltName= ExitPupil Element= Return

+----------------+-----------+-----------------------------------------+
| Surface=       | Conic     |                                         |
| KrElt= KcElt=  |           |                                         |
| psiElt=        | -1d22     |                                         |
| VptElt=        |           |                                         |
| RptElt=        | 0d0       |                                         |
| IndRef=        |           |                                         |
| Extinc= nObs=  | 0d0 0d0   |                                         |
| ApType= zElt=  | -1d0 0d0  |                                         |
| PropType=      | 0d0 0d0   |                                         |
| nECoord=       | 0d0 0d0   |                                         |
|                | 0d0       |                                         |
|                | 1         |                                         |
|                | .0488d+00 |                                         |
|                |           |                                         |
|                | 1d22      |                                         |
|                |           |                                         |
|                | 0         |                                         |
|                |           |                                         |
|                | None      |                                         |
|                | 1.d+01    |                                         |
|                |           |                                         |
|                | FarField  |                                         |
|                |           |                                         |
|                | -6        |                                         |
+================+===========+=========================================+
| iElt=          | 17        |                                         |
+----------------+-----------+-----------------------------------------+
| EltName=       | FP        |                                         |
| Element=       |           |                                         |
| Surface=       | F         |                                         |
| KrElt= KcElt=  | ocalPlane |                                         |
| psiElt=        | Flat      |                                         |
| VptElt=        |           |                                         |
| RptElt=        | -1d22     |                                         |
| IndRef=        |           |                                         |
| Extinc= nObs=  | 0d0       |                                         |
| ApType= zElt=  |           |                                         |
| PropType=      | 0d0 0d0   |                                         |
| nECoord=       | -1d0 0d0  |                                         |
|                | 0d0 1d3   |                                         |
|                | 0d0 0d0   |                                         |
|                | 0d0 1d0   |                                         |
|                |           |                                         |
|                | 1d22      |                                         |
|                |           |                                         |
|                | 0         |                                         |
|                |           |                                         |
|                | None 0d0  |                                         |
|                |           |                                         |
|                | Geometric |                                         |
|                |           |                                         |
|                | -6        |                                         |
+----------------+-----------+-----------------------------------------+
| nOutCord=      | 5         | 0d0 0d0 0d0 0d0                         |
|                |           |                                         |
| Tout=          | 1d0 0d0   |                                         |
|                | 0d0       |                                         |
+----------------+-----------+-----------------------------------------+
|                | 0d0 1d0   | 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 1d0 |
|                | 0d0 0d0   | 0d0 0d0                                 |
|                | 0d0 0d0   |                                         |
|                | 0d0 0d0   | 0d0 0d0 0d0 1d0                         |
|                | 0d0       |                                         |
|                |           |                                         |
|                | 0d0 0d0   |                                         |
|                | 0d0       |                                         |
+----------------+-----------+-----------------------------------------+

### Adaptive Optics Example

#### AOExample.jou

    % MACOS Adaptive Optics Example Macro
    % First show system response without atmosphere:
    old AOexample; % Load prescription
    draw 0 26 xz % Draw full system
    draw 6 26 xy % Draw AO optics only
    spo 1 b y % Pupil spot diagram
    opd 16 y % OPD in exit pupil shows deformable
    % mirror test pattern
    spo 26 Tout y % Hartman sensor spots on detector
    % Next show system response with atmosphere but no
    % deformable mirror:
    modify Surface(7)=Flat % Flatten deformable mirror q
    atmos 1 20 % Add atmospheric phase disturbance
    2.854331e-5 0 atm 12.5 % at entrance pupil
    opd 16 y % OPD in exit pupil
    spo 26 Tout y % Hartman sensor spots on detector

#### AOExample.in

    % MACOS Adaptive Optics Example Prescription

+------------------------+---------------------------------------------+
| ChfRayDir=             | 0d0 0d0 1d0                                 |
+========================+=============================================+
| ChfRayPos=             | 0d0 0d0 -2d+02                              |
+------------------------+---------------------------------------------+
| zSource=               | 1d22                                        |
+------------------------+---------------------------------------------+
| IndRef=                | 1d0                                         |
+------------------------+---------------------------------------------+
| Extinc=                | 0d0                                         |
+------------------------+---------------------------------------------+
| Wavelen=               | 2.854330709D-05                             |
+------------------------+---------------------------------------------+
| Flux=                  | 1d0                                         |
+------------------------+---------------------------------------------+
| GridType=              | Circular                                    |
+------------------------+---------------------------------------------+
| Aperture=              | 1.991873225632d02                           |
+------------------------+---------------------------------------------+
| Obscratn=              | 40.5d0                                      |
+------------------------+---------------------------------------------+
| nGridpts=              | 128                                         |
+------------------------+---------------------------------------------+
| xGrid=                 | 0d0 -1d0 0d0                                |
+------------------------+---------------------------------------------+
| yGrid=                 | -1d0 0d0 0d0                                |
+------------------------+---------------------------------------------+
| nElt=                  | 26                                          |
+------------------------+---------------------------------------------+
| iElt=                  | 1                                           |
+------------------------+---------------------------------------------+
| EltName=               | Pupil                                       |
+------------------------+---------------------------------------------+
| Element=               | Reference                                   |
+------------------------+---------------------------------------------+
| Surface=               | Conic                                       |
+------------------------+---------------------------------------------+
| KrElt=                 | -1d22                                       |
+------------------------+---------------------------------------------+
| KcElt=                 | 0d0                                         |
+------------------------+---------------------------------------------+
| psiElt=                | 0d0 0d0 -1d0                                |
+------------------------+---------------------------------------------+
| VptElt=                | 0d0 0d0 -6d0                                |
+------------------------+---------------------------------------------+
| RptElt=                | 0d0 0d0 0d0                                 |
+------------------------+---------------------------------------------+
| IndRef=                | 1d0                                         |
+------------------------+---------------------------------------------+
| Extinc=                | 0d0                                         |
+------------------------+---------------------------------------------+
| nObs=                  | 0                                           |
+------------------------+---------------------------------------------+
| ApType=                | None                                        |
+------------------------+---------------------------------------------+
| zElt=                  | 1d22                                        |
+------------------------+---------------------------------------------+
| PropType=              | Geometric                                   |
+------------------------+---------------------------------------------+
| nECoord=               | -6                                          |
+------------------------+---------------------------------------------+
| iElt=                  | 2                                           |
+------------------------+---------------------------------------------+
| EltName=               | PM                                          |
+------------------------+---------------------------------------------+
| Element=               | Reflector                                   |
+------------------------+---------------------------------------------+
| Surface=               | Conic                                       |
+------------------------+---------------------------------------------+
| KrElt=                 | -1.3357d+03                                 |
+------------------------+---------------------------------------------+
| KcElt=                 | -1d0                                        |
+------------------------+---------------------------------------------+
| psiElt=                | 0d0 0d0 -1d0                                |
+------------------------+---------------------------------------------+
| VptElt=                | 0d0 0d0 0d0                                 |
+------------------------+---------------------------------------------+
| RptElt=                | 0d0 0d0 0d0                                 |
+------------------------+---------------------------------------------+
| IndRef=                | 1d0                                         |
+------------------------+---------------------------------------------+
| Extinc=                | 0d0                                         |
+------------------------+---------------------------------------------+
| nObs=                  | 0                                           |
+------------------------+---------------------------------------------+
| ApType=                | None                                        |
+------------------------+---------------------------------------------+
| zElt=                  | 6.6785d+02                                  |
+------------------------+---------------------------------------------+
| PropType=              | Geometric                                   |
+------------------------+---------------------------------------------+
| nECoord=               | -6                                          |
+------------------------+---------------------------------------------+
| iElt=                  | 3                                           |
+------------------------+---------------------------------------------+
| EltName=               | SM                                          |
+------------------------+---------------------------------------------+

    Element= Reflector Surface= Conic
    KrElt= -3.184999999D+02 KcElt= -2.352756999D+00
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 -5.423064572D+02 RptElt= 0d0 0d0 -5.423064572D+02
    IndRef= 1d0
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1.256970591D+02
    PropType= Geometric nECoord= -6
    iElt= 4
    EltName= FM1 Element= Reflector Surface= Flat
    KrElt= -1d22 KcElt= 0d0
    psiElt= 7.071067812D-01 0d0 7.071067812D-01
    VptElt= 0d0 0d0 7.5375d+01 RptElt= 0d0 0d0 7.5375d+01 IndRef= 1d0
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 5
    EltName= OAP1 Element= Reflector Surface= Conic
    KrElt= -1.364d+02
    KcElt= -1d0
    psiElt= 9.668919598D-01 -2.551860851D-01 0d0 VptElt= -4.140665047D+01 1.740405528D+01 7.5375d+01 RptElt= -4.140665047D+01 1.740405528D+01 7.5375d+01
    IndRef= 1d0
    Extinc= 0d0
    nObs= 0
    ApType= None zElt= 6.82d+01
    PropType= Geometric nECoord= -6
    iElt= 6
    EltName= FSM Element= Reflector Surface= Flat
    KrElt= -1d22 KcElt= 0d0
    psiElt= 9.999929341D-01 3.759223637D-03 0d0 VptElt= 1.318856255D+01 -1.530759269D+01 7.5375d+01 RptElt= 1.318856255D+01 -1.530759269D+01 7.5375d+01
    IndRef= 1d0
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric
    nECoord= -6
    iElt= 7
    EltName= DM Element= Reflector Surface= UserDefined
    KrElt= -1d22

+-----------------------+--------------+-------------------------------+
| KcElt= 0d0            | -1.9         | \% 1 or 5                     |
|                       | 51098758D+01 |                               |
| UDSrfType= 1          |              | 7.5375d+01                    |
|                       |              |                               |
| UDSrfFile=            |              |                               |
| AOmirror.test         |              |                               |
|                       |              |                               |
| pMon=                 |              |                               |
| -2.249426003D+00      |              |                               |
+=======================+==============+===============================+
| xMon=                 | 8.8          | 0d0                           |
| -4.697040833D-01      | 28239202D-01 |                               |
+-----------------------+--------------+-------------------------------+
| yMon= 0d0 0d0 1d0     | 4.6          | 0d0                           |
|                       | 97040833D-01 |                               |
| zMon= 8.828239202D-01 |              |                               |
+-----------------------+--------------+-------------------------------+
| lMon= 0d0             | 4.6          | 0d0                           |
|                       | 97040833D-01 |                               |
| psiElt=               |              |                               |
| 8.828239202D-01       |              |                               |
+-----------------------+--------------+-------------------------------+
| VptElt=               | -1.9         | 7.5375d+01                    |
| -2.250697813D+00      | 50691308D+01 |                               |
+-----------------------+--------------+-------------------------------+
| RptElt=               | -1.9         | 7.5375d+01                    |
| -2.250697813D+00      | 50691308D+01 |                               |
|                       |              |                               |
| IndRef= 1d0           |              |                               |
|                       |              |                               |
| Extinc= 0d0           |              |                               |
|                       |              |                               |
| nObs= 0               |              |                               |
|                       |              |                               |
| ApType= None          |              |                               |
|                       |              |                               |
| zElt= 1d22 PropType=  |              |                               |
| Geometric nECoord= -6 |              |                               |
+-----------------------+--------------+-------------------------------+
| iElt= 8               | 9.9          | 0d0                           |
|                       | 59082335D-01 |                               |
| EltName= FM2 Element= |              |                               |
| Reflector Surface=    |              |                               |
| Flat                  |              |                               |
|                       |              |                               |
| KrElt= -1d22 KcElt=   |              |                               |
| 0d0                   |              |                               |
|                       |              |                               |
| psiElt=               |              |                               |
| 9.037029675D-02       |              |                               |
+-----------------------+--------------+-------------------------------+
| VptElt=               | 2.5          | 7.5375d+01                    |
| 2.331918268D+01       | 76404533D+00 |                               |
+-----------------------+--------------+-------------------------------+
| RptElt=               | 2.5          | 7.5375d+01                    |
| 2.331918268D+01       | 76404533D+00 |                               |
|                       |              |                               |
| IndRef= 1d0           |              |                               |
|                       |              |                               |
| Extinc= 0d0           |              |                               |
|                       |              |                               |
| nObs= 0               |              |                               |
|                       |              |                               |
| ApType= None          |              |                               |
|                       |              |                               |
| zElt= 1d22 PropType=  |              |                               |
| Geometric nECoord= -6 |              |                               |
+-----------------------+--------------+-------------------------------+
| iElt= 9               | 7.7          | 0d0                           |
|                       | 91729129D-01 |                               |
| EltName= OAP2         |              |                               |
| Element= Reflector    |              |                               |
| Surface= Conic        |              |                               |
|                       |              |                               |
| KrElt= -1.364d+02     |              |                               |
|                       |              |                               |
| KcElt= -1d0           |              |                               |
|                       |              |                               |
| psiElt=               |              |                               |
| -6.268090394D-01      |              |                               |
+-----------------------+--------------+-------------------------------+
| VptElt=               | -2.8         | 7.5375d+01                    |
| 2.538458250D+01       | 22361844D+01 |                               |
+-----------------------+--------------+-------------------------------+
| RptElt=               | -2.8         | 7.5375d+01                    |
| 2.538458250D+01       | 22361844D+01 |                               |
|                       |              |                               |
| IndRef= 1d0           |              |                               |
|                       |              |                               |
| Extinc= 0d0           |              |                               |
|                       |              |                               |
| nObs= 0               |              |                               |
|                       |              |                               |
| ApType= None zElt=    |              |                               |
| 6.82d+01              |              |                               |
|                       |              |                               |
| PropType= Geometric   |              |                               |
|                       |              |                               |
| nECoord= -6           |              |                               |
+-----------------------+--------------+-------------------------------+

+----------------+------------+-----------+---------------------------+
| iElt= EltName= | 10         | 1.1146    | 0d0                       |
| Element=       |            | 89322D-01 |                           |
| Surface=       | SSM1       |           |                           |
| KrElt= KcElt=  |            |           |                           |
|                | Reflector  |           |                           |
| psiElt=        | Flat       |           |                           |
|                |            |           |                           |
|                | -1d22      |           |                           |
|                |            |           |                           |
|                | 0d0        |           |                           |
|                |            |           |                           |
|                | -9.937     |           |                           |
|                | 679192D-01 |           |                           |
+================+============+===========+===========================+
| VptElt=        | 5.421      | 8.1171    | 7.5375d+01                |
|                | 565944D+00 | 54092D+00 |                           |
+----------------+------------+-----------+---------------------------+
| RptElt=        | 5.421      | 8.1171    | 7.5375d+01                |
| IndRef=        | 565944D+00 | 54092D+00 |                           |
| Extinc= nObs=  |            |           |                           |
| ApType= zElt=  | 1d0 0d0 0  |           |                           |
| PropType=      |            |           |                           |
| nECoord=       | None 1d22  |           |                           |
|                |            |           |                           |
|                | Geometric  |           |                           |
|                |            |           |                           |
|                | -6         |           |                           |
+----------------+------------+-----------+---------------------------+
| iElt= EltName= | 11         | 7.6719    | 0d0                       |
| Element=       |            | 02813D-02 |                           |
| Surface=       | SSM2       |           |                           |
| KrElt= KcElt=  |            |           |                           |
| psiElt=        | Reflector  |           |                           |
|                | Flat       |           |                           |
|                |            |           |                           |
|                | -1d22      |           |                           |
|                |            |           |                           |
|                | 0d0        |           |                           |
|                |            |           |                           |
|                | -9.970     |           |                           |
|                | 527522D-01 |           |                           |
+----------------+------------+-----------+---------------------------+
| VptElt=        | 1.733      | 1.3321    | 7.5375d+01                |
|                | 428143D+01 | 69152D+01 |                           |
+----------------+------------+-----------+---------------------------+
| RptElt=        | 1.733      | 1.3321    | 7.5375d+01                |
| IndRef=        | 428143D+01 | 69152D+01 |                           |
| Extinc= nObs=  |            |           |                           |
| ApType= zElt=  | 1d0 0d0 0  |           |                           |
| PropType=      |            |           |                           |
| nECoord=       | None 1d22  |           |                           |
|                |            |           |                           |
|                | Geometric  |           |                           |
|                |            |           |                           |
|                | -6         |           |                           |
+----------------+------------+-----------+---------------------------+
| iElt=          | 12         |           |                           |
+----------------+------------+-----------+---------------------------+
| EltName=       | FS         | 2.8234    | 0d0                       |
| Element=       |            | 14568D-01 |                           |
| Surface=       | Reflector  |           |                           |
| KrElt= KcElt=  | Flat       |           |                           |
|                |            |           |                           |
| psiElt=        | -1d22      |           |                           |
|                |            |           |                           |
|                | 0d0        |           |                           |
|                |            |           |                           |
|                | -9.593     |           |                           |
|                | 139745D-01 |           |                           |
+----------------+------------+-----------+---------------------------+
| VptElt=        | 4.387      | 2.1538    | 7.5375d+01                |
|                | 255193D+00 | 12421D+01 |                           |
+----------------+------------+-----------+---------------------------+
| RptElt=        | 4.387      | 2.1538    | 7.5375d+01                |
| IndRef=        | 255193D+00 | 12421D+01 |                           |
| Extinc= nObs=  |            |           |                           |
| ApType= zElt=  | 1d0 0d0 0  |           |                           |
| PropType=      |            |           |                           |
| nECoord=       | None 1d22  |           |                           |
|                |            |           |                           |
|                | Geometric  |           |                           |
|                |            |           |                           |
|                | -6         |           |                           |
+----------------+------------+-----------+---------------------------+
| iElt= EltName= | 13         | -6.9812   | 0d0                       |
| Element=       |            | 60298D-03 |                           |
| Surface=       | WFSRL1A    |           |                           |
| KrElt= KcElt=  |            |           |                           |
| psiElt=        | Refractor  |           |                           |
|                | Conic      |           |                           |
|                |            |           |                           |
|                | -1.550     |           |                           |
|                | 872636D+00 |           |                           |
|                |            |           |                           |
|                | 0d0        |           |                           |
|                |            |           |                           |
|                | 9.999      |           |                           |
|                | 756307D-01 |           |                           |
+----------------+------------+-----------+---------------------------+
| VptElt=        | 8.496      | 2.1509    | 7.5375d+01                |
|                | 015971D+00 | 43918D+01 |                           |
+----------------+------------+-----------+---------------------------+
| RptElt=        | 8.496      | 2.1509    | 7.5375d+01                |
| IndRef=        | 015971D+00 | 43918D+01 |                           |
| Extinc= nObs=  |            |           |                           |
|                | 1.512      |           |                           |
|                | 425970D+00 |           |                           |
|                |            |           |                           |
|                | 0d0 0      |           |                           |
+----------------+------------+-----------+---------------------------+

+----------+--------------------+------------------+------------------+
| ApType=  | None               |                  |                  |
| zElt=    | 5.288796394D+00    |                  |                  |
| P        |                    |                  |                  |
| ropType= | Geometric          |                  |                  |
| nECoord= |                    |                  |                  |
|          | -6                 |                  |                  |
+==========+====================+==================+==================+
| iElt=    | 14                 | 6.981260298D-03  | 0d0              |
| EltName= |                    |                  |                  |
| Element= | WFSRL1B            |                  |                  |
| Surface= |                    |                  |                  |
| KrElt=   | Refractor Conic    |                  |                  |
| KcElt=   |                    |                  |                  |
| psiElt=  | -4.419109771D+00   |                  |                  |
|          |                    |                  |                  |
|          | 0d0                |                  |                  |
|          |                    |                  |                  |
|          | -9.999756307D-01   |                  |                  |
+----------+--------------------+------------------+------------------+
| VptElt=  | 8.996003786D+00    | 2.150594855D+01  | 7.5375d+01       |
+----------+--------------------+------------------+------------------+
| RptElt=  | 8.996003786D+00    | 2.150594855D+01  | 7.5375d+01       |
| IndRef=  |                    |                  |                  |
| Extinc=  | 1.789618121D+00    |                  |                  |
| nObs=    |                    |                  |                  |
| ApType=  | 0d0 0              |                  |                  |
| zElt=    |                    |                  |                  |
| P        | None               |                  |                  |
| ropType= | 3.411018511D+00    |                  |                  |
| nECoord= |                    |                  |                  |
|          | Geometric          |                  |                  |
|          |                    |                  |                  |
|          | -6                 |                  |                  |
+----------+--------------------+------------------+------------------+
| iElt=    | 15                 | -6.981260298D-03 | 0d0              |
| EltName= |                    |                  |                  |
| Element= | WFSRL1C            |                  |                  |
| Surface= |                    |                  |                  |
| KrElt=   | Refractor Conic    |                  |                  |
| KcElt=   |                    |                  |                  |
| psiElt=  | -2.771554824D+01   |                  |                  |
|          |                    |                  |                  |
|          | 0d0                |                  |                  |
|          |                    |                  |                  |
|          | 9.999756307D-01    |                  |                  |
+----------+--------------------+------------------+------------------+
| VptElt=  | 9.495991602D+00    | 2.150245792D+01  | 7.5375d+01       |
+----------+--------------------+------------------+------------------+
| RptElt=  | 9.495991602D+00    | 2.150245792D+01  | 7.5375d+01       |
| IndRef=  |                    |                  |                  |
| Extinc=  | 1d0 0d0 0          |                  |                  |
| nObs=    |                    |                  |                  |
| ApType=  | None               |                  |                  |
| zElt=    | 2.198859945D+01    |                  |                  |
| P        |                    |                  |                  |
| ropType= | Geometric          |                  |                  |
| nECoord= |                    |                  |                  |
|          | -6                 |                  |                  |
+----------+--------------------+------------------+------------------+
| iElt=    | 16                 | -6.981321726D-03 | 3.414481680D-15  |
| EltName= |                    |                  |                  |
| Element= | LensArrayReference |                  |                  |
| Surface= | Reference          |                  |                  |
| KrElt=   |                    |                  |                  |
| KcElt=   | Conic              |                  |                  |
| psiElt=  |                    |                  |                  |
|          | -1d22              |                  |                  |
|          |                    |                  |                  |
|          | 0d0                |                  |                  |
|          |                    |                  |                  |
|          | 9.999756303D-01    |                  |                  |
+----------+--------------------+------------------+------------------+
| VptElt=  | 1.272900024D+01    | 2.147992068D+01  | 7.537500000D+01  |
+----------+--------------------+------------------+------------------+
| RptElt=  | 1.272900024D+01    | 2.147992068D+01  | 7.537500000D+01  |
| IndRef=  |                    |                  |                  |
| Extinc=  | 1d0 0d0 0          |                  |                  |
| nObs=    |                    |                  |                  |
| ApType=  | None 1d22          |                  |                  |
| zElt=    |                    |                  |                  |
| P        | Geometric          |                  |                  |
| ropType= |                    |                  |                  |
| nECoord= | -6                 |                  |                  |
+----------+--------------------+------------------+------------------+
| iElt=    | 17                 |                  |                  |
| EltName= |                    |                  |                  |
| Element= | LensArray1A        |                  |                  |
|          |                    |                  |                  |
|          | LensArray          |                  |                  |
+----------+--------------------+------------------+------------------+

    Surface= Conic
    KrElt= -0.47244094488189d+00
    KcElt= 0d0
    psiElt= 9.999756303D-01 -6.981321725D-03 0d0
    VptElt= 12.7291897290d0 21.4720113661d0 75.3671259843d0 RptElt= 1.272924470D+01 2.147988519D+01 7.5375d+01
    LensArrayType= 2
    LensArrayWidth= 0.01574803149606
    pMon= 1.272924470D+01 2.147988519D+01 7.5375d+01
    xMon= 0d0 0d0 -1d0
    yMon= -6.981321725D-03 -9.999756303D-01 0d0 zMon= -9.999756303D-01 6.981321725D-03 0d0
    IndRef= 1.512425970D+00
    Extinc= 0d0
    nObs= 0
    ApType= None zElt= 5.17d+00
    PropType= Geometric nECoord= -6
    iElt= 18
    EltName= LensArray1B Element= Refractor Surface= Flat
    KrElt= -1d22 KcElt= 0d0
    psiElt= -9.999756307D-01 6.981260298D-03 0d0 VptElt= 1.297923860D+01 2.147813988D+01 7.5375d+01 RptElt= 1.297923860D+01 2.147813988D+01 7.5375d+01
    IndRef= 1d0
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 19
    EltName= FieldLensA Element= Refractor Surface= Conic
    KrElt= -2.122047244D+00
    KcElt= 0d0
    psiElt= 9.999756307D-01 -6.981260298D-03 0d0 VptElt= 1.358583579D+01 2.147390496D+01 7.5375d+01 RptElt= 1.358583579D+01 2.147390496D+01 7.5375d+01 IndRef= 1.677871170D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1.359427521D+00
    PropType= Geometric nECoord= -6
    iElt= 20
    EltName= FieldLensB Element= Refractor Surface= Conic
    KrElt= -3.224409449D+00
    KcElt= 0d0
    psiElt= 9.999756307D-01 -6.981260298D-03 0d0 VptElt= 1.374331196D+01 2.147280555D+01 7.5375d+01 RptElt= 1.374331196D+01 2.147280555D+01 7.5375d+01
    IndRef= 1d0

+-----------+---------------------+---------------------+-------------+
| Extinc=   | 0d0                 |                     |             |
+===========+=====================+=====================+=============+
| nObs=     | 0                   |                     |             |
| ApType=   |                     |                     |             |
| zElt=     | None                |                     |             |
| PropType= | 1.359427521D+00     |                     |             |
| nECoord=  |                     |                     |             |
|           | Geometric           |                     |             |
|           |                     |                     |             |
|           | -6                  |                     |             |
+-----------+---------------------+---------------------+-------------+
| iElt=     | 21                  | -6.981260298D-03    | 0d0         |
| EltName=  |                     |                     |             |
| Element=  | WFSLAR1A            |                     |             |
| Surface=  |                     |                     |             |
| KrElt=    | Refractor Conic     |                     |             |
| KcElt=    |                     |                     |             |
| psiElt=   | -8.055118110D-01    |                     |             |
|           |                     |                     |             |
|           | 0d0                 |                     |             |
|           |                     |                     |             |
|           | 9.999756307D-01     |                     |             |
+-----------+---------------------+---------------------+-------------+
| VptElt=   | 1.880962705D+01     | 2.143743542D+01     | 7.5375d+01  |
+-----------+---------------------+---------------------+-------------+
| RptElt=   | 1.880962705D+01     | 2.143743542D+01     | 7.5375d+01  |
| IndRef=   |                     |                     |             |
| Extinc=   | 1.611736488D+00     |                     |             |
| nObs=     |                     |                     |             |
| ApType=   | 0d0 0               |                     |             |
| zElt=     |                     |                     |             |
| PropType= | None                |                     |             |
| nECoord=  | 1.359427521D+00     |                     |             |
|           |                     |                     |             |
|           | Geometric           |                     |             |
|           |                     |                     |             |
|           | -6                  |                     |             |
+-----------+---------------------+---------------------+-------------+
| iElt=     | 22                  | 6.981260298D-03     | 0d0         |
| EltName=  |                     |                     |             |
| Element=  | WFSLAR1B            |                     |             |
| Surface=  |                     |                     |             |
| KrElt=    | Refractor Conic     |                     |             |
| KcElt=    |                     |                     |             |
| psiElt=   | -5.090551181D-01    |                     |             |
|           |                     |                     |             |
|           | 0d0                 |                     |             |
|           |                     |                     |             |
|           | -9.999756307D-01    |                     |             |
+-----------+---------------------+---------------------+-------------+
| VptElt=   | 1.896631614D+01     | 2.143634151D+01     | 7.5375d+01  |
+-----------+---------------------+---------------------+-------------+
| RptElt=   | 1.896631614D+01     | 2.143634151D+01     | 7.5375d+01  |
| IndRef=   |                     |                     |             |
| Extinc=   | 1.677871170D+00     |                     |             |
| nObs=     |                     |                     |             |
| ApType=   | 0d0 0               |                     |             |
| zElt=     |                     |                     |             |
| PropType= | None                |                     |             |
| nECoord=  | 6.038543863D-01     |                     |             |
|           |                     |                     |             |
|           | Geometric           |                     |             |
|           |                     |                     |             |
|           | -6                  |                     |             |
+-----------+---------------------+---------------------+-------------+
| iElt=     | 23                  | 6.981260298D-03     | 0d0         |
| EltName=  |                     |                     |             |
| Element=  | WFSLAR1C            |                     |             |
| Surface=  |                     |                     |             |
| KrElt=    | Refractor Conic     |                     |             |
| KcElt=    |                     |                     |             |
| psiElt=   | -2.980708661D+00    |                     |             |
|           |                     |                     |             |
|           | 0d0                 |                     |             |
|           |                     |                     |             |
|           | -9.999756307D-01    |                     |             |
+-----------+---------------------+---------------------+-------------+
| VptElt=   | 1.901946445D+01     | 2.143597046D+01     | 7.5375d+01  |
+-----------+---------------------+---------------------+-------------+
| RptElt=   | 1.901946445D+01     | 2.143597046D+01     | 7.5375d+01  |
| IndRef=   |                     |                     |             |
| Extinc=   | 1d0 0d0 0           |                     |             |
| nObs=     |                     |                     |             |
| ApType=   | None                |                     |             |
| zElt=     | 2.027987141D+00     |                     |             |
| PropType= |                     |                     |             |
| nECoord=  | Geometric           |                     |             |
|           |                     |                     |             |
|           | -6                  |                     |             |
+-----------+---------------------+---------------------+-------------+
| iElt=     | 24                  |                     |             |
+-----------+---------------------+---------------------+-------------+

+----------------+------------+-----------+---------------------------+
| EltName=       | CCDWinA    |           |                           |
| Element=       |            |           |                           |
| Surface=       | Refractor  |           |                           |
|                |            |           |                           |
|                | Conic      |           |                           |
+================+============+===========+===========================+
| KrElt= KcElt=  | -1d22      | 6.9812    | 0d0                       |
|                |            | 60298D-03 |                           |
| psiElt=        | 0d0        |           |                           |
|                |            |           |                           |
|                | -9.999     |           |                           |
|                | 756307D-01 |           |                           |
+----------------+------------+-----------+---------------------------+
| VptElt=        | 2.027      | 2.1427    | 7.5375d+01                |
|                | 861775D+01 | 17977D+01 |                           |
+----------------+------------+-----------+---------------------------+
| RptElt=        | 2.027      | 2.1427    | 7.5375d+01                |
| IndRef=        | 861775D+01 | 17977D+01 |                           |
| Extinc= nObs=  |            |           |                           |
| ApType= zElt=  | 1.512      |           |                           |
| PropType=      | 425970D+00 |           |                           |
| nECoord=       |            |           |                           |
|                | 0d0 0      |           |                           |
|                |            |           |                           |
|                | None 1d22  |           |                           |
|                |            |           |                           |
|                | Geometric  |           |                           |
|                |            |           |                           |
|                | -6         |           |                           |
+----------------+------------+-----------+---------------------------+
| iElt= EltName= | 25         | 6.9812    | 0d0                       |
| Element=       |            | 60298D-03 |                           |
| Surface=       | CCDWinB    |           |                           |
| KrElt= KcElt=  |            |           |                           |
| psiElt=        | Refractor  |           |                           |
|                | Conic      |           |                           |
|                |            |           |                           |
|                | -1d22      |           |                           |
|                |            |           |                           |
|                | 0d0        |           |                           |
|                |            |           |                           |
|                | -9.999     |           |                           |
|                | 756307D-01 |           |                           |
+----------------+------------+-----------+---------------------------+
| VptElt=        | 2.040      | 2.1426    | 7.5375d+01                |
|                | 361470D+01 | 30711D+01 |                           |
+----------------+------------+-----------+---------------------------+
| RptElt=        | 2.040      | 2.1426    | 7.5375d+01                |
| IndRef=        | 361470D+01 | 30711D+01 |                           |
| Extinc= nObs=  |            |           |                           |
| ApType= zElt=  | 1d0 0d0 0  |           |                           |
| PropType=      |            |           |                           |
| nECoord=       | None 1d22  |           |                           |
|                |            |           |                           |
|                | Geometric  |           |                           |
|                |            |           |                           |
|                | -6         |           |                           |
+----------------+------------+-----------+---------------------------+
| iElt=          | 26         |           |                           |
+----------------+------------+-----------+---------------------------+
| EltName=       | CCD        | 6.8571    | 0d0                       |
| Element=       |            | 53362D-03 |                           |
| Surface=       | FocalPlane |           |                           |
| KrElt= KcElt=  | Flat       |           |                           |
|                |            |           |                           |
| psiElt=        | -1.000     |           |                           |
|                | 000000D+22 |           |                           |
|                | 0.000      |           |                           |
|                | 000000D+00 |           |                           |
|                |            |           |                           |
|                | -9.999     |           |                           |
|                | 764894D-01 |           |                           |
+----------------+------------+-----------+---------------------------+
| VptElt=        | 2.049      | 2.1425    | 7.537500000D+01           |
|                | 655742D+01 | 65820D+01 |                           |
+----------------+------------+-----------+---------------------------+
| RptElt=        | 2.049      | 2.1425    | 7.537500000D+01           |
| IndRef=        | 655742D+01 | 65820D+01 |                           |
| Extinc= nObs=  |            |           |                           |
| ApType= zElt=  | 1.000      |           |                           |
| PropType=      | 000000D+00 |           |                           |
| nECoord=       |            |           |                           |
|                | 0.000      |           |                           |
|                | 000000D+00 |           |                           |
|                |            |           |                           |
|                | 0          |           |                           |
|                |            |           |                           |
|                | None       |           |                           |
|                | 1.000      |           |                           |
|                | 000000D+22 |           |                           |
|                |            |           |                           |
|                | Geometric  |           |                           |
|                |            |           |                           |
|                | -6         |           |                           |
+----------------+------------+-----------+---------------------------+
| nOutCord=      | 5          | -9.9997   | 0d0 0d0 0d0 0d0 0d0       |
|                |            | 56307D-01 |                           |
| Tout=          | -6.981     |           |                           |
|                | 260298D-03 |           |                           |
+----------------+------------+-----------+---------------------------+

    0d0 0d0 -1d0 0d0 0d0 0d0 0d0 0d0
    0d0 0d0 0d0 -6.981260298D-03 -9.999756307D-01 0d0 0d0
    0d0 0d0 0d0 0d0 0d0 -1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0

#### AOmirror.test

    5.5000000e+01 4e-05
    7.6000000e+01 -4e-05
    7.7000000e+01 -4e-05
    9.5000000e+01 -4e-05
    9.6000000e+01 -4e-05

+--------------------------------------------+-------------------------+
| 1.0200000e+02                              | 3e-05                   |
+============================================+=========================+
| 1.0300000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 1.0400000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 1.0500000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 1.0600000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 1.0700000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 1.2300000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 1.2400000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 1.2500000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 1.2600000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 1.2700000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 1.2800000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 2.0700000e+02                              | 5e-06                   |
+--------------------------------------------+-------------------------+
| 2.2800000e+02                              | 5e-06                   |
+--------------------------------------------+-------------------------+
| 2.0800000e+02                              | 1e-05                   |
+--------------------------------------------+-------------------------+
| 2.2900000e+02                              | 1e-05                   |
+--------------------------------------------+-------------------------+
| 2.4900000e+02                              | 1e-05                   |
+--------------------------------------------+-------------------------+
| 2.6800000e+02                              | 1e-05                   |
+--------------------------------------------+-------------------------+
| 2.0900000e+02                              | 1.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.3000000e+02                              | 1.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.5000000e+02                              | 1.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.6900000e+02                              | 1.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.8700000e+02                              | 1.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.1000000e+02                              | 2e-05                   |
+--------------------------------------------+-------------------------+
| 2.3100000e+02                              | 2e-05                   |
+--------------------------------------------+-------------------------+
| 2.5100000e+02                              | 2e-05                   |
+--------------------------------------------+-------------------------+
| 2.7000000e+02                              | 2e-05                   |
+--------------------------------------------+-------------------------+
| 2.8800000e+02                              | 2e-05                   |
+--------------------------------------------+-------------------------+
| 3.0400000e+02                              | 2e-05                   |
+--------------------------------------------+-------------------------+
| 2.1100000e+02                              | 2.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.3200000e+02                              | 2.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.5200000e+02                              | 2.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.7100000e+02                              | 2.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.8900000e+02                              | 2.5e-05                 |
+--------------------------------------------+-------------------------+
| 3.0500000e+02                              | 2.5e-05                 |
+--------------------------------------------+-------------------------+
| 3.1900000e+02                              | 2.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.1200000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 2.3300000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 2.5300000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 2.7200000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 2.9000000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 3.0600000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 3.2000000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 3.3200000e+02                              | 3e-05                   |
+--------------------------------------------+-------------------------+
| 2.1300000e+02                              | 3.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.3400000e+02                              | 3.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.5400000e+02                              | 3.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.7300000e+02                              | 3.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.9100000e+02                              | 3.5e-05                 |
+--------------------------------------------+-------------------------+
| 3.0700000e+02                              | 3.5e-05                 |
+--------------------------------------------+-------------------------+
| 3.2100000e+02                              | 3.5e-05                 |
+--------------------------------------------+-------------------------+
| 3.3300000e+02                              | 3.5e-05                 |
+--------------------------------------------+-------------------------+
| 2.1400000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 2.3500000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 2.5500000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 2.7400000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 2.9200000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 3.0800000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+
| 3.2200000e+02                              | 4e-05                   |
+--------------------------------------------+-------------------------+

+----------------------+-----------------------------------------------+
| 3.3400000e+02        | 4e-05                                         |
+======================+===============================================+
| 3.4300000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.4200000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.4300000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.4400000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.4500000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.2000000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.2100000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.2200000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.2300000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.2400000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 1.9900000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.0000000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.0100000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.0200000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 2.2000000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 1.7800000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 1.7900000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 1.8000000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 1.8100000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 1.8200000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+
| 1.7500000e+02        | 4e-05                                         |
+----------------------+-----------------------------------------------+

### Segmented Telescope Example

#### SegDemo.jou

    % SegmentDemo.jou -- a MACOS demonstration macro that illustrates the gener-ation
    % of a linear model of a segmented mirror system. The linear model is used to
    % generate an image and OPD numbers. These are compared with the same num-bers
    % computed using the full nonlinear model. They agree very closely.
    % Run this macro with macos256
    old SegmentDemo; % Load prescription file
    %
    draw 0 10 xz; % Sketch optics layout
    %
    map 0 % Displays the ray mapping in the entrance pupil; also
    % shows the segment mapping
    %
    intensity 10 % Generate simulated image at the focal plane -- no
    % aberrations
    %
    log 10 % Display simulated unaberrated image with log10 stretch
    %
    %
    perturb 3 % Perturb the 3rd element using the full nonlinear model NO % Use global coordinates
    0,1e-4,0 % Rotational perturbation: rotate 100 urad about y 0,0,0 % No translational perturbation
    %
    opd 9 % Display wavefront error using linear model
    %
    log 10 % Display image (log10 stretch) using full nonlinear
    % model
    %
    old SegmentDemo; % Reload prescription file
    %
    build 9 % Compute linear optical model
    %
    lperturb 3 % Perturb the 3rd element using the linear model NO % Use global coordinates
    0,1e-4,0 % Rotational perturbation: rotate 100 urad about y 0,0,0 % No translational perturbation
    %
    lopd % Display wavefront error using linear model
    %
    llog % Display image (log10 stretch) using linear model

#### SegDemo.in

    % SegmentDemo.in -- prescription file for SegmentDemo.jou, a MACOS demon-stration
    % macro.

+---------------------+--------------------------------------+--------+
| ChfRayDir=          | 0d0 0d0 -1d0                         |        |
+=====================+======================================+========+
| ChfRayPos=          | 0d0 0d0 1.5d+00                      |        |
+---------------------+--------------------------------------+--------+
| zSource=            | 1.d+19                               |        |
+---------------------+--------------------------------------+--------+
| IndRef=             | 1d0                                  |        |
+---------------------+--------------------------------------+--------+
| Extinc=             | 0d0                                  |        |
+---------------------+--------------------------------------+--------+
| Wavelen=            | 1.d-04                               |        |
+---------------------+--------------------------------------+--------+
| Flux=               | 1d0                                  |        |
+---------------------+--------------------------------------+--------+
| GridType=           | Pie                                  |        |
+---------------------+--------------------------------------+--------+
| Aperture=           | 3.65d+00                             |        |
+---------------------+--------------------------------------+--------+
| Obscratn=           | 3.d-01                               |        |
+---------------------+--------------------------------------+--------+
| nGridpts=           | 25                                   |        |
+---------------------+--------------------------------------+--------+
| xGrid=              | 1d0 0d0 0d0                          |        |
+---------------------+--------------------------------------+--------+
| yGrid=              | 0d0 1d0 0d0                          |        |
+---------------------+--------------------------------------+--------+
| nElt=               | 10                                   |        |
+---------------------+--------------------------------------+--------+
| nSeg=               | 7                                    |        |
+---------------------+--------------------------------------+--------+
| width=              | 1.2d+00                              |        |
+---------------------+--------------------------------------+--------+
| gap=                | 0d0                                  |        |
+---------------------+--------------------------------------+--------+
| SegXgrid=           | 1d0 0d0 0d0                          |        |
+---------------------+--------------------------------------+--------+
| SegCoord=           | 0 0                                  | 0      |
+---------------------+--------------------------------------+--------+
|                     | 1 -1                                 | -2     |
+---------------------+--------------------------------------+--------+
|                     | 2 1                                  | -1     |
+---------------------+--------------------------------------+--------+
|                     | 1 2                                  | 1      |
+---------------------+--------------------------------------+--------+
|                     | -1 1                                 | 2      |
+---------------------+--------------------------------------+--------+
|                     | -2 -1                                | 1      |
+---------------------+--------------------------------------+--------+
|                     | -1 -2                                | -1     |
+---------------------+--------------------------------------+--------+

    iElt= 1
    EltName= S1 Element= Segment Surface= Conic
    KrElt= -2.92d+00
    KcElt= -1d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 0d0 RptElt= 0d0 0d0 0d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 1.46d+00
    PropType= Geometric
    nECoord= -6
    iElt= 2
    EltName= S2 Element= Segment Surface= Conic
    KrElt= -2.92d+00
    KcElt= -1d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 0d0
    RptElt= 6.124997330D-01 -1.060881274D+00 2.569563356D-01
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 1.46d+00
    PropType= Geometric nECoord= 6
    TElt= 4.610697370D-01 8.660255296D-01 -1.934282287D-01 0d0 0d0 0d0
    -7.985966745D-01 4.999997821D-01 3.350277145D-01 0d0 0d0 0d0 3.868566261D-01 -1.885624059D-18 9.221398760D-01 0d0 0d0 0d0
    0d0 0d0 0d0 4.610697370D-01 8.660255296D-01 -1.934282287D-01
    0d0 0d0 0d0 -7.985966745D-01 4.999997821D-01 3.350277145D-01
    0d0 0d0 0d0 3.868566261D-01 -1.885624059D-18 9.221398760D-01
    iElt= 3
    EltName= S3 Element= Segment Surface= Conic
    KrElt= -2.92d+00
    KcElt= -1d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 0d0
    RptElt= 1.225d+00 -3.699397078D-07 2.569563356D-01
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 1.46d+00
    PropType= Geometric nECoord= 6
    TElt= 9.221398760D-01 3.019915982D-07 -3.868566261D-01 0d0 0d0 0d0
    -2.784784949D-07 1d0 1.168274508D-07 0d0 0d0 0d0 3.868566261D-01 3.238991466D-24 9.221398760D-01 0d0 0d0 0d0
    0d0 0d0 0d0 9.221398760D-01 3.019915982D-07 -3.868566261D-01
    0d0 0d0 0d0 -2.784784949D-07 1d0 1.168274508D-07
    0d0 0d0 0d0 3.868566261D-01 3.238991466D-24 9.221398760D-01
    iElt= 4
    EltName= S4 Element= Segment Surface= Conic
    KrElt= -2.92d+00
    KcElt= -1d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 0d0
    RptElt= 6.125003738D-01 1.060880904D+00 2.569563356D-01
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 1.46d+00
    PropType= Geometric nECoord= 6

+-----+-------------------------------------------------+---+---+---+
| TE  | 4.610702194D-01 -8.660252276D-01                | 0 | 0 | 0 |
| lt= | -1.934284311D-01                                | d | d | d |
|     |                                                 | 0 | 0 | 0 |
+=====+=================================================+===+===+===+
|     | 7.985963960D-01 5.000003051D-01                 | 0 | 0 | 0 |
|     | -3.350275976D-01                                | d | d | d |
|     |                                                 | 0 | 0 | 0 |
+-----+-------------------------------------------------+---+---+---+
|     | 3.868566261D-01 2.270467800D-17 9.221398760D-01 | 0 | 0 | 0 |
|     |                                                 | d | d | d |
|     |                                                 | 0 | 0 | 0 |
+-----+-------------------------------------------------+---+---+---+

    0d0 0d0 0d0 4.610702194D-01 -8.660252276D-01 -1.934284311D-01
    0d0 0d0 0d0 7.985963960D-01 5.000003051D-01 -3.350275976D-01
    0d0 0d0 0d0 3.868566261D-01 2.270467800D-17 9.221398760D-01
    iElt= 5
    EltName= S5 Element= Segment Surface= Conic
    KrElt= -2.92d+00
    KcElt= -1d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 0d0
    RptElt= -6.124995728D-01 1.060881366D+00 2.569563356D-01
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 1.46d+00
    PropType= Geometric nECoord= 6
    TElt= -4.610696164D-01 -8.660256051D-01 1.934281781D-01 0d0 0d0 0d0 7.985967441D-01 -4.999996513D-01 -3.350277437D-01 0d0 0d0 0d0
    3.868566261D-01 -2.095137844D-19 9.221398760D-01 0d0 0d0 0d0
    0d0 0d0 0d0 -4.610696164D-01 -8.660256051D-01 1.934281781D-01
    0d0 0d0 0d0 7.985967441D-01 -4.999996513D-01 -3.350277437D-01
    0d0 0d0 0d0 3.868566261D-01 -2.095137844D-19 9.221398760D-01
    iElt= 6
    EltName= S6 Element= Segment Surface= Conic
    KrElt= -2.92d+00
    KcElt= -1d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 0d0
    RptElt= -1.225d+00 5.549095617D-07 2.569563356D-01
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 1.46d+00
    PropType= Geometric nECoord= 6
    TElt= -9.221398760D-01 -4.529873973D-07 3.868566261D-01 0d0 0d0 0d0 4.177177424D-07 -1d0 -1.752411762D-07 0d0 0d0 0d0
    3.868566261D-01 -1.071601397D-23 9.221398760D-01 0d0 0d0 0d0
    0d0 0d0 0d0 -9.221398760D-01 -4.529873973D-07 3.868566261D-01
    0d0 0d0 0d0 4.177177424D-07 -1d0 -1.752411762D-07
    0d0 0d0 0d0 3.868566261D-01 -1.071601397D-23 9.221398760D-01
    iElt= 7
    EltName= S7 Element= Segment Surface= Conic
    KrElt= -2.92d+00
    KcElt= -1d0

+----------------+-------------------------+----------+---------------+
| psiElt=        | 0d0 0d0 1d0 0d0 0d0 0d0 | 2.56956  |               |
| VptElt=        |                         | 3356D-01 |               |
| RptElt=        | -6.125005340D-01        |          |               |
| IndRef=        | -1.060880811D+00        | 1.93428  |               |
| Extinc= nObs=  |                         | 4817D-01 |               |
| ApType= zElt=  | 1d0 1d22                | 0d0 0d0  |               |
| PropType=      |                         | 0d0      |               |
| nECoord=       | 0                       |          |               |
|                |                         |          |               |
| TElt=          | None 1.46d+00           |          |               |
|                |                         |          |               |
|                | Geometric 6             |          |               |
|                |                         |          |               |
|                | -4.610703399D-01        |          |               |
|                | 8.660251521D-01         |          |               |
+================+=========================+==========+===============+
|                | -7.985963264D-01        | 3.35027  |               |
|                | -5.000004359D-01        | 5684D-01 |               |
|                |                         | 0d0 0d0  |               |
|                | 3.868566261D-01         | 0d0      |               |
|                | -1.153428518D-17        |          |               |
|                |                         | 9.22139  |               |
|                |                         | 8760D-01 |               |
|                |                         | 0d0 0d0  |               |
|                |                         | 0d0      |               |
+----------------+-------------------------+----------+---------------+
| 0d0 0d0 0d0    |                         |          |               |
| -4             |                         |          |               |
| .610703399D-01 |                         |          |               |
| 8              |                         |          |               |
| .660251521D-01 |                         |          |               |
| 1              |                         |          |               |
| .934284817D-01 |                         |          |               |
+----------------+-------------------------+----------+---------------+
|                | 0d0 0d0 0d0             |          | 3.            |
|                | -7.985963264D-01        |          | 350275684D-01 |
|                | -5.000004359D-01        |          |               |
|                |                         |          | 9.            |
|                | 0d0 0d0 0d0             |          | 221398760D-01 |
|                | 3.868566261D-01         |          |               |
|                | -1.153428518D-17        |          |               |
+----------------+-------------------------+----------+---------------+
| iElt= EltName= | 8                       |          |               |
| Element=       |                         |          |               |
| Surface=       | SecMir Reflector Conic  |          |               |
| KrElt= KcElt=  |                         |          |               |
| psiElt=        | -2.229620353D-01        |          |               |
| VptElt=        |                         |          |               |
| RptElt=        | -1.174529970D+00        |          |               |
| IndRef=        |                         |          |               |
| Extinc= nObs=  | 0d0 0d0 1d0             |          |               |
| ApType= zElt=  |                         |          |               |
| PropType=      | 0d0 0d0 1.353d+00       |          |               |
| nECoord=       |                         |          |               |
|                | 0d0 0d0 1.353d+00       |          |               |
|                |                         |          |               |
|                | 1d0 1d22                |          |               |
|                |                         |          |               |
|                | 0                       |          |               |
|                |                         |          |               |
|                | None 1.07d-01           |          |               |
|                |                         |          |               |
|                | Geometric               |          |               |
|                |                         |          |               |
|                | -6                      |          |               |
+----------------+-------------------------+----------+---------------+
| iElt= EltName= | 9                       |          |               |
| Element=       |                         |          |               |
| Surface=       | RefSurf Reference Conic |          |               |
| KrElt= KcElt=  |                         |          |               |
| psiElt=        | -2.6d+00                |          |               |
| VptElt=        |                         |          |               |
| RptElt=        | 0d0                     |          |               |
| IndRef=        |                         |          |               |
| Extinc= nObs=  | 0d0 0d0 -1d0            |          |               |
| ApType= zElt=  |                         |          |               |
| PropType=      | 0d0 0d0 1.291d+00       |          |               |
| nECoord=       |                         |          |               |
|                | 0d0 0d0 1.291d+00       |          |               |
|                |                         |          |               |
|                | 1d0 1d22                |          |               |
|                |                         |          |               |
|                | 0                       |          |               |
|                |                         |          |               |
|                | None 2.6d+00            |          |               |
|                |                         |          |               |
|                | FarField                |          |               |
|                |                         |          |               |
|                | -6                      |          |               |
+----------------+-------------------------+----------+---------------+
| iElt= EltName= | 10                      |          |               |
| Element=       |                         |          |               |
| Surface=       | FclPlane FocalPlane     |          |               |
| KrElt= KcElt=  | Flat                    |          |               |
| psiElt=        |                         |          |               |
| VptElt=        | 0d0                     |          |               |
| RptElt=        |                         |          |               |
| IndRef=        | -1.d+38                 |          |               |
| Extinc= nObs=  |                         |          |               |
| ApType=        | 0d0 0d0 -1d0            |          |               |
|                |                         |          |               |
|                | 0d0 0d0 -1.309d+00      |          |               |
|                |                         |          |               |
|                | 0d0 0d0 -1.309d+00      |          |               |
|                |                         |          |               |
|                | 1d0 1d22                |          |               |
|                |                         |          |               |
|                | 0                       |          |               |
|                |                         |          |               |
|                | None                    |          |               |
+----------------+-------------------------+----------+---------------+
| **196**        | **Modeling and Analysis |          |               |
|                | for Controlled Optical  |          |               |
|                | Systems**               |          |               |
+----------------+-------------------------+----------+---------------+

    zElt= 1.d+19 PropType= Geometric nECoord= -6
    nOutCord= 5
    Tout= 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0

### Near-Field Propagation Example

#### coroExample.jou

    % MACOS demonstration macro
    % This system is a simple coronograph. The demo
    % generates a bright on-axis (occulted) image, then
    % adds a dim, nearby but slightly off-axis image.
    % Run this with macos512!

+-------------------------------+---+----------------------------------+
| old coroExample;              | \ | sketch system                    |
|                               | % |                                  |
| draw 0 16 XZ                  |   |                                  |
+===============================+===+==================================+
| stretch sqrt                  | \ | undo occulting mask              |
|                               | % |                                  |
| mod ngridpts=128              |   |                                  |
|                               |   |                                  |
| nObs(6)=0                     |   |                                  |
+-------------------------------+---+----------------------------------+
| q                             | \ | pupil function                   |
|                               | % |                                  |
| int 1                         |   |                                  |
+-------------------------------+---+----------------------------------+
| int 6                         | \ | unocculted PSF                   |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| pix 6 64 ;                    | \ | blowup of PSF                    |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| mod ngridpts=128              | \ | redo occulting mask              |
|                               | % |                                  |
| nObs(6)=1                     |   |                                  |
+-------------------------------+---+----------------------------------+
| q                             | \ | blowup of occulted PSF           |
|                               | % |                                  |
| pix 6 64 ;                    |   |                                  |
+-------------------------------+---+----------------------------------+
| stretch linear                | \ | pupil                            |
|                               | % |                                  |
| int 10                        |   |                                  |
+-------------------------------+---+----------------------------------+
| int 11                        | \ | pupil at Lyot stop               |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| stretch sqrt int 16           | \ | on-axis PSF at detector          |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| win Tout 1.392542d-6 0,0 0,0  | \ | set window for composite image   |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| compose 16 64 1.392542E-06    | \ | build composite image            |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| add yes                       | \ | blowup of on-axis PSF            |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| stop elt 11 0,0               | \ | set system stop at Lyot stop     |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| ffp 6 2d-5,2d-5 yes           | \ | move new image off occulting     |
|                               | % | mask                             |
+-------------------------------+---+----------------------------------+
| ors 5 yes                     | \ | reset return surface             |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| ors 7 yes                     | \ | reset return surface             |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| fex 15 yes                    | \ | reset return surface             |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+
| mod flux=0.01                 | \ | 100x dimmer object               |
|                               | % |                                  |
+-------------------------------+---+----------------------------------+

    q stretch log10
    int 6 % new dim image at occulting mask stretch sqrt
    add yes % add dim image to composite image

#### coroExample.in

    ChfRayDir= 0d0 0d0 1d0 ChfRayPos= 0d0 0d0 -2.d+01 zSource= 1.d+30
    IndRef= 1d0
    Extinc= 0d0 Wavelen= 1.d-06
    Flux= 1d0 GridType= Circular Aperture= 4d0 Obscratn= 0d0 nGridpts= 512
    xGrid= 1d0 0d0 0d0 yGrid= 0d0 1d0 0d0 nElt= 16
    iElt= 1
    EltName= SMObs_1 Element= Obscuring Surface= Conic
    KrElt= -1.d+30 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 -5.4d+00 RptElt= 0d0 0d0 -5.4d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 3
    ObsType= Circle
    ObsVec= 5d-01 0d0 0d0 ObsType= Rectangle
    ObsVec= -5d-02 5d-02 -2d0 2d0 ObsType= Rectangle
    ObsVec= -2d0 2d0 -5d-02 5d-02 xObs= 1d0 0d0 0d0
    ApType= Circular ApVec= 2d0 0d0 0d0 zElt= 1.d+30
    PropType= Geometric nECoord= -6
    iElt= 2
    EltName= PM_1 Element= Reflector Surface= Conic
    KrElt= -1.08d+01
    KcElt= -1d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 0d0 RptElt= 0d0 0d0 0d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 5.4d+00
    PropType= Geometric nECoord= -6
    iElt= 3
    EltName= SM_1 Element= Reflector
    Surface= Conic
    KrElt= -3.526787501D+00 KcElt= -2.670556022D+00
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 -4.061145902D+00
    RptElt= 0d0 0d0 -5.4d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1.338854098D+00
    PropType= Geometric nECoord= -6
    iElt= 4
    EltName= ret1_1 Element= Return Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 1.5d+00 RptElt= 0d0 0d0 1.5d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 0d0 PropType= Geometric nECoord= -6
    iElt= 5
    EltName= ret2_1 Element= Return Surface= Conic
    KrElt= -6.329334641D+00
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 -4.829334641D+00 RptElt= 0d0 0d0 -4.829334641D+00
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 6.329334641D+00
    PropType= NF1 nECoord= -6
    iElt= 6
    EltName= CoroMask Element= Obscuring Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 1.5d+00 RptElt= 0d0 0d0 1.5d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 1
    ObsType= Circle
    ObsVec= 1.5d-05 0d0 0d0 xObs= 1d0 0d0 0d0
    ApType= Circular ApVec= 2d0 0d0 0d0 zElt= 1.d+30
    PropType= NF2 nECoord= -6
    iElt= 7
    EltName= ret2_2 Element= Return Surface= Conic
    KrElt= -6.329334641D+00
    KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 7.829334641D+00 RptElt= 0d0 0d0 7.829334641D+00
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= -6.329334641D+00
    PropType= Geometric nECoord= -6
    iElt= 8
    EltName= ret1_2 Element= Return Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 1.5d+00 RptElt= 0d0 0d0 1.5d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 0d0 PropType= Geometric nECoord= -6
    iElt= 9
    EltName= SM_2 Element= Reflector Surface= Conic
    KrElt= -3.526787501D+00 KcElt= -2.670556022D+00
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 7.061145902D+00
    RptElt= 0d0 0d0 8.4d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1.338854098D+00
    PropType= Geometric nECoord= -6
    iElt= 10
    EltName= PM_2 Element= Reflector
    Surface= Conic KrElt= -1.08d+01
    KcElt= -1d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 3d0 RptElt= 0d0 0d0 3d0 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 5.4d+00
    PropType= Geometric nECoord= -6
    iElt= 11
    EltName= LyotStop Element= Obscuring Surface= Conic
    KrElt= -1.d+30 KcElt= 0d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 8.4d+00 RptElt= 0d0 0d0 8.4d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 3
    ObsType= Circle
    ObsVec= 5.d-01 0d0 0d0 ObsType= Rectangle
    ObsVec= -2.5d-01 2.5d-01 -2d0 2d0
    ObsType= Rectangle
    ObsVec= -2d0 2d0 -2.5d-01 2.5d-01
    xObs= 1d0 0d0 0d0 ApType= Circular ApVec= 2d0 0d0 0d0
    zElt= 1.d+30 PropType= Geometric nECoord= -6
    iElt= 12
    EltName= PM_3 Element= Reflector Surface= Conic
    KrElt= -1.08d+01
    KcElt= -1d0
    psiElt= 0d0 0d0 -1d0 VptElt= 0d0 0d0 1.38d+01 RptElt= 0d0 0d0 1.38d+01 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None zElt= 5.4d+00
    PropType= Geometric nECoord= -6
    iElt= 13
    EltName= SM_3 Element= Reflector
    Surface= Conic
    KrElt= -3.526787501D+00 KcElt= -2.670556022D+00
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 9.738854098D+00
    RptElt= 0d0 0d0 8.4d+00 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1.338854098D+00
    PropType= Geometric nECoord= -6
    iElt= 14
    EltName= ret1_3 Element= Return Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 1.53d+01 RptElt= 0d0 0d0 1.53d+01 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 0d0 PropType= Geometric nECoord= -6
    iElt= 15
    EltName= ret2_3 Element= Return Surface= Conic
    KrElt= -6.329334641D+00
    KcElt= 0d0
    psiElt= 0d0 0d0 1d0
    VptElt= 0d0 0d0 8.970665359D+00 RptElt= 0d0 0d0 8.970665359D+00
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 6.329334641D+00
    PropType= FarField nECoord= -6
    iElt= 16
    EltName= foc_pln Element= FocalPlane Surface= Flat
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 1d0 VptElt= 0d0 0d0 1.53d+01 RptElt= 0d0 0d0 1.53d+01 IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 0d0 PropType= Geometric nECoord= -6
    nOutCord= 5
    Tout= 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0

### Image Simulation Example

#### ImageDemo.jou

1.  **ImageDemo.in**

```{=html}
<!-- -->
```
    ChfRayDir= 0d0 0d0 -1d0 ChfRayPos= 0d0 0d0 6.d+02 zSource= 1d22
    IndRef= 1d0
    Extinc= 0d0 Wavelen= 7.d-04
    Flux= 1.093342085D+00
    GridType= Circular Aperture= 3.8d+02 Obscratn= 0d0 nGridpts= 33
    xGrid= -1d0 2.033814221D-09 0d0 yGrid= -2.033814221D-09 -1d0 0d0
    nElt= 13
    iElt= 1
    EltName= SEC_OBS Element= Obscuring Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 4.002557400D+02 RptElt= 0d0 0d0 4.002557400D+02
    IndRef= 1d0
    Extinc= 1d22
    nObs= 2
    ObsType= Circle
    ObsVec= 4.67036d+01 0d0 0d0
    ObsType= Rectangle
    ObsVec= 0d0 1.75d+02 2.39d+01 2.69d+01 xObs= 9.896513868D-01 1.434926220D-01 0d0
    ApType= Circular
    ApVec= 1.75d+02 0d0 0d0 zElt= 1d22
    PropType= Annulus nECoord= -6
    iElt= 2
    EltName= SPIDER_2 Element= Obscuring Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 4.002557399D+02 RptElt= 0d0 0d0 4.002557400D+02
    IndRef= 1d0
    Extinc= 1d22
    nObs= 1
    ObsType= Rectangle
    ObsVec= 0d0 1.75d+02 2.39d+01 2.69d+01 xObs= -6.190939493D-01 7.853169309D-01 0d0
    ApType= Circular
    ApVec= 1.75d+02 0d0 0d0 zElt= 1d22
    PropType= Annulus nECoord= -6
    iElt= 3
    EltName= SPIDER_3 Element= Obscuring Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 4.002557398D+02 RptElt= 0d0 0d0 4.002557400D+02
    IndRef= 1d0
    Extinc= 1d22
    nObs= 1
    ObsType= Rectangle
    ObsVec= 0d0 1.75d+02 2.39d+01 2.69d+01 xObs= -3.705574375D-01 -9.288095529D-01 0d0
    ApType= Circular
    ApVec= 1.75d+02 0d0 0d0 zElt= 1d22
    PropType= Annulus nECoord= -6

+--------+----------------------------------------------+--------------+
| iElt=  | 4                                            |              |
| El     |                                              |              |
| tName= | PRI_MIR                                      |              |
| El     |                                              |              |
| ement= | Reflector                                    |              |
|        |                                              |              |
| Su     | Zernike                                      |              |
| rface= |                                              |              |
+========+==============================================+==============+
| KrElt= | -9.337886403D+02                             |              |
+--------+----------------------------------------------+--------------+
| KcElt= | -1.009724001D+00                             |              |
+--------+----------------------------------------------+--------------+
| p      | 0d0 0d0 1d0                                  |              |
| siElt= |                                              |              |
+--------+----------------------------------------------+--------------+
| V      | 0d0 0d0 0d0                                  |              |
| ptElt= |                                              |              |
+--------+----------------------------------------------+--------------+
| R      | 0d0 0d0 0d0                                  |              |
| ptElt= |                                              |              |
+--------+----------------------------------------------+--------------+
| I      | 1d0                                          |              |
| ndRef= |                                              |              |
+--------+----------------------------------------------+--------------+
| E      | 1d22                                         |              |
| xtinc= |                                              |              |
+--------+----------------------------------------------+--------------+
| Zer    | 0d0 0d0 0d0 -1.415805330D-04                 | 5.4          |
| nCoef= | -3.404653098D-05                             | 12244128D-06 |
+--------+----------------------------------------------+--------------+

    1.464334396D-05 -5.046596565D-05 1.587508743D-05 -3.536721061D-
    05 -1.951889059D-06 -2.971194787D-05
    1.561541879D-05 -8.027617601D-06 8.191919142D-06 0d0 0d0 0d0
    0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0
    pMon= 0d0 0d0 0d0 xMon= 1d0 0d0 0d0 yMon= 0d0 1d0 0d0 zMon= 0d0 0d0 1d0 lMon= 1.75d+02 nObs= 0
    ApType= None
    zElt= 1d22

+----------+----------------------------------------+------------------+
| P        | Geometric                              |                  |
| ropType= |                                        |                  |
| nECoord= | -6                                     |                  |
+==========+========================================+==================+
| iElt=    | 5                                      | 4.002602400D+02  |
| EltName= |                                        |                  |
| Element= | SEC_MIR                                |                  |
| Surface= |                                        |                  |
| KrElt=   | Reflector Conic                        |                  |
| KcElt=   |                                        |                  |
| psiElt=  | -1.565935400D+02                       |                  |
| VptElt=  |                                        |                  |
|          | -1.932283999D+00                       |                  |
|          |                                        |                  |
|          | 0d0 0d0 1d0                            |                  |
|          |                                        |                  |
|          | -7.887151838D-04 -6.647047308D-05      |                  |
+----------+----------------------------------------+------------------+
| RptElt=  | -7.887151838D-04 -6.647047308D-05      | 4.002602400D+02  |
| IndRef=  |                                        |                  |
| Extinc=  | 1d0 1d22                               |                  |
| nObs=    |                                        |                  |
| ApType=  | 0                                      |                  |
| zElt=    |                                        |                  |
| P        | None 1d22                              |                  |
| ropType= |                                        |                  |
| nECoord= | Geometric                              |                  |
|          |                                        |                  |
|          | -6                                     |                  |
+----------+----------------------------------------+------------------+
| iElt=    | 6                                      |                  |
| EltName= |                                        |                  |
| Element= | COR_S1                                 |                  |
| Surface= |                                        |                  |
| KrElt=   | Refractor Conic                        |                  |
| KcElt=   |                                        |                  |
| psiElt=  | -4.241546000D+01                       |                  |
| VptElt=  |                                        |                  |
| RptElt=  | 0d0                                    |                  |
| IndRef=  |                                        |                  |
| Extinc=  | 0d0 0d0 1d0                            |                  |
| nObs=    |                                        |                  |
| ApType=  | 0d0 0d0 -3.988054000D+01               |                  |
| zElt=    |                                        |                  |
| P        | 0d0 0d0 -4.260596000D+01               |                  |
| ropType= |                                        |                  |
| nECoord= | 1.61783d+00                            |                  |
|          |                                        |                  |
|          | 0d0 0                                  |                  |
|          |                                        |                  |
|          | None 1d22                              |                  |
|          |                                        |                  |
|          | Geometric                              |                  |
|          |                                        |                  |
|          | -6                                     |                  |
+----------+----------------------------------------+------------------+
| iElt=    | 7                                      |                  |
| EltName= |                                        |                  |
| Element= | COR_S2                                 |                  |
| Surface= |                                        |                  |
| KrElt=   | Refractor Conic                        |                  |
| KcElt=   |                                        |                  |
| psiElt=  | -3.15595d+01                           |                  |
| VptElt=  |                                        |                  |
| RptElt=  | 0d0                                    |                  |
| IndRef=  |                                        |                  |
| Extinc=  | 0d0 0d0 -1d0                           |                  |
| nObs=    |                                        |                  |
| ApType=  | 0d0 0d0 -4.138168000D+01               |                  |
| zElt=    |                                        |                  |
| P        | 0d0 0d0 -4.260596000D+01               |                  |
| ropType= |                                        |                  |
| nECoord= | 1.657863000D+00                        |                  |
|          |                                        |                  |
|          | 0d0 0                                  |                  |
|          |                                        |                  |
|          | None 1d22                              |                  |
|          |                                        |                  |
|          | Geometric                              |                  |
|          |                                        |                  |
|          | -6                                     |                  |
+----------+----------------------------------------+------------------+
| iElt=    | 8                                      |                  |
| EltName= |                                        |                  |
| Element= | COR_S3                                 |                  |
| Surface= |                                        |                  |
| KrElt=   | Refractor Conic                        |                  |
|          |                                        |                  |
|          | -5.533237600D+02                       |                  |
+----------+----------------------------------------+------------------+

    KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 -4.533138000D+01 RptElt= 0d0 0d0 -4.260596000D+01
    IndRef= 1d0
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 9
    EltName= FLT_S1 Element= Refractor Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 -4.763008000D+01
    RptElt= 0d0 0d0 -4.9d+01 IndRef= 1.535018000D+00
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 10
    EltName= FLT_S2 Element= Refractor Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 -5.163058000D+01
    RptElt= 0d0 0d0 -4.9d+01 IndRef= 1d0
    Extinc= 0d0
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 11
    EltName= IMG_RET Element= Return Surface= Conic
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 0d0 0d0 -5.175758000D+01 RptElt= 0d0 0d0 -5.175758000D+01
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    iElt= 12
    EltName= EXIT_RET Element= Return Surface= Conic
    KrElt= -6.874662920D+01
    KcElt= 0d0
    psiElt= 1.095041158D-04 -5.786697499D-06 -9.999999940D-01 VptElt= -7.904014798D-05 -6.661258889D-06 1.698904900D+01 RptElt= -7.904014798D-05 -6.661258889D-06 1.698904900D+01
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 6.874662920D+01
    PropType= FarField nECoord= -6
    iElt= 13
    EltName= IMG_PL Element= FocalPlane Surface= Flat
    KrElt= -1d22 KcElt= 0d0
    psiElt= 0d0 0d0 -1d0
    VptElt= 7.421487214D-01 -4.232318881D+00 -5.175758000D+01 RptElt= 0d0 0d0 -5.175758000D+01
    IndRef= 1d0
    Extinc= 1d22
    nObs= 0
    ApType= None
    zElt= 1d22 PropType= Geometric nECoord= -6
    nOutCord= 5
    Tout= 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 0d0 1d0

### S-MACOS Example Files

#### smacosvars.cmn

    C*********************************************************************** C Begin file smacosvars.cmn
    C This file defines call-line variables and common blocks that give C user programs direct access to MACOS variables. The form of the C SMACOS call line is:

+---+----+------------------------------------------------------------+
| C |    | CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,              |
+===+====+============================================================+
| C | &  | OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat,          |
+---+----+------------------------------------------------------------+
| C | &  | PixArray,Cmatrix)                                          |
+---+----+------------------------------------------------------------+

    C The param.cmn file defines array sizes for SMACOS call-line and
    C internal MACOS variables. The param.cmn file should be copied from C param128.cmn, param256.cmn, etc. -- whichever size of SMACOS you C are using:
    INCLUDE 'param.cmn'
    C These are the SMACOS call-line variables:
    CHARACTER*132 command CHARACTER*32 CARG(9) REAL*8 DARG(9) INTEGER IARG(9) LOGICAL LARG
    REAL*4 RARG(9), OutRaySpot(mRay, 4), NomRayDir(2, mRay), & RMSWFE, NomOPDMat(mpts, mpts), PertOPDMat(mpts, mpts), & PixArray(mPix, mPix), Cmatrix(7, mCm, bRay)
    C These are other SMACOS common blocks, that you may use or comment C out as you wish. The ifEcho flag, if set to .TRUE., provides for C full printout of SMACOS operations. If set to .FALSE., reduced
    C printout is provided.
    LOGICAL ifEcho COMMON /SCIO/ifEcho
    INTEGER npts,NoiseSeed(2) REAL*4 CntrSpot(3)
    COMMON /UserCommon/npts,NoiseSeed,CntrSpot
    C These includes provide access to internal MACOS variables. They can C be commented out if you do not plan to use them:
    INCLUDE 'elt.cmn' INCLUDE 'src.cmn'
    C These includes provide access to internal MACOS scratch arrays. They C should be commented out if you do not plan to use them:
    LOGICAL L1(md2)
    INTEGER DrawEltVec(mDrawElt,mDrawRay) REAL*4 R1(mdttl,mdttl),RV1(md2),RV2(md2) REAL*8 D1(mdttl,mdttl),D2(mdttl,mdttl) REAL*8 DV1(md2),DV2(md2)
    COMMON /SCRATCH/ D1,D2
    EQUIVALENCE (D1(1,1),RV1(1)),(D1(1,1),DV1(1)),(D2(1,1),DV2(1)), & (D1(1,1),R1(1,1)),(D2(1,1),RV2(1))
    c END of smacosvars.cmn C***********************************************************************

#### smacosExample.f

    C*********************************************************************** C Begin file smacosExample.f
    IMPLICIT NONE
    C SMACOS variable declarations INCLUDE 'smacosvars.cmn'
    C These are local variables:
    INTEGER REAL*4 REAL*8
    C Load prescription file
    command='OLD' CARG(1)='cassExitPupil'
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    C Set plot type
    command='GRAY'
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    C Set system stop
    command='STOP' CARG(1)='ELT' IARG(1)=3 DARG(1)=0d0 DARG(2)=0d0
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    C Set exit pupil
    command='FEX' IARG(1)=5 CARG(1)='YES'
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    C Plot OPD at exit pupil after resetting ray grid size
    nGridPts=128 command='RESET'
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    command='OPD' IARG(1)=5
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    C Plot detector intensity after resetting ray grid size
    nGridPts=96 command='RESET'
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    command='INTENSITY' IARG(1)=6
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    C Plot detector using a log10 stretch
    command='STRETCH' CARG(1)='LOG10'
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    command='INTENSITY' IARG(1)=6
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    C Quit!
    command='QUIT'
    CALL SMACOS(command,CARG,DARG,IARG,LARG,RARG,
    & OutRaySpot,NomRayDir,RMSWFE,NomOPDMat,PertOPDMat, & PixArray,Cmatrix)
    STOP END
    C End file smacosExample.f C***********************************************************************
