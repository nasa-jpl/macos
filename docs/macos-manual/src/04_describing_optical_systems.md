## SECTION 4: Describing Optical Systems

*Systems*

A first MACOS session typically begins by invoking the NEW command. This begins a dialog whose function is to define the data used by MACOS in its computations. This data comprises the optical prescription, and is saved by MACOS as an “.in-file,” which can be edited and reused later using the OLD command. The contents of the .in-file can be altered within a MACOS session using MODIFY, PERTURB or other commands, and saved as a new .in-file using the SAVE command. The .in-files can also be created and edited off-line using a text editor. The commands covered in this section are:

                    NEW
                    DRAw
                    OLD (LOAd)
                    SUMmarize
                    STatus
                    REset
                    SHOw
                    MODify
                    BEAm
                    SAVe
                    VALidate

### Specifying the Light Source

The NEW dialog begins by asking the user for information about the light source, such as where it is located, in what direction it is oriented, and what its optical properties are. After the source is defined, similar questions are asked about each optical element in turn. This section is intended to help the user provide the appropriate answers to these questions. It is also intended to define MACOS variable names, which are used in the

.in-files to specify aspects of the prescription.

MACOS represents light sources as *beams*, consisting of a ray bundle and an associated complex amplitude matrix (or matrices). Beams can be collimated or focused, as per Figure 11. In addition, the illumination profile of a beam can be specified as uniform, gaussian, dipole or cosine-to-a-power, using the BEAM command discussed in a later section. Collimated beams are composed of a bundle of parallel rays, originating at a plane in front of the entrance aperture. Collimated beams are used to represent the light received from a point-source at an infinite distance from the system. They are also used when collimated light sources, such as some lasers, are to be modeled. Focused beams are bundles of rays converging to or diverging from a point in space and represent the light from a nearby point source. Focused beams are defined on a spherical surface whose center of curvature is at the source point.

The properties that must be specified by the user include:

-   Source position and direction, defined by the position and orientation of the “chief” ray (really the gut ray – the ray that lies in the center of the beam)

-   Aperture size (collimated source) or angle (point source)

```{=html}
<!-- -->
```
-   Aperture shape (circular, square, hexagonal, pie, flower)

-   Illumination profile (flat, gaussian, dipole or cosine-to-a-power)

-   Wavelength

-   Flux through the defined aperture

-   Polarization state (if using polarization ray-trace analysis functions)

-   Segmentation (unsegmented, hexagonal segments, or “pie” segments)

-   Coordinates to define the orientation of the beam to the optical elements

#### Source Position and Orientation

In MACOS, a *light source* is a surface but not an element. The surface is defined by a single wavefront emitted from a single point source. The zSource parameter is the distance from the origin of the light to the light source surface, as defined by the location of the chief ray, ChfRayPos (see Figure 11). Collimated sources are specified by zSource=1d22. For point sources, the convention is that the sign of zSource is positive for a converging beam and negative for an expanding beam.

![](media/fig_rescued_11.png)

**FIGURE 11** Light sources

The location and orientation of the light source is determined from the *chief ray*. The MACOS “chief ray” is the central or gut ray of the beam, as discussed in [Section ?](#geometric-optics-tracing-rays). It is the center ray of the incident beam bundle. It is normally a true chief ray as well, passing through the center of the system stop. This condition is enforced by the user, using the STOP or CENTER commands, as described later.

To define the source, the user must specify in the global coordinate system:

-   Chief ray position vector

-   Chief ray direction vector

These are set by the 3-vector inputs ChfRayPos and ChfRayDir. If the source is collimated, then every ray in the beam starts out with its direction equal to ChfRayDir. The position of each ray is calculated from ChfRayPos, the shape and type of aperture, and the specified ray density. If the source is a point source, appropriate position and direction is computed so that each ray starts from the surface of a sphere centered at the source point.

NOTE: Converging or diverging sources **should not be located at the source point**, or, for example, divide-by-zero errors will result when performing diffraction analysis.

#### Refractive Index and Extinction Coefficient at the Source

The refractive index is, in general, a complex number. The real part is the conventional index, giving the speed of light in the particular medium, and so defining the direction

changes a ray experiences at refractive interfaces. The imaginary part is the extinction coefficient, which defines the attenuation of the field in the medium. Only the index (real part) is significant in MACOS scalar ray-tracing calculations. The extinction coefficient is used only when polarization ray-trace is selected, by use of the POL command.

-   The real part of the refractive index of the source medium is specified by IndRef. If ChfRayPos is in a vacuum, then IndRef=1. For many applications it is sufficient to approximate the index of refraction of air as 1, although the exact value depends on the pressure and temperature.

-   The negative of the imaginary part of the refractive index, the extinction coefficient, of the source medium is specified by Extinc. For lossless media such as vacuum and air, set ExtInc=0. For perfect metals, set ExtInc=1d22. This parameter is used only when polarization ray-tracing is turned on.

The user must specify the index explicitly for the source only. For other elements, particular glasses can be specified, and the index will be varied as a function of the wavelength of the beam. The source IndRef and Extinc are not changed automatically with changes in wavelength. The user can implement such changes as needed using the MODIFY command.

#### Wavelength

The wavelength is specified by the wavelen variable. For most calculations, MACOS uses one wavelength at a time. Multiple-wavelength images can be generated using the COMPOSE command to define detector parameters, using use the MODify command to change the wavelen (and Flux) variables (see [Section ?](#modify)), then using the ADD command to separately propagate and accumulate different wavelength and flux level images at the detector. This procedure can be automated by using “.filt” filter files and the MULtispec command.

#### Flux

The total flux of the source illumination through the aperture area (defined by the Aperture, Obscuration and GridType parameters) is specified by Flux.This parameter is used in MACOS to facilitate the radiometric modeling of systems. It is useful when constructing composite images of multiple sources or at multiple wavelengths. This procedure can be automated by using “.filt” filter files and the MULtispec command.

#### Aperture

The aperture size is set by the Aperture parameter, whose interpretation is different for collimated or point sources. If the source is collimated (zSource=1d22), Aperture is the *diameter* of the input beam. If it is a point source(zSource\<1d22), Aperture is the full *angle* (entered in radians like all angles in MACOS) subtended by the two marginal rays (see Figure 12). Thus, the numerical aperture of the system is

$$
\mathrm{N.A.} \;=\; n \,\sin(\mathrm{Aperture})
\qquad (4.1)
$$

where n is the index of refraction.

Obscurations can be defined in the source beam if desired. Obscuration size is specified in the same way as the aperture, with the parameter Obscratn setting the diameter or angle of a *central* obscuration.

![](media/fig_rescued_12.png)

**FIGURE 12** Aperture sizing

NOTE: This aperture obscuration option can be useful, but is not the best way to represent most beam-train obscurations. It is better to explicitly place obscurations on the appropriate elements as discussed in [Section ?](#obscurations-and-apertures).

The entrance aperture can have several shapes, as shown in Figure 13, depending on whether the system is monolithic – in which all optics are of one piece – or segmented – in which at least one mirror is split into separate segments. For monolithic systems, the aperture shape can be round or square, as set by the GridType parameter. For segmented systems the options are hexagonal array, round with 6 pie-shaped outer segments, or round with an arbitrary number of outer segments. For computational efficiency, the segmentation is defined at the light source.

![](media/fig_rescued_13.png)

**FIGURE 13** Aperture shapes

It is not necessary to segment the aperture to use lenslet arrays, as discussed later.

NOTE: Registration of a segmented aperture on a segmented mirror may shift when the beam incidence angle changes, or when intervening optics are perturbed. This can cause rays that nominally strike one mirror segment to move to another. In these conditions, the user should check to ensure proper registration.

#### Ray Grids

The density of rays across the beam and the orientation of the ray grid are set by three

parameters:

-   Number of ray grid points across the aperture is set by nGridPts

-   Ray grid coordinate vectors xGrid and yGrid determine the global orientation of the ray grid and so the beam coordinates

The ray grid is laid out on the surface that defines the light source, according to a coordinate system for which zGrid equals the negative chief ray direction ChfRayDir and xGrid is as specified by the user (see Figure 14). The number of rays across the aperture is set by nGridPts.

![](media/fig_rescued_14.png)

**FIGURE 14** Ray grid layout

### Example Light Source

MACOS obtains the information it needs to build an optical prescription (an “.in-file” dataset) through the NEW file definition dialog. The first questions in this dialog are about the light source and the latter questions pertain to the optical elements.

For the first example, we start the prescription for a Cassegrain telescope, Cassegrain.in, shown in Figure 15. The light source is a point source at infinity so the wave-front incident on the telescope is collimated. The example dialog is MACOS's own output, run from the answers shown after each prompt; the variable each prompt sets is listed in [Section ?](#summary-of-prescription-variables).

![](media/fig_rescued_15.png)

**FIGURE 15** Cassegrain telescope layout

<!-- transcript: s04_new_source -->

After the elements are entered (Example Telescope, below), SHOW with element 0 prints the source block as the prescription file Cassegrain.in holds it:

<!-- transcript: s04_new_source_show -->

The Cassegrain.in example is finished in [Section ?](#example-telescope).

Now consider the case of a non-collimated input beam, zSource not equal to infinity. A Cassegrain telescope has been turned into a beam projector. The light source is now located where the detector was previously located as shown in Figure 16.

<!-- transcript: s04_zsourc -->

![](media/fig_rescued_16.png)

**FIGURE 16** Non-collimated light source example

After the six elements are entered, SHOW with element 0 prints the source block as the prescription file zsourc.in holds it:

<!-- transcript: s04_zsourc_show -->

Notice that the Aperture entry has changed significantly, since Aperture now represents an angle, instead of the diameter of the input beam. The ChfRayPos is now a small distance in front of the source, so zSource is now negative to indicate an diverging beam.

### Specifying Optics

The next part of the new dialog is the specification of the optical elements. The first step in defining the beam train is actually the last step in the source dialog: specifying the number nElt of optical elements in the beam train. The number of elements is equal to the number of active optical surfaces, plus any reference (dummy) surfaces, plus 2 return surfaces for an exit pupil, if desired.

The new dialog continues by prompting the user for all data needed to define each optic. Optics are specified one at a time, in the sequence in which they are encountered by the light from the source (except nonsequential and segmented optics). Some optics, such as mirrors, are realized by a single *element*, or optically active surface. Other optics require

two or more elements. Lenses require two refracting elements, one for each surface (see the example in [Section ?](#refractors)). Corner cubes require three non-sequential reflecting elements. Each element in multiple-element optics is entered separately.

The first step in specifying an element is to define its *type* (e.g., Reflector) and its *surface type* (e.g., Conic). Table 5 lists the elements currently supported by MACOS, and Table 6 lists the surface types. Further data required to define the surfaces may be requested, depending on specified surface type.

**TABLE 5.** Element types

+---------------------+-----------------+-----------------------------+
| **Element Type      | **Allowed       | **Description**             |
| (**Element=..**)**  | Surface Types** |                             |
+=====================+=================+=============================+
| Reflector           | all             | reflective surfaces (e.g.,  |
|                     |                 | mirrors)                    |
+---------------------+-----------------+-----------------------------+
| Refractor           | all             | refractive surfaces (e.g.,  |
|                     |                 | lenses)                     |
+---------------------+-----------------+-----------------------------+
| LensArray           | conics or       | lens arrays                 |
|                     | aspheres        |                             |
+---------------------+-----------------+-----------------------------+
| Grating             | all             | reflective diffraction      |
|                     |                 | gratings                    |
+---------------------+-----------------+-----------------------------+
| HOE                 | all             | reflective holographic      |
|                     |                 | optical elements            |
+---------------------+-----------------+-----------------------------+
| NSReflector         | conics          | non-sequential reflective   |
|                     |                 | surfaces                    |
+---------------------+-----------------+-----------------------------+
| NSRefractor         | conics          | non-sequential refractors   |
+---------------------+-----------------+-----------------------------+
| Segment             | all             | segmented mirrors           |
+---------------------+-----------------+-----------------------------+
| ObsSrf              | conics          | obscuring surfaces          |
+---------------------+-----------------+-----------------------------+
| FocalPln            | flats           | focal planes                |
+---------------------+-----------------+-----------------------------+
| RefSrf              | conics          | reference surfaces          |
+---------------------+-----------------+-----------------------------+
| ReturnSrf           | conics          | return surfaces             |
+---------------------+-----------------+-----------------------------+

A non-sequential group (`NSReflector`/`NSRefractor` elements) may carry an
optional `NSCount=` keyword on each member: the **hit budget**, i.e. the
number of surface hits a ray may take inside that group before it is forced
on to the next sequential element. `NSCount=` absent or `0` means
**unlimited** (needed where a ray legitimately bounces among the members,
as in a corner cube); `NSCount= 1` is the single-hit case of a segmented
mirror. The budget is authoring intent — set it to `1` on a grid-surfaced
non-sequential group, where an unbounded search can send a ray's
surface-solve bracket to divergence.

**TABLE 6.** Surface types

+--------------+-------------------------------------------------------+
| **Surface    | **Surface Description**                               |
| Type         |                                                       |
| (**Sur       |                                                       |
| face=..**)** |                                                       |
+==============+=======================================================+
| Flat         | Planar                                                |
+--------------+-------------------------------------------------------+
| Conic        | Conicoid of revolution. Precise shape is defined by   |
|              | conic constant                                        |
|              |                                                       |
|              | KcElt and radius KrElt parameters (or fElt and eElt)  |
+--------------+-------------------------------------------------------+
| Toric        | Toric or cylindrical shapes, defined by separate      |
|              | radius parameters in x and y directions.              |
+--------------+-------------------------------------------------------+
| Anamorphic   | “Saddle” shapes, defined by separate conic and radius |
|              | parameters in x and y directions.                     |
+--------------+-------------------------------------------------------+
| Aspheric     | Base conicoid, with 10th-order aspheric deviation.    |
+--------------+-------------------------------------------------------+
| Monomial     | Base conicoid, with 14th-order cartesian polynomial   |
|              | deviation.                                            |
+--------------+-------------------------------------------------------+
| Zernike      | Base conicoid, with deviation defined by Zernike      |
|              | polynomials.                                          |
+--------------+-------------------------------------------------------+
| Interpolated | Base conicoid, with deviation defined by X-Y-Z data   |
|              | triplets                                              |
+--------------+-------------------------------------------------------+
| GridData     | Base conicoid, with deviation defined by square       |
|              | matrix pixels                                         |
+--------------+-------------------------------------------------------+
| User-defined | Deformable mirrors defined using influence functions. |
+--------------+-------------------------------------------------------+

The second step is to define element location, orientation and rotation point, as discussed in [Section ?](#geometry-and-coordinate-systems):

-   Element *location* is specified by the VptElt vector. This is a three-vector in global coordinates giving the location of the vertex of the element.

-   Element *orientation* is specified by the psiElt vector. This is a unit magnitude three-vector in global coordinates giving the direction of the principal axis of the element.

-   Element *rotation point*, RptElt, defines the pivot point for rotations of the element implemented using the Perturb command. It also defines the coupling between rotation and translation for the Build command. Together with any element coordinates, Telt, RptElt governs the element response to perturbations.

The third step is to define any aperture or obscurations that might occur at the element. Apertures define the margins of the element, and can be circular or rectangular. Obscurations are circular or rectangular regions on the surface of the element where light is completely absorbed. Negative obscurations “un-absorb” the light, countering the effect of an obscuration in regions where they overlap. Obscurations and apertures can be placed anywhere on the surface.

The fourth step is to define the zElt and PropType diffraction propagation parameters, which are needed if the element is at the start or end of a diffraction propagation.

The final step is to specify element coordinates (Telt), if they are desired.

In the next several subsections, each of these steps is discussed more carefully, and examples are shown.

### Element Types

Element type and surface shape are each most critical in the determination of the properties of an optic. In this subsection we define the MACOS element types and illustrate their use through a sequence of specific examples.

#### Reflectors

Mirrors are entered as type Reflector. Reflectors can have any surface type or shape; the user will be prompted for different parameters depending on which surface type is chosen. No element-specific parameters are required for reflectors unless polarization ray-tracing is turned on. In that case, Extinc and IndRef should be carefully specified, as they determine the transmission of the beam past the reflective interface.

This example shows how to enter a flat mirror in MACOS:

<!-- transcript: s04_reflector -->

After the last prompt for an element, NEW prints the element back as the prescription file will hold it and asks for confirmation.

#### Refractors

A refracting element is the interface between two media of differing index of refraction, such as air and glass. Rays traversing the interface are, of course, bent at the interface. lenses, windows and other refracting optics are actually modeled using 2 optical elements, both of type Refractor. The first element defines the front surface of the optic and sets the index for the ray traveling inside the optic. The second element defines the back surface, and resets the index to the value outside the optic.

Some refracting elements, such as prisms, use more than two elements, and may require use of nonsequential reflectors and refractors as well – see [Section ?](#reflective-non-sequential-surfaces). Arrays of lenslets are covered in [Section ?](#lens-arrays).

As with reflectors, refractors in MACOS may have any surface shape. The index of refraction IndRef can be set directly in the prescription, or indirectly, by specifying a particular glass for the element. To set the index directly, select FixedIndex glass type and then specify the index. No other element-specific parameters are required, unless polarization ray-tracing is turned on. In that case, Extinc must be specified as well.

This example lens has a spherical front surface and a planar rear surface, as illustrated in Figure 17:

<!-- transcript: s04_refractor -->

A Refractor may carry a gradient-index model (`GradInd=` with its coefficients). Only model 5, the radial profile n(r) = n0 + n2 r² + n4 r⁴ integrated with a Runge–Kutta step, is implemented; a ray meeting any other model is lost by name (`GRIN model not implemented`, counted in the WARN summary) rather than traced with undefined values. Model 5 cannot presently be entered from a prescription file; this is an open item.

NOTE: It is not uncommon to “lose” rays, if they miss the properly defined area of a lens, that is, if the rays miss the region where the 2 surfaces of the lens are sequential. If rays extend beyond the intersection point of the two elements, the rays become undefined and are not traced further as shown in Figure 17.

This ray hits the elements out of sequence, becomes undefined

![](media/fig_rescued_17.png)

**FIGURE 17** Ray bundle must pass within surface intersection points of the refracting elements for all rays to be traced properly

#### Glass Tables

Version 2.8 of MACOS comes with a single glass table, which incorporates the Schott glass types as of May 1999 (from [http://schottglass.com).](http://schottglass.com/) This table is a text file, macos.glass, located in the MACOS data directory. See [Section ?](#macos-data-directory) for details on how to specify this directory.

The user may specify the index of glasses or other optical materials in several ways:

-   By selecting from the materials in the provided macos.glass file, as in the first part of the example in [Section ?](#refractors).

-   By using a different macos.glass file, located in the MACOS data directory. The user can easily add other materials to a private copy of macos.glass, using the MACOS_DATA environment variable to define its location.

-   By directly specifying the 6 Sellmeier glass coefficients in the prescription. These coefficients are B1, B2, B3, C1, C2, C3 from the manufacturers catalog. Note that the MACOS SAVE command saves the coefficients of each glass in the prescription file, even if the glass data was originally derived from the MACOS glass table.

-   By declaring a material to be FixedIndex, and then specifying the index, in the prescription. This approach uses a single value of the index for any wavelength of light.

The form of the glass table is simple. Each glass is defined in a single row of the file, with the first entry being the name of the glass, followed by the 6 Sellmeier coefficients. The following shows the first few entries from macos.glass:

+-----------------+------------+------------+------------+----------+
| Air 0d0 0d0 0d0 | 0d0 0d0    |            |            |          |
|                 | 0d0        |            |            |          |
+=================+============+============+============+==========+
| BaF13           | 1.51       | 1.07       | 9.39       | 4.49     |
| 1.57068590E+00  | 295152E-01 | 638026E+00 | 377487E-03 | 479866E- |
| 02              |            |            |            |          |
| 1.13934126E+02  |            |            |            |          |
+-----------------+------------+------------+------------+----------+
| BaF3            | 1.33       | 8.85       | 8.87       | 4.20     |
| 1.32064267E+00  | 572683E-01 | 521821E-01 | 798715E-03 | 290346E- |
+-----------------+------------+------------+------------+----------+
| 02              |            |            |            |          |
| 1.11729167E+02  |            |            |            |          |
+-----------------+------------+------------+------------+----------+
| BaF4            | 1.46       | 9.12       | 9.29       | 4.40     |
| 1.37492590E+00  | 354500E-01 | 643359E-01 | 554343E-03 | 426654E- |
+-----------------+------------+------------+------------+----------+
| 02              |            |            |            |          |
| 1.16607899E+02  |            |            |            |          |
+-----------------+------------+------------+------------+----------+

#### Lens Arrays

Lens arrays (a.k.a. lenticular arrays or lenslet arrays) are refracting optics. Usually, a lens array will have one flat surface and one surface figured into a mosaic of regularly distributed lenslets (Figure 18). Lens arrays are used most often in Hartmann wavefront sensors. An example adaptive optics wavefront sensor is provided in Appendix A.4.

MACOS uses a LensArray element type to model the bumpy surface of a lens array. LensArray elements are specified much as Refractor elements are, with the addition of a few new parameters. LensArray elements can have conic or aspheric lenslet shapes. Enter shape parameters that describe a single lenslet, and all the other lenslets will be automatically replicated. The location and orientation parameters VptElt and psiElt define the plane containing all of the lenslet vertices, and specify the orientation of each lenslet.

The first new parameter is LensArrayType, which can have two values: LensArray-Typ=1 for hexagonal arrays; and LensArrayTyp=2 for rectangular arrays. The second new parameter is LensArrayWidth, which is the width of a single lenslet (flat to flat for hexagonal lenslets) (Figure 18).

![](media/fig_rescued_18.png)

**FIGURE 18** Hexagonal lens array

LensArray elements use a surface coordinate system to define the location and orientation of the lenslet grid on the underlying plane. The variable pMon is the center of the coordinate system, defining the vertex position of the “central” lenslet. pMon can differ

from VptElt, though only the projection of pMon onto the plane is significant. The variables xMon, yMon, and zMon define the surface coordinate frame as shown in Figure 27.

An additional example of the use of lenslet arrays is provided in Appendix A.4.

<!-- transcript: s04_lensarray -->

#### Reflective Gratings

Diffraction gratings are surfaces that have been ruled into a highly regular array of grooves. The rulings scatter light in all directions. The periodicity of the ruling causes there to be coherence in the light scattered in particular discrete directions. These directions correspond to various “diffraction orders,” which are functions of the ruling direction, the ruling width, and the wavelength of the incident beam (ref. Hecht). The amount of power scattered into each order is a more subtle matter, depending on the the detailed shape of the rulings. MACOS provides grating elements that model the beam diffracted into any specified order, but does not predict the intensity of the diffracted beams.

MACOS gratings can be defined on any surface shape. Enter surface parameters to define the underlying element, just as you would to define a regular mirror. Three additional element-specific parameters are required to define the grating geometry. These are:

-   OrderHOE: The diffraction order you wish to trace off of the grating surface.

-   h1HOE: A unit 3-vector perpendicular to the ruling direction and the psiElt vector. This vector is rotated and translated appropriately if the element is perturbed using the MACOS Perturb command.

-   RulingWidth: The fixed distance between rules as projected to a flat plane underlying the surface. Distance between rules *along the curved surface* can vary if the surface shape is curved, which will be the case with a conic or aspheric surface type.

This is the convention of a straight-ruled grating: equidistant parallel groove planes cutting the surface, so on a curved grating the period measured along the surface grows where the surface tilts (the chord model, as in CODE V and Zemax); a flat grating is unaffected. Two consequences for analysis. The optical path difference across a grating is defined only modulo the wavelength (the engine's smooth order-m phase equals the physical groove staircase mod λ), so compare OPDs across a grating as wrapped phase or as the complex field, never as unwrapped lengths. And a grating's aperture and obscurations vignette rays as those of any other element do.

![](media/image46.png){width="2.609722222222222in" height="1.3715277777777777in"}The following example shows the new dialog used to specify a grating on a flat surface. The complete prescription, with a collimated source illuminating the grating at normal incidence, is provided in Appendix A.2. Figure 19 shows the diffracted beam at a couple of orders.

*Order=1 Order=-2*

**FIGURE 19** Diffraction grating showing diffraction orders 1 and -2.

<!-- transcript: s04_grating -->

#### Transmissive Gratings

The transmissive grating elements, like the reflective gratings of the previous section, use surfaces that have been ruled into a highly regular array of grooves to diffract light. The periodicity of the ruling causes there to be coherence in the light scattered in particular discrete directions. These directions correspond to various “diffraction orders,” which are functions of the ruling direction, the ruling width, and the wavelength of the incident beam (ref. Hecht).

MACOS grating elements model the beam diffracted into any specified order, but do not predict the intensity of the diffracted beams.

MACOS gratings can be defined on any surface shape. Enter surface parameters to define the underlying element, just as you would to define a regular mirror. Three additional element-specific parameters are required to define the grating geometry. These are:

-   OrderHOE: The diffraction order you wish to trace off of the grating surface.

-   h1HOE: A unit 3-vector perpendicular to the ruling direction and the psiElt vector. This vector is rotated and translated appropriately if the element is perturbed using the MACOS Perturb command.

-   RulingWidth: The fixed distance between rules as projected to a flat plane underlying the surface. Distance between rules *along the curved surface* can vary if the surface shape is curved, which will be the case with a conic or aspheric surface type.

#### Reflective Holographic Optical Elements (HOEs)

Holgraphic optical elements (HOEs) are diffractive elements that can be used to produce a more complex output beam than is possible with straight-ruled gratings. They are typically formed by interfering two coherent beams on the surface of the underlying optic. The surface is covered with photo-sensitive material, so the interference is recorded as areas of exposed and unexposed fringes. The piece is then etched to turn the interference photograph into permanent features on the surface of the optic. Subsequent illumination by one of the 2 reference beams produces a diffracted beam in the form and direction of the second beam. The result is a diffraction grating. MACOS currently supports reflective HOEs only.

MACOS HOEs can be defined on any surface shape. Enter surface parameters to define the underlying element, just as you would to define the shape of a regular mirror. The HOE grating is then defined using additional element-specific parameters. These parameters set up the original HOE interference geometry, so the grating on any point of the surface can be computed from the intended reference beams. The reference beams are rotated and translated appropriately if the element is perturbed using the MACOS Perturb command. The parameters are:

-   h1HOE: The first is the input reference beam point source location in x,y,z. A collimated reference beam is located at infinity (1d22). Note that the reference beams have spherical (or flat) wavefronts only.

-   h2HOE: The second is the ouput reference beam point source location in x,y,z.

-   OrderHOE: The third is the diffraction order to be used to define the grating and to be ray-traced through the rest of the optical system. MACOS currently does not allow you to trace a different order than was used to define the grating.

```{=html}
<!-- -->
```
-   WaveHOE: The last is the construction or reference wavelength. The beam is traced using the source wavelength (Wavelen) which can be different.

![](media/image48.jpeg){width="3.3716535433070867in" height="1.213779527559055in"}![](media/image49.jpeg){width="3.3791338582677164in" height="1.6716535433070867in"}The following example shows the new dialog used to specify a HOE on a parabolic surface. The complete prescription, with a collimated source illuminating the HOE at normal incidence, is provided in Appendix A.2. Figure 20 shows the output beam for both the diffracted beam and the reflected beam. Both will be present in an actual system, illustrating one use of the HOE: a means for splitting light off the same surface.

*Diffracted beam Reflected beam*

**FIGURE 20** Holographic optical element (HOE), showing diffracted and reflected beams.

<!-- transcript: s04_hoe -->

NEW prints the element back as the prescription file will hold it (Appendix A.2 lists the complete prescription).

#### Reference Surfaces or Dummy Elements

Reference surfaces are “dummy” elements, with no direct effect on the light. They are used in MACOS to define the start and stop points of diffraction propagations. They are also used to provide a screen on which to project the beam for viewing using ray-trace or diffraction commands. They have many other uses as well.

A Reference may be given `Surface= Zernike` or `Surface= Aspheric` so that it carries a basis definition (modes, coefficients, aperture) for a program to read; the surface terms are stored but never applied to the rays. The one active exception is a Reference with grid data, which is a deliberate phase screen (Section 6). Obscurations on a Reference clip rays but do not touch the diffraction-grid field; a mask for a diffraction propagation must be an `Element= Obscuring` ([Section ?](#setting-up-diffraction-propagations)).

As discussed in Section , reference surfaces can be aligned with the nominally spherical (or flat) wavefront of a particular beam. When done correctly, the reference surface then defines the nominal reference wavefront, or the wavefront the beam would have if it were perfect. The OPD of the actual beam wavefront on the surface ([Section ?](#ray-trace-commands)) then represents the deviation from perfection of the beam wavefront. These differences drive the MACOS diffraction propagators. This approach greatly reduces the sampling density required to avoid undersampling a wavefront for diffraction (ref. Siegman and Szicklas).

MACOS provides functions that help assure that reference surfaces used for diffraction calculations are correctly set up. These are:

-   ORS: Optimizes a reference surface to provide best sampling of a beam ([Section ?](#ors)).

-   SRS: Slaves one reference surface to another to properly set up a near-field diffraction propagation ([Section ?](#srs)).

Figure 21 illustrates a typical scenario, with a reference surface inserted in the focused beam of a lens. The reference surface has its center located at the focus of the lens, and so is aligned with the ideal (spherical) wavefront. The Cassegrain example (see Appendix A.1) provides another example; indeed all of the examples in the Appendix use reference surfaces.

![](media/fig_rescued_21.png)

**FIGURE 21** Reference surface correctly aligned with wavefront.

<!-- transcript: s04_reference -->

NEW prints the element back as the prescription file will hold it.

#### Focal Planes

The focal plane is a reference element, but one that always has a flat surface. The only difference between MACOS FocalPlane and Reference elements is in the way that linear models are computed for them. Use focal planes to end a prescription that is to be used for linear models! As a a reference element, a focal plane element has no direct optical effect. Focal planes are commonly used to terminate far-field diffraction propagations.

<!-- transcript: s04_focalplane -->

NEW prints the element back as the prescription file will hold it.

#### Return Surfaces

The return surface is another element which is used for calculations. It has no effect on wavefront propagation. Specifically, two return surfaces in sequence are used to find the position of the exit pupil of the system in the FEX command. This discussed more thoroughly in [Section ?](#beam-set-up-commands).

<!-- transcript: s04_return -->

NEW prints the element back as the prescription file will hold it.

#### Polarizers and Waveplates

Two element types apply a polarization transformation to the beam without
bending it: `TrPolarizer`, an ideal linear polarizer, and `WavePlate`, a
linear retarder. Both are transmissive and both are active only when
polarization is enabled (the `POLARIZATION` command, or `macos.polarization`
/ `pymacos` equivalents). With polarization off they behave exactly as
`Reference` surfaces, so adding one to a prescription cannot change a
scalar result.

Geometrically they are `Reference` surfaces: the surface may be `Flat` or
`Conic`, the ray is intersected and its path length accumulated in the usual
way, and the ray direction is unchanged. What they add is a 2×2 Jones matrix
applied to the ray's electric field in the plane transverse to the ray.

Two new keywords describe them:

| Keyword | Applies to | Meaning |
|---|---|---|
| `PolAxis=` | `TrPolarizer`, `WavePlate` | Reference axis as a 3-vector in global coordinates: the TRANSMISSION axis of a polarizer, the FAST axis of a waveplate. Required. |
| `Retardance=` | `WavePlate` | Retardance in waves at the prescription `Wavelen`. Required for a waveplate; 0.25 is a quarter-wave plate, 0.5 a half-wave plate. |

`PolAxis=` is a direction, not a position, and it does not have to lie in the
element's surface — MACOS derives the element's material axis from it (see
below) and projects that into each ray's transverse plane, normalizing as it
goes. The projection must not degenerate: an axis that leaves the element with
no transverse direction extinguishes the ray rather than having MACOS choose
an arbitrary basis. Both keywords
are required on the element types that use them: a defaulted axis would
orient a polarizer arbitrarily, and a defaulted retardance would make a
waveplate a silent no-op.

`Retardance=` follows the same convention as `Coating=` layer thickness. The
prescription value is in waves at the `Wavelen` current when the file is
parsed, and MACOS converts it to a physical retardance
(*n*<sub>slow</sub> − *n*<sub>fast</sub>)·*d* at load. The trace then divides
by the wavelength in force at the time. A plate is therefore a fixed piece of
glass and a wavelength sweep is chromatic, as a real plate is: a quarter-wave
plate at 1 µm is an eighth-wave plate at 2 µm.

The retardance sign follows MACOS's propagation convention
(`exp(+iωt)` time dependence, so a field accumulates
`exp(−i·2π·L·N/λ)` through a medium). The slow axis therefore takes the more
negative phase, and the element's Jones matrix in its own (fast, slow)
eigenbasis is diag(1, exp(−iδ)) with δ = 2π·`Retardance`. Fast axis leads.

Both elements are thin idealizations. There is no ray splitting, no
walk-off, no Fresnel loss at the plate faces and no substrate thickness, and
the output field is purely transverse. In particular a polarizer discards the
rejected component rather than emitting it, so a polarizing beamsplitter
needs two traces — one per output port. For a beamsplitter at a substantial
angle of incidence, a coated `Reflector` is usually the better model: there
the s and p directions are set by the physical plane of incidence, and the
thin-film recursion gives real diattenuation and retardance rather than an
idealization.

A waveplate's retardance is likewise independent of the angle of incidence,
where a real crystal plate's is not — the field-of-view effect that drives
compound and Pancharatnam designs. Bounding that needs a birefringent-plate
model with ordinary and extraordinary indices and a thickness, which MACOS
does not provide.

Off normal incidence, what MACOS projects into the transverse plane is the
element's *material* axis — the direction fixed in the substance of the
element. This matters because projection does not preserve orthogonality, so
projecting the declared pass axis is not the same operation as projecting the
block axis and transmitting the complement; the two differ by
acos(2 cos *a*/(1 + cos² *a*)) of transmitted-axis orientation, which is 3.56°
at 20° of incidence and identically zero at normal incidence. For a waveplate
the material axis is the declared (fast) axis. For a polarizer it is the
*absorbing* direction — the wires of a wire grid, the aligned dipole chains of
a dichroic sheet — which is the in-element complement of the declared pass
axis; MACOS forms it, projects it, extinguishes it, and transmits its
orthogonal partner. The keyword still declares the pass axis, so only its
component perpendicular to the element normal has any effect, and a pass axis
parallel to the normal extinguishes the ray. This is the dipole model, which
Korger et al. (*Opt. Express* **21**, 27032, 2013) measured on a tilted
dichroic polarizer.

The reflective polarizer type `RfPolarizer` is reserved but not implemented:
a reflective wire grid carries grid reflection efficiency and the substrate's
own s/p response beyond the axis rule.

    iElt= 3
    EltName= QuarterWavePlate
    Element= WavePlate
    Surface= Flat
    KrElt= -1.000000000D+20 KcElt= 0.000000000D+00
    psiElt= 0.000000000D+00 0.000000000D+00 -1.000000000D+00
    VptElt= 0.000000000D+00 0.000000000D+00 2.000000000D-01
    RptElt= 0.000000000D+00 0.000000000D+00 2.000000000D-01
    PolAxis= 1.000000000D+00 1.000000000D+00 0.000000000D+00
    Retardance= 2.500000000D-01
    IndRef= 1.000000000D+00
    Extinc= 0.000000000D+00
    nCoat= 0
    nObs= 0
    ApType= None
    PropType= Geometric
    zElt= 2.000000000D-01
    nECoord= -6

#### Hex and Pie Segmented Mirrors

Segmented systems (e.g. telescopes with segmented primary mirrors) are defined at the light source using GridType=hex or GridType=pie. This section describes entering the segmented surfaces. Figures 22 and 23 show the geometry for a segmented pie-shaped aperture and a large array of hexagonal segments, respectively. In both cases, the pannel coordinates are given by three variables (X,L,R). The segments must be entered in as elements in the same order they are listed in SegCoord. Segmented surfaces also require

-   nSeg: number of segments

-   gap: spacing between segments

-   width: width of segments (flat to flat)

![](media/fig_rescued_22.png)

**FIGURE 22** Hexagonal coordinates X,L,R for a pie-shaped segmented mirror.

More than 1 segmented mirror can be used in a single prescription. In the case of 2 segmented mirrors, the segments are required to be lined up, so that the segmentation defined at the source is valid for both beams.

An example of a segment mirror .in-file (seg.in) is in Appendix A.5. The panels are numbered in the same order as they are in the Figure 23.

Input data:

<!-- transcript: s04_segment_1 -->
![](media/fig_rescued_23.png){width="5.0in"}

**FIGURE 23** Hexagonal coordinates X,L,R for a large array

<!-- transcript: s04_segment_2 -->

Element data:

<!-- transcript: s04_segment_3 -->

Segments 3 to 6 are entered the same way, each with its own rotation point; NEW takes the vertex, axis and shape from segment 1.

<!-- transcript: s04_segment_4 -->

#### Generating segmented prescriptions: SegMirMaker

Hand-editing a .in file to define every segment of a 7- or 19-segment mirror is tedious. MACOS v4.00 ships with a separate tool, SegMirMaker, in MACOS_resources/segmirmaker/, that generates a .presc (suitable for inclusion in a .in file) and an Hx.m edge-sensor measurement matrix from a parent prescription.

SegMirMaker accepts:

-   A MACOS .in parent prescription. The parent element can be any conic or FreeForm surface.

-   The number of rings (nRing = 1 for 7-segment, 2 for 19-segment, etc.).

-   The axis orientation psi, segment size or aperture diameter, inter-segment gap, and segment standoff.

-   3-DOF (piston/tip/tilt) or 6-DOF segments.

For FreeForm parents, SegMirMaker copies the parent's lFF, FFCoef, FF coordinate frame, grid data, and grid frame verbatim into every segment, so all segments share the parent's underlying FreeForm shape. Each segment's Mon slot is left empty (MonZernCoef = 0) but with lMon = L2 (half the segment width); the user can later add per-segment figure error via the Mon Zernike coefficients without re-running SegMirMaker.

SegMirMaker is the modernized successor to the SMPGe ('Segmented Mirror Prescription Generator') VAX Fortran tool from 1992; the original is preserved in MACOS_resources/segmirmaker/Archive/SMPGe.for. See MACOS_resources/segmirmaker/README.md for full usage.

#### Flower Segmented Mirrors

Flower segmented mirrors consist of a central segment surrounded by “petal” segments. As such, they are like the “pie” segmented mirrors of the previous section, except flower segmentation is not limited to a hexagonal geometry. Flower segmentation is defined at the light source using GridType=Flower. Two other parameters are used to define the segmentation:

-   nPetals: number of petals surrounding the central segment

-   radCtr: radius of the largest inscribed circle defining the size of the center segment

Only 1 segment coordinate is needed per segment, to specify its location in sequence around the central segment, with the 0th segment being the center segment, the 1st segment having one edge along the SegXgrid axis and the other along a line rotated 360°/

nPetals towards the SegYgrid axis, and the 2nd and subsequent segments located in the same sequence around the center segment. See Fig. 24 for an illustration.

![](media/fig_rescued_24.png)

**FIGURE 24**Flower segmentation, showing segment coordinates

#### Reflective Non-Sequential Surfaces

In most optical systems, all rays strike one element before passing on to the next element. There are exceptions, such as corner cube reflectors, where some rays will strike an element before others. These are modeled in MACOS using non-sequential surfaces. There is no limit on the number of non-sequential elements and the order that the non-sequential elements are entered into MACOS does not matter. However, there must be sequential surfaces before and after the group of non-sequential surfaces. For the purpose of the following discussion, we call them the first and second surfaces, respectively. After the ray passes the first sequential surface, MACOS solves for the next surface that the ray strikes which can be any of the group of non-sequential surfaces or the second sequential surface. If the ray hits a non-sequential surface, the procedure continues until the second sequential surface is encountered. After the second sequential surface is encountered the ray trace proceeds normally. Figure 25 shows the corner cube reflector example. The .in-file is in Appendix A.3

NOTE: Geometric propagation only is supported through non-squential elements.

![](media/fig_rescued_25.png)

**FIGURE 25** Non-sequential surface: Ray A has Face2 as first element, Ray B’s first element is Face1.

Corner cube example

<!-- transcript: s04_nonseq -->

#### Refractive Non-Sequential Surfaces

Refractive non-sequential surfaces are very similiar to reflective non-sequential surface (see [Section ?](#reflective-non-sequential-surfaces)). An example illustrating their application is a Luneberg Lens Antenna is provided in Appendix A.3. The Luneberg Lens is a stepped-gradient index lens consisting of concentric spherical shells, which grow denser as they decrease in radius. A point source on the surface produces rays that follow nearly elliptical trajectories, emerging on the other side of the lens in a collimated beam. Using multiple point sources distributed on the surface, the Luneberg Lens can provide instantaneous, solid-state beam coverage over nearly 4π steradians.

### Element Surface Types

#### Flat Surfaces

Flat surfaces are completely defined by a point (the “vertex” VptElt) on the surface, and the unit vector normal to the surface at that point (psiElt). Technically, the *vertex* of a plane or sphere is undefined; the user is free to pick any point on the surface as the vertex. In MACOS, a flat surface is a spherical surface with an infinite radius of curvature, KrElt=1d22. See [Section ?](#reflectors) for an example of entering flat surface data.

#### Conic Surfaces

Conic surfaces – really conicoids of revolution – include the most common optical shapes: spheres, ellipsoids, paraboloids, hyperboloids. Flat surfaces are spheroids with an infinite radius of curvature. MACOS specifies the shape of a conic surface using two surface-specific parameters: the conic constant KcElt and the “radius,” KrElt. Table 7 lists the common names associated with special values of the conic constant.

A hyperboloid or an ellipsoid has two sheets. A ray is intersected with the sheet the vertex lies on (the two are told apart by the sign of the surface normal's axial component), and the far sheet is never taken even when a ray crosses it first; on a Return or Reference the root nearest the element's reference point RptElt is taken, which lets the exit-pupil sphere of FEX be traced back from the focus. An intersection choice that differs from the pre-2026 rule is counted and reported once per run in the WARN summary; such a count on an old prescription means it was tracing onto the wrong sheet.

**TABLE 7.** Conic coefficients

+---------------------+------------------------------------------------+
| **Conic             | **Conic section**                              |
| coefficient**       |                                                |
+=====================+================================================+
| KcElt\<-1           | hyperboloid                                    |
+---------------------+------------------------------------------------+
| KcElt=-1            | paraboloid                                     |
+---------------------+------------------------------------------------+
| -1\<KcElt\<0        | ellipsoid with major axis on the principal     |
|                     | axis (also known as prolate spheroid), psiElt  |
+---------------------+------------------------------------------------+
| KcElt=0             | sphere                                         |
+---------------------+------------------------------------------------+
| KcElt\>0            | oblate ellipsoid                               |
+---------------------+------------------------------------------------+

MACOS also accepts input of conic parameters in terms of the geometric focal length (which is *not* the same as optical focal length of a spherical mirror) and eccentricity which are denoted fElt and eElt, respectively. Eccentricity and focal length are related to the conic constant and radius by:

$$
eElt \;=\; \sqrt{-KcElt}
\qquad (4.1)
$$

$$
fElt \;=\; \frac{-KrElt}{1 + \sqrt{-KcElt}}
\qquad (4.2)
$$

The “sag equation” for a conic surface, which is to say, the equation which defines all points on the surface, is:

$$
Z(h) \;=\; \frac{h^2}{KrElt + \sqrt{(KrElt)^2 - (1 + KcElt)\,h^2}}
\qquad (4.3)
$$

Here h is the radial distance to the projection of the point of incidence onto a plane perpendicular to the principal axis and including the surface vertex; and Z is the height of the surface above that plane (Figure 26).

![](media/fig_rescued_26.png)

**FIGURE 26** Conic mirror parameters

Note that the *vertex* of a conic refers to the geometric vertex except in the case of an oblate ellipsoid. For an oblate ellipsoid (KcElt\>0), the “vertex” is taken as the point where the minor axis intersects the surface.

See [Section ?](#refractors) for an example of entering conic surface data.

#### Generalized (10th-Order) Aspheric Surfaces

Generalized aspheric surfaces use a conic surface as a base, and then impose a deviation on it, defined by additional, higher order, symmetric terms in radius h. In MACOS, generalized aspheric surfaces take the form

$$
Z(h) \;=\; \frac{h^2}{KrElt + \sqrt{(KrElt)^2 - (1 + KcElt)\,h^2}}
 + A h^4 + B h^6 + C h^8 + D h^{10}
\qquad (4.1)
$$

where h is the radial pupil coordinate (the magnitude projection of the of the point of incidence vector onto the plane perpendicular to psiElt as shown in Figure 29), KrElt is the radius of curvature of the element, KcElt is the conic constant, and A, B, C, D are the aspheric coefficients, AsphCoef(i), with AsphCoef(1) = A, AsphCoef(2) = B, AsphCoef(3) = C, and AsphCoef(4) = D.

Note: The A,B,C,D parameters are used by optical design codes as well, except that the signs of the coefficients may differ in some cases. The sign of the coefficients in some codes is a function of the direction of the incident ray. MACOS references surfaces to a global coordinate system and does not assume a particular direction of the beam in defining the optical train. In MACOS the coefficients are positive in the direction of the element psiElt vector.

This example shows how to enter a higher order aspheric reflector.

<!-- transcript: s04_asphere -->

An `AsphCoef=` (or `AnaCoef=`) line that carries fewer values than the surface needs is accepted: the missing higher-order coefficients are zero, and one line naming the element and the count is printed at load.

#### Zernike Surfaces

Zernike surfaces add a deformation—described by a Zernike polynomial—to the surface of a conicoid. The deformation is entered as the “excess sag” of the surface (*not* the wavefront error). MACOS allows up to 66 Zernike terms (37 in the Fringe orderings, 45 in `NormHex`) to be added to a conicoid base surface shape. Table 8 gives the terms through radial order 4 in each ordering MACOS supports. Appendix E gives all of them, with the normalized types, and says how the table is measured from the engine.

The ordering and normalization are selected per element with the `ZernType=` keyword: `ANSI` (the legacy name `Malacara` is accepted), `BornWolf`, `Fringe`, `Noll`, `ExtFringe`, and the RMS-normalized `NormANSI`, `NormBornWolf`, `NormFringe`, `NormNoll`, and `NormHex` (orthonormal over a hexagon inscribed in the lMon circle, vertices on ±xMon; Appendix E, Table E.2). Appendix E describes each ordering. `NormAnnularNoll` is parsed but refused at load until an annular converter exists (it traced as circular `NormNoll` at any ratio), and an `ExtFringe` coefficient in modes 47–49 or 58–64 is refused at load, naming the mode, because those terms (radial orders 11–14) are not converted; Appendix E has the measurements. **If no `ZernType=` is given, the ANSI ordering is assumed** — for every surface that carries Zernike terms (`Zernike`, `ZrnGrData`), whatever the element type. **A prescription converted from Code V must declare `ZernType= Fringe`** (or `Noll`) explicitly; the converter's ordering is not inferred. A `nZernCoef=` line on a surface that carries no Zernike terms is refused at load, naming the element. Fewer than `nZernCoef` coefficients may be entered; the rest are zero. The same keyword family applies to the Mon and FF channels of FreeForm surfaces (`MonZernType`, `FFZernType`; see “FreeForm Surfaces” below).

<!-- transcript: zernike_table8 -->

KcElt and KrElt (or eElt and fElt) are the parameters that describe the conicoid. The deviations occur only in the z-direction of a special coordinate frame. The variable pMon is the center of the coordinate frame. The variables xMon, yMon, and zMon define the special coordinate frame as shown in Figure 27.The variable lMon is the semi-diameter of the circle to which the Zernike deviations are scaled (remember that Zernikes are defined on a unit circle).

![](media/fig_rescued_27.png)

**FIGURE 27** Specialized coordinate frame

Zernike polynomials, listed in Table 8 and Appendix E, are orthogonal over a unit circle. The polar coordinates (R, A) are related to the coordinates of the projection of the point of incidence vector onto the plane perpendicular defined by xMon and yMon (X, Y) by X=R lMon cos(A) and Y=R lMon sin(A). The more common Seidel aberrations can be calculated from the Zernike coefficients (ref. Malacara).An example using a Zernike surface to simulate deformations on the primary mirror of a small telescope is provided in Appendix A.7. The following dialog illustrates the responses the NEW command expects in entering a Zernike reflector:

<!-- transcript: s04_zernike -->

#### Sparse Zernike input

In addition to the dense MonCoef form (120 coefficients in Cartesian monomial order), MACOS reads Zernike coefficients in a sparse 'modes + coefs' form:

    MonZernType= BornWolf
    nMonZernCoef= 4
    MonZernModes= 4 9 16 25
    MonZernCoef= -4.819E+00 4.001E+00 1.197E+00 2.269E-03

MonZernType selects the Zernike ordering and normalization. Supported types: ANSI, BornWolf, Fringe, NormANSI, NormBornWolf, NormFringe, NormNoll, Noll, NormHex, ExtFringe; NormAnnularNoll is parsed and refused at load until an annular converter exists. The Norm\* variants use RMS normalization; the unprefixed variants are unnormalized. The orderings are in Appendix E.

Conversion to the internal monomial form (MonCoef) happens automatically at trace time. The same scheme applies to FF Zernikes: FFZernType, nFFZernCoef, FFZernModes, FFZernCoef. Modes are 1-indexed; a modes line may be written on one line or wrapped over several.

#### Monomial Surfaces

MACOS has a second polynomial-based deformed-surface type, evaluated as a sum of Cartesian monomials. The generating function is:

$$
Z(x,y) \;=\; \sum_{i=0}^{14}\;\sum_{j=0}^{i} c_k\, \bar{x}^{\,j}\, \bar{y}^{\,i-j},
\qquad \bar{x} = \frac{x}{lMon},\; \bar{y} = \frac{y}{lMon},
\qquad (4.1)
$$

where (x, y) are the coordinates of the point of incidence in the special frame (pMon, xMon, yMon, zMon) of Figure 27, the deviation is along zMon, and the 120 coefficients MonCoef(k) are taken in that order: k = 1 is the constant term, then for each total degree i = 1 … 14 the terms j = 0 … i (so MonCoef(2) multiplies ȳ, MonCoef(3) multiplies x̄, MonCoef(4) ȳ², MonCoef(5) x̄ȳ, MonCoef(6) x̄², …).  (The printed 3.2 manual left this equation blank; the form is the engine's MonomialEval.)

The user must specify 120 coefficients in the prescription file.

#### Toric Surfaces

Toric surfaces are surfaces that have different shapes in some specified orthogonal x-and y-axes.Cylindrical and toroidal surfaces are two examples. Toric surfaces are distinguished from the anamorphic surfaces by there being a non-zero conic constant in ony 1 of the axes. MACOS toric surfaces are entered using the standard KrElt and KcElt format, except that two sets of these conic parameters are entered: one set for each of the x-and y-axis directions. The variable AnaCoef holds the values in the following order: Krx, Kry, Kcx. In addition, a specialized coordinate frame must be used to define the local x- and y-axes for these variables. This is the same specialized coordinate frame used in [Section ?](#zernike-surfaces) and is defined by pMon, xMon, yMon, and zMon as shown in Figure 27. An example of an (cylindrical) anamorphic surface as defined in a .in-file is shown below

#### Anamorphic Surfaces

Anamorphic surfaces are surfaces that have different (conic) shapes in some specified orthogonal x- and y-axes. Cylindrical optics and astigmatic plates are two examples. MACOS anamorphic surfaces are entered using the standard KrElt and KcElt format, except that two sets of these conic parameters are entered: one set for each of the x- and y-axis directions. The variable AnaCoef holds the values in the following order: Krx, Kry, Kcx, Kcy. In addition, a specialized coordinate frame must be used to define the local x- and y-axes for these variables. This is the same specialized coordinate frame used in [Section ?](#zernike-surfaces) and is defined by pMon, xMon, yMon, and zMon as shown in Figure 27. An example of an (cylindrical) anamorphic surface as defined in a .in-file is shown below.

    iElt= 1
    EltName= ana_surf Element= Refractor Surface= Anamorph
    psiElt= 0.000000000D+00 0.000000000D+00 -1.000000000D+00 VptElt= 0.000000000D+00 0.000000000D+00 0.000000000D+00 RptElt= 0.000000000D+00 0.000000000D+00 0.000000000D+00 IndRef= 1.000000000D+00
    Extinc= 1.000000000D+22
    AnaCoef= -1.700000000D+02 -1.000000000D+00 -1.000000000D+22 0.000000000D+00 pMon= 0.000000000D+00 0.000000000D+00 0.000000000D+00
    xMon= -1.000000000D+00 0.000000000D+00 0.000000000D+00 yMon= 0.000000000D+00 1.000000000D+00 0.000000000D+00 zMon= 0.000000000D+00 0.000000000D+00 -1.000000000D+00
    nObs= 0
    ApType= None
    zElt= 1.000000000D+22
    PropType= Geometric nECoord= -6

#### User-Defined Surfaces: Deformable Mirrors

MACOS comes with 5 built-in user-defined surfaces. These surfaces use influence functions to model standard deformable mirror configurations for adaptive optics applications. They assume a grid of push-pull actuators spread out in regular patterns across the surface of the element. Coefficients that define the extension of each actuator are read from a file (or are directly modified by a calling program, if you are using SMACOS). The mirror coefficient file is specified in the prescription (or through the MODIFY command) as:

    UDSrfFile=’filename’

The surface types, in order of UDSrfType (as `usersub.F` dispatches them):

1.  1313-actuator deformable mirror.

2.  349-actuator deformable mirror.

3.  397-actuator triangular-array deformable mirror.

4.  349-actuator deformable mirror with squared-out influence function.

5.  349-actuator deformable mirror with ideal actuator influence functions (bilinear approximation).

6.  349-actuator deformable mirror with squared-out influence function.

7.  149-by-149 x-y grid (written for the Hubble primary mirror).

8.  115-by-115 x-y grid (written for the Hubble secondary mirror).

9.  Radial-actuator deformable mirror.

10. SQAna-by-SQAna deformable mirror with damped-sine influence functions.

11. x-y grid surface.

Each type takes its parameters from

    UDSrfParam= p1 p2 ... p6

up to six reals on ONE line, given BEFORE `UDSrfFile=` (the file is read with them; with no `UDSrfParam=` the engine sets the type's defaults at load and says so).  Six is the width everywhere -- the prescription parser, MODIFY, the element-entry prompt and SAVE; values past the sixth are ignored, and a short line is padded with 0 with a one-line note.  The models read at most three:

-   types 1-6 and 10 (influence-function mirrors): p1 = beta, the influence-function shape parameter; p2 = dAct, the actuator spacing; type 10 also p3 = SQAna, the actuators per side;

-   type 9 (radial): p1 = drad, the radial actuator spacing; p2 = dq, the angular spacing (rad); p3 = minrad, the inner radius;

-   types 7, 8 and 11 (x-y grids): p2 = the grid spacing; p1 is not read.

(The element arrays hold ten per element, `mUDSrfParam`; only six are ever read or written.  Measured 2026-10-17, PLAN_CONSOLIDATION 1.2 (h).)

An example of the use of deformable mirrors, as well as lenslet arrays and the ATMOS

atmospheric phase disturbance command, is provided in Appendix A.4.

#### Interpolated Surfaces Using Regularly Gridded Data

MACOS provides interpolated surfaces of 2 types. The first type uses regular square gridded data, in the form of a matrix, whose rows and columns are filled with values giving the surface height deviation from a base conicoid of revolution. This form of interpolated surface, described in this section, is compatible with the high-density information typically generated by optical test interferometers. It uses linear interpolation to determine the intercept point of each ray at the optical surface, rather than the more compute-intensive 3rd-order interpolation of the more general X-Y-Z interpolated surface discussed in the next section.

Grid-data surfaces are defined using 3 additional parameters. These are:

-   mGridMat: the dimension of the gridded data matrix (nGridMat by nGridMat)

-   GridSrfdx: spacing of the data grids in base units

-   GridFile: name of the text file containing the nGridMat by nGridMat surface matrix

The maximum number of grid-data surfaces, and the maximum number of points for grid-data surfaces, are limited at compile time, by setting the mGridSrf and mGridMat parameters in the param.cmn file.

The surface data file GridFile should be in ASCII text, with nGridMat rows and

nGridMat columns. The data is read in using the following Fortran READ sequence:

    ((GridMat(i,j,jGridSrf),i=1,nGridMat),j=1,nGridMat)

The i and j indices correspond to the xMon and yMon coordinates, respectively:. The matrix is assumed centered at the pMon point. Increasing index values represent steps of size GridSrfdx along that axis. The maximum extent of the surface is ±(nGridMat/ 2)\*GridSrfdx in either direction. Rays that fall outside the covered area are passed, but using only the base conic to determine direction and pathlength.

The grid spans (nGridMat − 1)·GridSrfdx base units in each direction, centered on the grid frame origin pData; a ray outside the span sees no figure (this is not an error), so a GridSrfdx too small for the element leaves figure only near its center. A grid declared larger than the compile-time mGridMat is not traced: the file is declined by name at load and the grid term is inert. When an element carries both a grid and Zernike terms, the grid frame xData, yData should equal the Zernike frame xMon, yMon, or a figure applied to one segment will not localize.

The data is read in when the prescription is read. No setup computations are required (no SINT command). Interpolation data can also be provided to SMACOS directly from a calling program that includes the appropriate common statements (see [Section ?](#sharing-data-via-fortran-modules)).

#### Interpolated Surfaces Using X-Y-Z Data

Interpolated surfaces using X-Y-Z data triplets allow arbitrary perturbations to be applied to a base conic-of-revolution surface, using irregularly sampled data. Compared to the gridded-data surface of the previous section, the X-Y-Z interpolated surfaces are slower to compute, using a 3rd-order, triangular, continuous surface interpolation algorithm. They also require large amounts of RAM.

The maximum number of interpolated data surfaces, and the maximum number of points for interpolated data surfaces, are limited at compile time, by setting the mIntSrf and mDP parameters in the param.cmn file. These settings are shown in the command window as MACOS is started up.

The interpolation data is entered from an input file, which is not read until the user executes a SINT (Surface INTerpolator) command. The SINT command reads the data files defining the surface(s) and precomputes interpolation coefficients used in the ray trace.

If rays are traced prior to executing the SINT command, the interpolated surfaces are converted to conic surfaces and the interpolation data is not used.

Interpolation data can also be provided to SMACOS directly from a calling program that includes the appropriate common statements (see [Section ?](#sharing-data-via-fortran-modules)).

The input data file can be in either binary or text data format. A binary surface data file is expected to have the name “filnam.srf*N*.bin,” where filnam is the loaded prescription's file name *including* its `.in` and *N* is the element number (for `Cassegrain.in`, element 3: `Cassegrain.in.srf3.bin`). A text surface data file is expected to have the name “filnam.srf*N*.txt.” If no file with the expected name is found, the user is prompted for a different file name. If none is provided, the surface is not initialized and will be converted to a conic when rays are traced.

The first line of the file contains two integer numbers: nDP, the number of nodes for the optical element; and iComputeDZ, which indicates whether slope data is provided or not. If iComputeDZ=1, slopes are to be estimated from the height data; iComputeDZ=0 indicates that slopes are provided in the data file.

The nDP subsequent lines of the file contain the data for the interpolation. The data for one node is in one row in the file. Each row has 3 or 5 real numbers, giving the x and y coordinates of the data point (in the element's local frame), the surface height at that point — the deviation from the element's base conic, positive along psiElt, toward the incoming beam — and then, if iComputeDZ=0, the x and y slopes at that point. The grid defining the data does not need to be a regular grid. The grid points need not be in any particular order in the data file.

#### FreeForm Surfaces

FreeForm (SrfType = 14, Surface = FreeForm) is a single MACOS surface type combining up to four independent contributions:

    z(rho) = z_conic(rho) + z_Mon(rho) + z_FF(rho) + z_grid(rho)

where:

z_conic is the conic of revolution defined by KrElt, KcElt, and the optical axis psiElt through the vertex VptElt.

z_Mon is a Cartesian-monomial sag in the local frame (pMon, xMon, yMon, zMon), normalized by aperture lMon, with coefficients MonCoef (or sparse MonZernCoef). Active only when lMon \> 0.

z_FF is a free-form Zernike sag in an independent frame (pFF, xFF, yFF, zFF), normalized by lFF, with coefficients FFCoef (or sparse FFZernCoef). Active only when lFF \> 0.

z_grid is an interpolated sag from a regularly-gridded data array, on an independent frame (pData, xData, yData, zData) with grid spacing GridSrfdx. Active only when nGridMat \> 0.

Mon, FF, and grid contributions can be enabled in any combination; setting the corresponding sentinel to zero (or omitting it from the prescription) disables that term.

Why three independent frames. In segmented mirrors generated by SegMirMaker (§4.4.11), the FF coefficients describe the parent surface and the FF frame is parent-aligned (so all segments share the same FF Zernike set), while the Mon frame is per-segment (face-aligned) so each segment can carry its own figure error via MonZernCoef. The grid frame is parent-aligned so the data array represents a global perturbation map.

Independence is enforced in the surface evaluation. The FreeFormSrf routine treats Mon, FF, and grid sags as displacements along their own z-axes (zMon, zFF, zData). Earlier versions assumed all three shared a single z-axis, which produced incorrect intersections when, for example, the segment Mon frame was face-aligned with zMon = -psi while zFF = +psi (a 180-deg mismatch).

Example. Conic + Mon + FF + grid:

    iElt= 1
    EltName= primary
    Element= Reflector
    Surface= FreeForm
    KrElt= -1.000D+02
    KcElt= -1.0
    psiElt= 0 0 1
    VptElt= 0 0 0
    RptElt= 0 0 0
    lMon= 50.0
    pMon= 0 0 0
    xMon= 1 0 0
    yMon= 0 1 0
    zMon= 0 0 1
    MonZernType= BornWolf
    nMonZernCoef= 16
    MonZernModes= 1 3 4 5 9 10 11 12 13 19 20 21 22 23 24 25
    MonZernCoef= 1.197E+0 ...
    lFF= 100.0
    pFF= 0 0 0
    xFF= 1 0 0
    yFF= 0 1 0
    zFF= 0 0 1
    FFZernType= BornWolf
    nFFZernCoef= 16
    FFZernModes= ...
    FFZernCoef= ...
    nGridMat= 256
    GridFile= flat.txt
    GridSrfdx= 31.25
    pData= 0 0 0
    xData= 1 0 0
    yData= 0 1 0
    zData= 0 0 1

### Coatings and Polarization Data

The polarization ray trace (Section 5, *Polarized Ray Tracing*) reads
three kinds of prescription data beyond the scalar trace: the complex
index of each element, an optional thin-film coating stack per element,
and the source polarization state.

**Element index and extinction.** Every element carries `IndRef=` and
`Extinc=`, forming the complex index N = n − iκ (κ > 0 is loss, under
the engine's exp(+iωt) time convention). In the scalar trace only the
real index of refractors matters; with polarization on, a bare
(uncoated) surface applies the Fresnel s/p coefficients computed from
its complex index. The standard mirror idiom

    IndRef=   1.0E+00
    Extinc=   1.0E+22

is a perfect conductor: it reflects with RP = RS = −1 and is
polarization-neutral (no diattenuation, no retardance). To model a real
metal mirror, either give the element the metal's complex index
directly, or (preferred) add a coating stack whose outer layer is an
optically thick metal film — thickness much greater than the skin depth
reproduces the bare metal's Fresnel coefficients regardless of
substrate.

**Thin-film coating stacks (`Coating=`).** An element may carry a
multilayer coating, applied by a recursive thin-film calculation at
each ray's actual angle of incidence when polarization is on. The
keyword gives the layer count, followed by one line per layer, ordered
**outermost to innermost**, each line holding the layer's real index,
extinction coefficient, and thickness:

    IndRef=   1.0E+00
    Extinc=   1.0E+22
    Coating=  2
      1.38  0.0   0.25
      2.30  0.10  0.125

Two rules deserve emphasis:

- **Layer thickness in the prescription is in units of waves** — waves
  of the `Wavelen=` value current at the point the parser reaches the
  `Coating=` line. It is converted to physical thickness
  (thickness·Wavelen/IndRef) at load time, so the loaded stack is
  physically fixed and wavelength sweeps after loading are handled
  correctly. (The programmatic interface, `coating`/`coat_set` in the
  bindings, takes **physical** thickness in BaseUnits directly and has
  no wavelength coupling.)
- **`Coating=` must appear after the element's `IndRef=`/`Extinc=`
  lines** — the parser snapshots the boundary media when it reads the
  coating block.

The incident medium for the coating calculation is the medium the ray
is actually traveling in when it reaches the surface (the running
index, e.g. air), and the substrate is the element's own complex index.

**Source polarization state (`Ex0Ey0=`).** The complex source field
amplitudes can be carried in the prescription:

    Ex0Ey0=   1.0  0.0  0.0  0.0

giving Ex real, Ex imaginary, Ey real, Ey imaginary, in the source grid
frame. This sets the *state* only; turning polarization tracing on
remains a command/API decision (POLarization, or `polarization('on')`
in the bindings, which also prompt for or accept the state directly).
The keyword round-trips through SAVE.

### Obscurations and Apertures

Apertures are the areas through which light can pass. MACOS supports circular and rectangular apertures as shown in Figure 28. When a wavefront is propagated, the intensity that does not go through the hole is set to zero. Likewise, obscurations block light that falls within their interior. Figure 29 shows the types of obscurations implemented in MACOS: circular, rectangular, annular, and triangular. Geometrical rays are unaffected by the obscurations; the rays are traced through the optical system regardless. However, the OBScuration command will allow one to view only the unobscured rays in the spot diagram.

When entering data for obscurations and apertures, fElt, eElt, and, psiElt are required inputs. This may seem unnatural, but including these properties allows curvature to be given to apertures and obscurations. If an aperture or obscuration is planar, then fElt=1d22 and eElt=0. In this case, psiElt is normal to the plane. However, if the obscuration or aperture is not planar, such as an obscuration caused by mirror support pads, fElt, eElt, and, psiElt are the same as they would be for that mirror. All obscuring apertures are defined on the plane defined by psiElt and VptElt and projected onto the surface.

This section has three different examples. The first example is of an element representing a square aperture in front of the primary mirror.

![](media/fig_rescued_28.png)

**FIGURE 28** Examples of apertures types.

![](media/fig_rescued_29.png)

**FIGURE 29** Examples of obscuration types and some composite obscurations.

<!-- transcript: s04_obscuration_1 -->

The second example shows an obscuration consisting of four circles (see Figure 29). This is similar to what was done to model the Hubble Space Telescope’s primary mirror which was drilled in three places for rubber retaining rings. This example is a Cassegrain telescope example with four support pads on the primary mirror. fElt and eElt are the same as the primary mirror. NGridpts was also increased so the obscurations could be seen in the focal plane. Note that the support pads are located just in front of the primary mirror.

<!-- transcript: s04_obscuration_2 -->

The third example shows an obscuration produced by a secondary mirror and it’s support spider: two rectangles and a circle (see Figure 29). These planar obscurations obscure the primary mirror.

<!-- transcript: s04_obscuration_3 -->

#### Obscuration types and the ObsVec= line

Every obscuration of an element is two prescription lines, `ObsType=` and
`ObsVec=`, repeated `nObs=` times; the shapes lie in the element's local
plane, x along `xObs=`, y = psiElt × xObs, origin at the vertex, in base
units.  The `ObsVec=` values by type (the interactive dialog asks for the
same numbers, in the same order):

| ObsType=                      | ObsVec= values                         | Meaning |
|-------------------------------|----------------------------------------|---------|
| Circle / NegCircle            | r, xc, yc                              | radius and center |
| Rectangle / NegRectangle      | xmin, xmax, ymin, ymax                 | axis-aligned box |
| Annulus / NegAnnulus          | rOuter, xc, yc, rInner                 | ring |
| Triangle / NegTriangle        | q1x, q1y, q2x, q2y, q3x, q3y           | three corners |
| Ellipse / NegEllipse          | xrad, yrad, xc, yc                     | axis-aligned ellipse |
| Matrix / NegMatrix            | px, py, dx, dy                         | a 0/1 matrix on a pixel grid; the dialog then asks the matrix dimension and file name |
| Tapered_Ellipse               | a, b, xc, yc, dr                       | ellipse whose transmission tapers over the outer fraction dr of a and b |
| Polygon / NegPolygon          | xc, yc, n — then n lines of x y        | convex polygon, vertices in order around it, relative to (xc, yc); see below |
| RotRectangle / RotNegRectangle | W, H, dx, dy, angle                   | box of width W and height H centered at (dx, dy), rotated by angle (radians) |
| RotEllipse / RotNegEllipse    | a, b, dx, dy, angle                    | ellipse of half-axes a, b centered at (dx, dy), rotated by angle (radians) |

A **Neg** shape is a hole: a ray is obscured when it falls inside a positive
shape and inside no negative one, so a NegCircle inside a Rectangle makes a
plate with a round opening.

**Polygon obscurations in three forms.**  A polygon must be convex, with at
most 128 vertices (`mPolySide`), listed in order around its perimeter.
Besides the 2-D form in the table, the vertices may be given in the GLOBAL
frame as x y z triplets, one per line, after a `PolyObsVec= n` (or
`Poly3DObsVec= n`) line in place of `ObsVec=`:

    ObsType= Polygon
    PolyObsVec= 4
    0.5 0.5 -4.2
    -0.5 0.5 -4.2
    -0.5 -0.5 -4.2
    0.5 -0.5 -4.2

MACOS takes the vertices' centroid as the polygon center and projects each
vertex onto the element's local plane; the projection needs the element's
`xObs=`, and the order of the lines no longer matters — a polygon written
before its element's `xObs=` is held until `xObs=` is read (or until the
load's default for it), where it used to be a hard error.  The 2-D form of
the same square is `ObsVec= 0 0 4` followed by `0.5 0.5`, `-0.5 0.5`,
`-0.5 -0.5`, `0.5 -0.5`.  On the Cassegrain example (4 m aperture, 51-point
grid, 1877 rays) the secondary obscuration — a circle of radius 0.5 and its
spider — loses 121 rays; this 1 m square in its place loses 145 in all three
forms, the ray density of the grid times its area (`cli_tests/cmd/OBS_POLYGON.jou`).
Apertures have the same polygon forms under `ApType= Polygonal` with
`PolyApVec=`, but there `xObs=` must still come first.

### Specifying Coatings

Multi-layer coatings are not supported in this release of MACOS.

### Example Telescope

In [Section ?](#example-light-source) the source data for the Cassegrain telescope example was entered. This section finishes the example adding the necessary optical elements: the primary mirror, secondary mirror, reference surface, and focal plane. Figure 15 shows the telescope. Table 9 lists the optical design.

**TABLE 9.** Cassegraing telescope design

+--------------+-------------+-------------+--------+---------------+
| **Element**  | **Focal     | **Ecc       | **Diam | **Distance    |
|              | length**    | entricity** | eter** | from          |
|              |             |             |        | primary**     |
+==============+=============+=============+========+===============+
| Primary      | 5.4         | 1.0         | 4.0    |               |
| Mirror       |             |             |        |               |
+--------------+-------------+-------------+--------+---------------+
| Secondary    | 1.338854098 | 1.634183595 | 1.0    | 4.061145901   |
| Mirror       |             |             |        |               |
+--------------+-------------+-------------+--------+---------------+
| Focal Plane  | 1d22        | 0.0         |        | 1.5           |
+--------------+-------------+-------------+--------+---------------+

Continuation of MACOS input from [Section ?](#example-light-source).

Element 1, the secondary mirror's obscuration (with the aperture):

<!-- transcript: s04_new_elt1 -->

Element 2, the primary mirror:

<!-- transcript: s04_new_elt2 -->

Element 3, the secondary mirror:

<!-- transcript: s04_new_elt3 -->

Element 4, the reference surface:

<!-- transcript: s04_new_elt4 -->

Element 5, the focal plane:

<!-- transcript: s04_new_elt5 -->

The output coordinate frame:

<!-- transcript: s04_new_tout -->

This completes the input example. The file is now saved as Cassegrain.in within the directory in which MACOS is running. The input file is included in Appendix A.1.

### Editing .In-Files

Once a prescription has been entered and stored as a .in-file, the user can edit and reuse it for further analyses. Editing an .in-file is a matter of opening the file in the text editor

of your choice. On Unix machines this can be vi, Emacs, TextEdit, or another text editor. On the Macintosh, Qued-M, MPW Editor, BBEdit, Alpha, or other text or word processors can be used.

Files can be opened in directories other than the one MACOS is running in, by using appropriate path descriptors in specifying the file name. On Unix systems, the full name takes the form: “subdirectory/subdirectory/filename.” On Macintosh systems, the full name takes the form: “:subdirectory:subdirectory:filename.”

There are some simple rules for .in-files that must be followed to avoid ambiguities in the data. These are:

-   The first data block must be for the source.

-   Each element’s data must be in a block starting with ‘iElt=’ followed by the number of that element. All data within that data block refers to that element. The end of the data block is indicated by the next ‘iElt=’ line or the ‘nOutCord=’ line.

-   Ordering within data blocks is generally not important, except that nObs must preceed the description of the obscurations.

-   The last element data block is followed by the ‘nOutCord=’ specifying the number of output coordinates and the output coordinate matrix (`Tout=`).  A prescription that ends without this block still loads: the loader proposes the beam frame at the last element (z along the beam from the previous element’s vertex to the last element’s vertex, or along the chief ray when there is one element; x from the element’s `xObs=`, or the global axis least aligned with z; y = z × x), prints it, and in the interactive program asks `Keep the proposed output coordinates? [YES]:` — NO reads `nOutCord` and the `Tout` rows at the prompt.  The MATLAB and Python bindings take the proposal and print the same lines.  SAVE writes the block either way.

-   A numeric `PropType=` (the legacy form, 0 through 15 in the order of the propagation-type table) is accepted and resolved to its name at load; `PropType= 0` is Geometric.  SAVE writes the name.  A value that is neither a name nor a number in the table is refused at load, naming it (`** Rx load refused: PropType= <value> (elt n): not a propagation type`).

-   Blank lines can be used to separate data for clarity.

-   Some optional keywords are never written by the NEW dialog and are added by hand (or by SAVE when the value differs from the default). In the source block: `UseChfRay4OPD= Y|N` selects the OPD reference (the chief ray, `Y`, is the default; `N` selects the aperture mean — [Section ?](#ray-trace-commands)) and `OPDRefRayLen=` fixes the reference pathlength directly. In an element block: `ZernType=` (ANSI when absent — see “Zernike Surfaces”) and `NSCount=` (the non-sequential hit budget — see “Reflective Non-Sequential Surfaces”).

-   Comments.  A whole line beginning with ‘%’ is a comment, and SAVE writes it back (up to 500 such lines; an overflow is reported once at SAVE).  A ‘%’ after a value (`Key= value  % note`) ends the value: the note is read past but is NOT written back by SAVE.  Block comments open with `/*` or `CommentBegin` and close with `*/` or `CommentEnd`, each marker first on its own line; SAVE writes the block back as it was.  A one-line `/* note */` is not supported (VALIDATE refuses it as a comment block never closed).

-   The design-optimization keywords (`OptTarget=`, `VarElt=`, `OptBeam*=`, …) are listed in Table 13 at the end of this section.  They are reset at every load; SAVE writes the ones whose values differ from their defaults, and a load that reads any of them prints one line naming the setup in force.

-   An unknown or misspelled keyword in the source block is reported once per load (`MBFile6: Unidentified header string = <key>`) and otherwise ignored; in an element block the message is `MBFile6: Unidentified string`.

-   Data strings can be entered in free format. Up to 16 decimals of precision will be recognized.

Data for a single element must be grouped together, but the order of the data is not important as long as the ‘iElt=’ statement is first. When a new ‘iElt=’ statement is encountered, MACOS checks to make sure all the necessary data has been provided for the previous element. If data for an element has not been provided, MACOS will do one of two things. For some elements, MACOS will provide “default” data for that element. Other data is “necessary” for an .in-file and MACOS will print a warning that the .in-file was not loaded and the name of the omitted variable. A listing of the defaults provided for variables is included in [Section ?](#summary-of-prescription-variables).

### Commands

MACOS provides the user several commands to help setup the optical prescription, to modify or perturb the optical prescription, and to list parameters of the optical system. Commands to help Some of these are described here; others are mentioned here but described elsewhere.

#### OLD and LOAD

Previously created .in-files are loaded using the OLD command. This command has one argument, the name of the input file (minus the “.in” suffix). For instance, to load the Cassegrain.in prescription:

<!-- transcript: s04_old -->

LOAd also work identically:

<!-- transcript: s04_loa -->

#### VALIDATE

The VALidate command syntax-checks a .in prescription file without loading it into MACOS. This is useful when sharing or hand-editing prescription files: rather than discovering a typo halfway through a load, the user can check the file first. Like OLD, the command takes the filename (with or without the .in suffix):

<!-- transcript: s04_validate -->

If the file is well-formed, MACOS prints "Cassegrain.in: OK". If a problem is detected, the message identifies the offending line and key, for example on a copy of Cassegrain.in with its first KcElt= value deleted:

<!-- transcript: s04_validate_bad -->

The same check runs automatically as part of the OLD and LOAD commands; if it fails, the load is aborted and the user is re-prompted for a different filename, leaving the previously loaded prescription untouched.

The validator reads the file as text only.  It refuses a key with no value (`EltName=`, `BaseUnits=` and `WaveUnits=` excepted, since older prescriptions leave them blank), a blank line inside a multi-row block such as a `Tout=` matrix, and a comment block that is never closed.  It accepts `%` and `!` comment lines and the block-comment forms described under “Editing .In-Files”; note that the loader itself treats only `%` as a comment character.  Files with CR or CRLF line endings load through a temporary LF copy and need no conversion.

#### MODIFY

The MODify command allows the user to change the current optical prescription without editing the .in-file. These changes are not saved in a .in-file unless the SAVe command is used. This allows the user to perturb elements, modify the prescription, and save the results if desired. Any of the .in-file variables listed in [Section ?](#summary-of-prescription-variables) can be changed.  Write the assignment without spaces around the equals sign (`mod nGridpts=33`); `mod nGridpts = 33` is rejected with a message rather than read.

<!-- transcript: s04_mod_help -->

To use the MODify command, just type MOD. The user will then be prompted for the modification. If a MODification is successful, it will be echoed at the terminal. For instance to change the eccentricity of element 2 to 1.03.

<!-- transcript: s04_mod_ecc -->

If a MOD is not successful, it will not be echoed.

<!-- transcript: s04_mod_fail -->

There are 2 ways to enter use MODify for vectors and arrays. The first is changing one element of the vector or array at a time:

<!-- transcript: s04_mod_vec1 -->

The second is to re-enter the entire vector:

<!-- transcript: s04_mod_vec2 -->

The SHOw command can be used to display the element’s data.

#### SAVE

The SAVe command writes a new .in-file containing the current optical prescription data to disk. Parameters changed as the result of a PERturb or MODify command will be written as modified. This command takes one argument: the new .in-file name (without the “.in” suffix). If the specified file already exists, the user is asked if it should be overwritten.

What SAVE writes is the optical system as loaded and modified: every element block, the whole-line and block comments of the original file (in-line `% notes` are dropped), the stop (`ApStop=`), and the optional keywords whose values differ from their defaults (`UseChfRay4OPD= N` when the aperture-mean OPD reference is in force, `ZernType=`, `NSCount=`, `OPDRefRayLen=`).  A ray trace never changes what SAVE writes.  SAVE also writes the design-optimization keywords of Table 13 whose values differ from their load defaults (the target, the variable elements and their bounds, the field list and weights, the beam targets, the iteration cap and tolerance), so a saved prescription carries the CALIB setup made on it, and the next load announces that setup (`LOAD: CALIB setup in force from the prescription: …`).  SAVE does not write the session options set by command (POLARIZATION, VECTOR, NOLNEG, CENTROID, CIR), which belong to the session rather than the system.  Loading a saved file and saving it again reproduces the file: the test suite checks this round trip on every prescription in the repository.

#### SUMMARIZE

The SUMmarize command displays a brief summary of the loaded optical system, listing each element and showing its type, conic constant KcElt, and radius of curvature KrElt.

#### STATUS

The STAtus command displays information about some of the current states.

<!-- transcript: s04_status -->

#### RESET

The RESet command reinitializes some of the current states (see below) to their default values.

<!-- transcript: s04_reset_status -->

#### SHOW

The SHOw command displays all the prescription data for a particular element. It has one argument, the element number.

### Summary of Prescription Variables

Tables 10 through 12 provide a summary of the MACOS prescription variables and their default values. When MACOS loads an .in-file, these defaults are used variables, if the user has not supplied values. When a default is used, the user is notified, e.g.:

    Default used for zSource

Not all variables have defaults, and these variables must be specified by the user. If one of these variables is not given a value in the .in-file, then the .in-file is not loaded and the user is notified; for instance:

    VptElt(1) must be specified by user Input file not properly loaded

The variables are shown with their index types in parentheses. These indices are:

-   iObs: which obscuration in an element (can be up to 6 obscurations per element).

-   iElt: number of the element in the .in-file.

-   iAxis: ranges over the x,y,z coordinates in the global coordinate system.

-   iPert: ranges over x,y,z rotations, x,y,z translations in the u-vector.

-   iLin: ranges over the x,y,z components of the ray direction perturbation, the x,y,z components of the ray transverse aberration and the pathlength perturbation in the linear state vector x.

-   iAsCoef: A, B, C, and D coefficients described in [Section ?](#generalized-10th-order-aspheric-surfaces).

-   iAp: either radius,x,y, or xmin,xmax,ymin,ymax depending on ApType.

-   iZern: the 36 Zernike coefficients.

-   iLayer: ranges over the coating layers.

All variables can be changed with the MODify command.

**TABLE 10.** Light source variables

+---------------------+-----------------+-----------------------------+
| **Source            | **Defaults**    | **Description**             |
| Variables**         |                 |                             |
+=====================+=================+=============================+
| Aperture            | 1.0D0           | input aperture              |
+---------------------+-----------------+-----------------------------+
| BeamType            | 1 (uniform)     | beam intensity              |
+---------------------+-----------------+-----------------------------+
| ChfRayDir(iAxis)    | NO DEFAULT      | chief ray direction         |
+---------------------+-----------------+-----------------------------+
| ChfRayPos(iAxis)    | (0.0            | chief ray location          |
|                     | D0,0.0D0,0.0D0) |                             |
+---------------------+-----------------+-----------------------------+
| CosPower            |                 |                             |
+---------------------+-----------------+-----------------------------+
| ExtInc(0)           | 0.0D0           | extinction coefficient      |
+---------------------+-----------------+-----------------------------+
| Ex0                 |                 | electric field strength of  |
|                     |                 | source along x-axis         |
+---------------------+-----------------+-----------------------------+
| Ey0                 |                 | electric field strength of  |
|                     |                 | source along x-axis         |
+---------------------+-----------------+-----------------------------+
| Flux                | 1.0D0           | total intensity received    |
|                     |                 | across the unobscured       |
|                     |                 | aperture                    |
+---------------------+-----------------+-----------------------------+
| gap                 | 0.0D0           | spacing between segments in |
|                     |                 | seg-mented apertures        |
+---------------------+-----------------+-----------------------------+
| GridType            | 1 (circular)    | type of input aperture      |
+---------------------+-----------------+-----------------------------+
| IndRef(0)           | 1.0D0           | real part of index of       |
|                     |                 | refraction of input medium  |
+---------------------+-----------------+-----------------------------+
| nElt                | NO DEFAULT      | number of elements          |
+---------------------+-----------------+-----------------------------+
| nGridPts            | IFIX(mpts/2)    | number of ray grid points   |
|                     |                 | across the aperture         |
+---------------------+-----------------+-----------------------------+
| nSeg                | 0               | number of segments          |
+---------------------+-----------------+-----------------------------+
| Obscratn            | 0.0D0           | input obscuration           |
+---------------------+-----------------+-----------------------------+
| rxBeam              |                 | beam diameter in            |
|                     |                 | x-direction                 |
+---------------------+-----------------+-----------------------------+
| ryBeam              |                 | beam diameter in            |
|                     |                 | y-direction                 |
+---------------------+-----------------+-----------------------------+
| SegCo-ord(xyz,iElt) | NO DEFAULT      | segment location in the     |
|                     |                 | entrance pupil              |
+---------------------+-----------------+-----------------------------+
| Wavelen             | 6.0D-7          | wavelength of light         |
+---------------------+-----------------+-----------------------------+
| width               | Aperture        | width of segments           |
+---------------------+-----------------+-----------------------------+
| xGrid(iAxis)        | NO DEFAULT      | orientation of input x-axis |
|                     |                 | vector                      |
+---------------------+-----------------+-----------------------------+
| yGrid(iAxis)        | cross product   | orientation of input y-axis |
|                     | of xGrid and    | vector                      |
|                     | ChRay-Dir       |                             |
+---------------------+-----------------+-----------------------------+
| zGrid(iAxis)        |                 | orientation of input z-axis |
|                     |                 | vector                      |
+---------------------+-----------------+-----------------------------+
| Zelt(0)             |                 |                             |
+---------------------+-----------------+-----------------------------+
| zSource             | 1.0D22          | distance of light source    |
|                     |                 | from point of origin        |
+---------------------+-----------------+-----------------------------+

**TABLE 11.** Element variables

+---------------------------+-----------------+-----------------------+
| **Element variables**     | **Default       | **Description**       |
|                           | values**        |                       |
+===========================+=================+=======================+
| AnaCoef(iAna,iElt)        |                 | anamorphic            |
|                           |                 | coefficients          |
+---------------------------+-----------------+-----------------------+
| ApType                    | 1 (circle)      | aperture type         |
+---------------------------+-----------------+-----------------------+
| ApVec(iAp,iElt)           | NO DEFAULT      | aperture definition   |
|                           |                 | informa-tion          |
+---------------------------+-----------------+-----------------------+
| AsphCoef(iAsCoef,iElt)    | (0.0            | asphere coefficients  |
|                           | D0,0.0D0,0.0D0, | (4)                   |
|                           | 0.0D0)          |                       |
+---------------------------+-----------------+-----------------------+
| CoatIndxElt(iLayer, iElt) |                 | index of coating      |
|                           |                 | layers                |
+---------------------------+-----------------+-----------------------+
| CoatThkElt(iLayer,iElt)   |                 | thickness of coating  |
|                           |                 | layers                |
+---------------------------+-----------------+-----------------------+
| eElt(iElt)                | 0.0D0           | element eccentricity  |
+---------------------------+-----------------+-----------------------+
| EltName(iElt)             | ' ' (blanks)    | element name          |
+---------------------------+-----------------+-----------------------+
| Element(iElt)             | NO DEFAULT      | element type          |
+---------------------------+-----------------+-----------------------+
| ExtInc(iElt)              | 0.0D0           | extinction            |
|                           |                 | coefficient           |
+---------------------------+-----------------+-----------------------+
| fElt(iElt)                | 1.0D22          | element focal length  |
+---------------------------+-----------------+-----------------------+
| h1HOE(iAxis,iElt)         |                 |                       |
+---------------------------+-----------------+-----------------------+
| h2HOE(iAxis,iElt          |                 |                       |
+---------------------------+-----------------+-----------------------+
| IndRef(iElt)              | 1.0D0           | real part of          |
|                           |                 | refactive index for   |
|                           |                 | element               |
+---------------------------+-----------------+-----------------------+
| KcElt(iElt)               | 0.0D0           | element conic         |
|                           |                 | constant              |
+---------------------------+-----------------+-----------------------+
| KrElt(iElt)               | -1.0D22         | element conic radius  |
+---------------------------+-----------------+-----------------------+
| LensArrayType(iElt)       | 1               | type of lens in array |
|                           |                 | (hex, square)         |
+---------------------------+-----------------+-----------------------+
| LensArrayWidth(iElt)      | 1               | width of a single     |
|                           |                 | lens in the array     |
+---------------------------+-----------------+-----------------------+
| lMon(iElt)                | Aperture        | scale length for      |
|                           |                 | monomial (Zernike)    |
|                           |                 | element calcula-tions |
+---------------------------+-----------------+-----------------------+
| MonCoef(iMon,iElt)        |                 | monomial coefficients |
+---------------------------+-----------------+-----------------------+
| nCoatElt(iElt)            |                 | number of coatings    |
+---------------------------+-----------------+-----------------------+
| nECoord(iElt)             | -6              | number of element     |
|                           |                 | coordi-nates          |
+---------------------------+-----------------+-----------------------+
| nLayer                    |                 | number of layers in   |
|                           |                 | coating, maximum 20   |
+---------------------------+-----------------+-----------------------+
| nObs(iElt)                | 0               | number of             |
|                           |                 | obscurations          |
+---------------------------+-----------------+-----------------------+
| ObsType(iObs,iElt)        | 1 (circle)      | obscuration type      |
+---------------------------+-----------------+-----------------------+
| ObsVec(iAp,iObs,iElt)     | (0.0            | obscuration           |
|                           | D0,0.0D0,0.0D0, | definition            |
|                           |                 | infor-mation          |
|                           | 0.0             |                       |
|                           | D0,0.0D0,0.0D0) |                       |
+---------------------------+-----------------+-----------------------+
| OrderHOE(iElt)            |                 | diffraction order of  |
|                           |                 | HOE to trace          |
+---------------------------+-----------------+-----------------------+
| PinHole()                 |                 | pin hole              |
+---------------------------+-----------------+-----------------------+
| pMon(iAxis,iElt)          | VptElt          | reference for         |
|                           |                 | monomial (Zernike)    |
|                           |                 | surface calcula-tions |
+---------------------------+-----------------+-----------------------+

**TABLE 11.** Element variables

+---------------------------+-----------------+-----------------------+
| **Element variables**     | **Default       | **Description**       |
|                           | values**        |                       |
+===========================+=================+=======================+
| PropType(iElt)            | 1 (Rays Only)   | type of propagation   |
|                           |                 | from element iElt to  |
|                           |                 | element iElt+1        |
+---------------------------+-----------------+-----------------------+
| psiElt(iAxis,iElt)        | NO DEFAULT      | element principal     |
|                           |                 | axis vec-tor          |
+---------------------------+-----------------+-----------------------+
| RptElt(iAxis,iElt)        | VptElt          | element rotation      |
|                           |                 | point                 |
+---------------------------+-----------------+-----------------------+
| RuleWidth(iElt)           |                 | grating ruling width  |
+---------------------------+-----------------+-----------------------+
| Surface(iElt)             |                 | surface type          |
+---------------------------+-----------------+-----------------------+
| TElt(iPert,iElt-Cord,iElt)| 6x6 Identity    | element local         |
|                           | Matrix          | coordinates           |
|                           |                 | transformation matrix |
|                           |                 |                       |
|                           | (0.0            |                       |
|                           | D0,0.0D0,0.0D0, |                       |
|                           |                 |                       |
|                           | 1.0             |                       |
|                           | D0,0.0D0,0.0D0) |                       |
|                           |                 |                       |
|                           | (0.0            |                       |
|                           | D0,0.0D0,0.0D0, |                       |
|                           |                 |                       |
|                           | 1.0             |                       |
|                           | D0,0.0D0,0.0D0) |                       |
|                           |                 |                       |
|                           | (0.0            |                       |
|                           | D0,0.0D0,0.0D0, |                       |
|                           |                 |                       |
|                           | 0.0             |                       |
|                           | D0,1.0D0,0.0D0) |                       |
|                           |                 |                       |
|                           | (0.0            |                       |
|                           | D0,0.0D0,0.0D0, |                       |
|                           |                 |                       |
|                           | 0.0             |                       |
|                           | D0,0.0D0,1.0D0) |                       |
+---------------------------+-----------------+-----------------------+
| UDSrfCoef(iCoef,iElt)     |                 | user defined surface  |
|                           |                 | coeffi-cients         |
+---------------------------+-----------------+-----------------------+
| UDSrfType(iElt)           |                 | user defined surface  |
|                           |                 | type                  |
+---------------------------+-----------------+-----------------------+
| VptElt(iAxis,iElt)        | NO DEFAULT      | element vertex point  |
+---------------------------+-----------------+-----------------------+
| xMon(iAxis,iElt)          | xGrid           | x-axis for monomial   |
|                           |                 | (Zernike) surface     |
|                           |                 | calcula-tions         |
+---------------------------+-----------------+-----------------------+
| xObs(iAxis,iElt)          | NO DEFAULT      | orientation of x-axis |
|                           |                 | for each element      |
+---------------------------+-----------------+-----------------------+
| yMon(iAxis,iElt)          | cross product   | y-axis for monomial   |
|                           | of              | (Zernike) surface     |
|                           |                 | calcula-tions         |
|                           | psiElt and xMon |                       |
+---------------------------+-----------------+-----------------------+
| WaveHOE(iElt)             |                 | recording wavelength  |
|                           |                 | for HOE               |
+---------------------------+-----------------+-----------------------+
| zElt(iElt)                | fElt            | Fresnel propagation   |
|                           |                 | distance              |
+---------------------------+-----------------+-----------------------+

**TABLE 11.** Element variables

+---------------------------+-----------------+-----------------------+
| **Element variables**     | **Default       | **Description**       |
|                           | values**        |                       |
+===========================+=================+=======================+
| ZernCoef(iZern,iElt)      | 45-element      | Zernike coefficients  |
|                           | vector of zeros |                       |
+---------------------------+-----------------+-----------------------+
| zMon(iAxis,iElt)          | psiElt          | z-axis for monomial   |
|                           |                 | (Zernike) surface     |
|                           |                 | calcula-tions         |
+---------------------------+-----------------+-----------------------+

**TABLE 12.** Output variables

+----------------------+--------------+--------------------------------+
| **Output Variables** | **Defaults** | **Description**                |
+======================+==============+================================+
| Tout(iLin,iPert)     | beam frame   | the output coordinate          |
|                      | at the last  | transformation matrix          |
|                      | element      | (proposed when absent; the     |
|                      |              | interactive program asks)      |
+----------------------+--------------+--------------------------------+
| nOutCord             | 5            | number of output coordinates   |
+----------------------+--------------+--------------------------------+

### Design Optimization (CALIB) Keywords

The CALib command runs the native multi-field design optimizer: a Levenberg–Marquardt least-squares solve over the degrees of freedom of the variable elements, scored at up to 12 fields of view and 6 wavelengths (Appendix C lists the command; the Command Reference carries its dialog and the related commands SETC, AVAR, MVAR, DVAR, VARS, AFOV, DFOV, SFOV, FOVS, WLENS).  Its setup is written into the prescription with the keywords of Table 13.  Header keywords belong in the source block, before the first `iElt=`; element keywords belong in the block of the element they vary.  Every keyword in the table is prescription state: it returns to its default at every OLD or LOAD.  SAVE writes the keywords whose values differ from their defaults (header keys in the source block, element keys in their element’s block), so the optimization set up on a prescription — interactively or by the bindings’ `calib_set_*` calls — transfers through a saved file; a load that reads such a setup says so in one line (`LOAD: CALIB setup in force from the prescription: target SPOT at element 2, 1 field(s), 1 wavelength(s), 1 variable element(s), OptMaxItrs 10 -- VARS/FOVS/WLENS list them; SAVE writes them`).  The same settings can be made interactively (SETC, AVAR, AFOV) and from the MATLAB and Python bindings (`calib_set_target`, `calib_set_var_elt`, `calib_set_iter`, `calib_set_tol`, `calib_set_beam`, `calib_set_beam_pos_fov`, `calib_set_beam_wt`; Command Reference Part II).

**TABLE 13.** Design-optimization (CALIB) keywords.  All are reset at load; SAVE writes those whose values differ from their defaults, and a load that reads them announces the setup.

| Keyword | Block | Value | Default | Meaning |
|---|---|---|---|---|
| `OptTarget=` | source | `WFE`, `SPOT`, `ZWF`, `Beam` | `WFE` | The merit.  `OPL` is no longer a target; a prescription naming it is refused at load. |
| `OptTgtElt=` (`OptWFElt=`) | source | element number | nElt − 1 | The element at which the merit is evaluated. |
| `OptTgtWF=` | source | file name; a name ending in `T` is read as text, any other as binary | none | A target wavefront to match instead of a null. |
| `OptSpotSize=` | source | size, base units | 0 | SPOT target: the spot size aimed for (a negative value reads as 0). |
| `OptZernModes=` | source | n m1 … mn | none | ZWF target: the Zernike modes scored (Noll ordering, 1-based; 11 = primary spherical). |
| `OptRayGrid=` | source | points across the grid | nGridpts/2 | The optimizer's own ray grid, same convention as `nGridpts=`; capped at the model size. |
| `ChfRayDir=` / `ChfRayPos=` (`OptChfRayDir=` / `OptChfRayPos=`) | source, repeated | one 3-vector each per field | the source's own | The field list: the first pair is the nominal field; up to 12 pairs. |
| `OptFOVWt=` | source, after the field pairs | w1 … wn | 1 each | Field weights (a row is divided by √w); FOVS lists them.  Before 2026-10-17 a prescription with no `OptFOVWt=` loaded with every weight 0, read by the solver as 1e-10: merits printed then are 1e-10 × today's for the same solve. |
| `Wavelen=`, then `OptWaveLen=` (repeated) or `ArrWaveLen=` | source | wavelengths | the source's | The optimization wavelengths, up to 6; the first is the source wavelength. |
| `OptMaxItrs=` (`OptMxItrs=`) | source | n | 10 | Iteration cap. |
| `OptTol=` | source | > 0 | 1e-12 | Convergence tolerance: the absolute change of the rms merit between accepted steps, in the target's units.  A value ≤ 0 is ignored with a note. |
| `OptFEX=` | source | `Y` / `N` | `N` | Evaluate the WFE at the exit pupil found by FEX at each field; needs a Return or Reference at element nElt − 1. |
| `OptAlg=` | source | `NonLin` / `Lin` | `NonLin` | Levenberg–Marquardt, or the linear (SVD) solve. |
| `OptSvdSvCut=` | source | cut | machine epsilon | Singular-value cut of the linear solve. |
| `OptNomSens=` | source | `Y` | off | Run the nominal sensitivity at load. |
| `OptSavePinv=`, `OptUseSavedPinv=`, `OptPinvFile=` | source | `Y`; `Y`; file name | off; off; `dady0_saved.txt` | Save, or reuse, the linear solve's pseudo-inverse. |
| `VarSrc=` | source | six flags | none | The source as a variable (TIP TILT CLOCK DX DY PIST); moving the source needs a STOP. |
| `OptSrcRpt=` | source | x y z | none | Rotation point of a variable source. |
| `VarElt=` | element | eight flags | none | The element's degrees of freedom: TIP TILT CLOCK DX DY PIST ROC CONIC, nonzero = vary.  A linked (double-pass) element cannot be variable; its master moves it. |
| `VarCons=` | element, after `VarElt=` | lower upper, per varied flag in order | none | Bounds on the varied degrees of freedom; any `VarCons=` switches the solve to the bound-constrained SLSQP driver. |
| `OptZern=` (`VarZern=`) | element | n t1 … tn | none | Surface Zernike terms as degrees of freedom; marks the element variable if `VarElt=` is absent. |
| `OptAsph=` (`OptAsphCoef=`) | element, after `VarElt=` | n t1 … tn | none | Aspheric terms as degrees of freedom.  Above its `VarElt=` it is ignored with a message. |
| `OptBeamDir=` | element | dx dy dz (unitized) | none | Chief-ray direction target at this element. |
| `OptBeamPos=` | element | x y z | none | Chief-ray position target at this element; also makes this element `OptTgtElt`. |
| `OptBeamPosFov=` | element, repeated | x y z per field, in field order | none | Per-field position targets (up to 12). |
| `OptBeamSize=` | element | size | none | Beam size target at this element. |
| `OptBeamRefRayDir=` | element | dx dy dz | none | Reference-ray direction target at this element. |
| `OptBeamWt=` | element | w | 1 | Weight of the beam rows against the WFE or SPOT rows (≤ 0 reads as 1). |
| `OptBeamCentroid=` | element | `Y` / `N` | `N` | The position row reads the centroid of the rays that pass the train to this element, not the chief ray. |

Conventions that go with the table:

-   **Rows and units.**  Every (wavelength, field) block of the least-squares problem holds the target's rows followed by any beam rows.  WFE rows are in the target's units; SPOT and beam rows are in base units.  Beam rows ride on any target, so a WFE solve can hold a pupil position or a chief-ray direction at the same time; weight them with `OptBeamWt=`.

-   **Stopping.**  The solve stops *converged* when two consecutive accepted steps each change the rms merit by no more than `OptTol` and the sum of squares by no more than one part in 10⁶; at the *round-off floor* after ten rejected trials in a row; at the *cap* of `OptMaxItrs`; *failed* when an evaluation fails (the chief ray lost) with no accepted step, or *stalled* when no step is accepted and none fails — the start is already a stationary point of the merit.  A failed evaluation inside a trial counts as a rejected trial.  The CLI prints `Optimization iterations =`, `Final LM merit (chisq) =` and `Stopped: …`; the bindings return the same three through `calib_stats`.  The stopping rule does not depend on an earlier CALIB in the same session.

-   **Aspheric degrees of freedom.**  The differential step of a zero aspheric coefficient (an asphere started from a conic) is sized from the sag at the element's circular aperture (`ApType= 1`).  Declare that aperture; without it the step is a fixed small number that is round-off on a prescription in meters, and one line says so.

-   **Fields.**  Repeat the chief-ray pairs in the source block, one pair per field, the nominal field first, and give `OptFOVWt=` after them (it reads as many weights as there are fields).  The FOVS command lists the result.

-   **Targets at the exit pupil.**  `OptFEX= Y` evaluates the WFE at the exit pupil FEX finds for each field (the chief-ray axis convention of Section 5).  The default is off because FEX is applied at element nElt − 1, which must be a Return or Reference surface for that to work.
