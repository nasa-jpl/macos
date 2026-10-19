## SECTION 8: Differential Ray-Tracing and Linear Optical Models

The MACOS differential ray-trace functions are used to compute optical sensitivity matrices. These provide linear models of the effects of optical element perturbations on optical performance. The linear models can be exercised internally to MACOS using functions such as LPErturb, LOPd and LINtensity which parallel commands introduced earlier. Linear optical models can also be EXPported for use external to MACOS.

    BUild DMBuild PArtials LPErturb LPRead LREset LSPot LOPd
    LPIxillate LINtensity EXPort

### Background

*linear models* of optical systems provide a linear matrix transformation between element perturbations and ray perturbations. These optical models can be used in place of the full-blown ray-trace calculations for systems that undergo only small changes. They can also be used in ways that the nonlinear models cannot, such as in linear statistical studies.

Linear models are generally valid where

-   The small-angle approximation is valid for changes in ray direction

-   The ray transverse aberration is “small” relative to the curvature of the optical elements

The limits are not absolute, but vary from system to system. The user can check the range of validity of a linear model for a particular system. This is done by using the MACOS LPErtrurb, LOPd, LSPot and LINtensity commands to exercise the linear model. Results are compared with exact results produced for the same cases using the MACOS PERturb, OPD, SPOt and INtensity commands. If, over the expected range of perturbation values, the linear models produce results that are within the accuracy requirements of the particular analysis, then they can be safely used for that analysis. Otherwise a nonlinear SMACOS model should be used.

#### Linear Ray States

Rays are defined at any particular point in a beam train by their location, direction and length. In MACOS, the current *exact ray state* is defined as the position, direction and path length of a ray immediately after the current element (see Figure 65). The *position* of the iRayth ray is a 3-vector giving the point of incidence with the current element. MACOS stores it in the array RayPos as RayPos(1:3,iRay). The *direction* of the iRayth ray is a unit-magnitude 3-vector stored in RayDir(1:3,iRay). Note that if the current element is a mirror or lens, the direction is of the reflected or refracted ray. The *optical path length* of the iRayth ray from the source to the current element is a scalar and is stored in CumRayL(iRay).

![](media/fig_rescued_65.png)

**FIGURE 65** Ray state

The effect of a *perturbation* in the optical system, such as a rotation or translation of an optical element, is to *change the ray state* at every downstream element. Such perturbations can be applied using the PERturb command, and the ray retraced to find the new ray state.

If a perturbation is sufficiently small, it will effect the rays *linearly*. That is, the new ray state at every downstream element is very nearly a linear function of the perturbation. This is true for some range of perturbation magnitudes that varies with each system, but is in general many times the wavelength of light.

MACOS linear modeling features use *analytic differential* ray-trace equations to compute the partial derivatives of the rays with respect to element perturbations. The differential ray-trace equations are derived in References 1, 8, and 9. They provide a set of analytical formulas giving the ray partial derivative matrices in a form that is easily computed.

The *linear ray state* differs slightly from the exact ray state. The linear ray state is a single 7-vector, composed of a 3-vector of *ray direction perturbation*, a 3-vector of *ray transverse aberration*, and a scalar for *ray optical path difference*. It is written *x*

$$
\hat{x} \;=\; \begin{bmatrix} d\hat{r} \\ d\hat{\gamma} \\ dL \end{bmatrix}
\qquad (8.1)
$$

where $d\hat{r}$ is the ray direction perturbation. The ray direction perturbation is a 3-vector whose magnitude is small compared to 1 and which is perpendicular to the ray direction, RayDir. The total ray direction is equal to the vector sum of RayDir and $d\hat{r}$ (to the accuracy of the small angle approximation). Figure 66 shows the linear ray state.

![](media/fig_rescued_66.png)

**FIGURE 66** Linear ray state

The transverse aberration γ is the lateral displacement of the ray. It is a vector perpendicular to the ray direction, from the nominal ray position to a point on the perturbed ray. It is not necessarily located on the surface of the element and is not, in general, the derivative of the ray position, RayPos. Exceptions are the focal-plane element type and any surface that is perpendicular to the nominal ray direction.

The path length perturbation, dL, is the difference between the optical path length of the perturbed and nominal rays. It is not the OPD, as that is the difference of the ray and the chief ray. However, the OPD is easily computed from dL as discussed below.

#### Linear Models

The linear ray state at a particular element is a function of perturbations of the upstream elements and of perturbations in the input ray. These perturbations are with respect to a “nominal state” (*xn* )*nom* , which may or may not be zero. For each ray at an element n, the linear ray state has the functional form:

$$
\hat{x}_n \;=\; (\hat{x}_n)_{nom}
 + \frac{\partial \hat{x}_n}{\partial \hat{u}_1}\,\hat{u}_1 + \cdots
 + \frac{\partial \hat{x}_n}{\partial \hat{u}_n}\,\hat{u}_n
 + \frac{\partial \hat{x}_n}{\partial \hat{x}_0}\,\hat{x}_0
\qquad (8.1)
$$

Here the *u* vectors are 6-vectors of element perturbations, composed of 3-vectors θ for element rotational perturbation and δ for element translation perturbation:

$$
\hat{u} \;=\; \begin{bmatrix} \hat{\theta} \\ \hat{\delta} \end{bmatrix}
\qquad (8.2)
$$

The individual optical sensitivity matrices can be printed using the PARtials command. The sensitivity of the nth ray to perturbations of the input (source) ray is a 7 by 7 matrix of the form:

$$
\frac{\partial \hat{x}_n}{\partial \hat{x}_0} \;=\;
\begin{bmatrix}
\dfrac{\partial \hat{r}_n}{\partial \hat{r}_0} & \dfrac{\partial \hat{r}_n}{\partial \hat{\gamma}_0} & \dfrac{\partial \hat{r}_n}{\partial L_0} \\[2ex]
\dfrac{\partial \hat{\gamma}_n}{\partial \hat{r}_0} & \dfrac{\partial \hat{\gamma}_n}{\partial \hat{\gamma}_0} & \dfrac{\partial \hat{\gamma}_n}{\partial L_0} \\[2ex]
\dfrac{\partial L_n}{\partial \hat{r}_0} & \dfrac{\partial L_n}{\partial \hat{\gamma}_0} & \dfrac{\partial L_n}{\partial L_0}
\end{bmatrix}
\qquad (8.3)
$$

The sensitivity of the nth ray to perturbations of the ith element is a 7 by 6 matrix of the form:

$$
\frac{\partial \hat{x}_n}{\partial \hat{u}_i} \;=\;
\begin{bmatrix}
\dfrac{\partial \hat{r}_n}{\partial \hat{\theta}_i} & \dfrac{\partial \hat{r}_n}{\partial \hat{\delta}_i} \\[2ex]
\dfrac{\partial \hat{\gamma}_n}{\partial \hat{\theta}_i} & \dfrac{\partial \hat{\gamma}_n}{\partial \hat{\delta}_i} \\[2ex]
\dfrac{\partial L_n}{\partial \hat{\theta}_i} & \dfrac{\partial L_n}{\partial \hat{\delta}_i}
\end{bmatrix}
\qquad (8.4)
$$

The complete expression for the linear ray state for a single ray in terms of element and input-ray perturbations is written in matrix form as:

$$
\hat{x} \;=\; (\hat{x}_n)_{nom} +
\begin{bmatrix}
\dfrac{\partial \hat{x}_n}{\partial \hat{x}_0} & \dfrac{\partial \hat{x}_n}{\partial \hat{u}_1} & \cdots & \dfrac{\partial \hat{x}_n}{\partial \hat{u}_n}
\end{bmatrix}
\begin{bmatrix} \hat{x}_0 \\ \hat{u}_1 \\ \cdots \\ \hat{u}_n \end{bmatrix}
\qquad (8.5)
$$

The same form is extended to give the linear ray state for multiple rays at the same element

$$
\begin{bmatrix} \hat{x}_{ray1} \\ \cdots \\ \hat{x}_{raym} \end{bmatrix}
\;=\;
\begin{bmatrix} \hat{x}_{ray1} \\ \cdots \\ \hat{x}_{raym} \end{bmatrix}_{nom}
+
\begin{bmatrix}
\left[\dfrac{\partial \hat{x}_n}{\partial \hat{x}_0} \;\; \dfrac{\partial \hat{x}_n}{\partial \hat{u}_1} \;\cdots\; \dfrac{\partial \hat{x}_n}{\partial \hat{u}_n}\right]_{ray1} \\[2ex]
\cdots \\[1ex]
\left[\dfrac{\partial \hat{x}_n}{\partial \hat{x}_0} \;\; \dfrac{\partial \hat{x}_n}{\partial \hat{u}_1} \;\cdots\; \dfrac{\partial \hat{x}_n}{\partial \hat{u}_n}\right]_{raym}
\end{bmatrix}
\begin{bmatrix} \hat{x}_0 \\ \hat{u}_1 \\ \cdots \\ \hat{u}_n \end{bmatrix}
\qquad (8.6)
$$

Both of these equations have the form *x* = *xnom* + *Cu* , where the *C-matrix* is the *optical sensitivity matrix*, giving the change in ray state due to perturbations of the optical elements. The C-matrix provides a linear model of the optical system. The MACOS BUild command assembles these sensitivities from the individual ray partial derivatives.

It can be extremely useful to incorporate local coordinates directly into the optical sensitivity C-matrices. If the perturbations for one element are to be applied in coordinates other than the global coordinate system–perhaps in an actuator coordinate system–then that local coordinate transformation matrix can be included directly within the C-matrix. The transform from element iElt coordinates to global coordinates is written $\text{ }^{G}T^{\,iElt}$.

Similarly, local coordinates can be used at the output end of the linear model. This allows for the creation of models for ray direction only, or transverse aberration or path length only. Focal-plane axes are easily incorporated. We write the transform for ray state from global coordinates to output coordinates as $\text{ }^{out}T^{G}$ . As an example, the 1x7 output transform matrix that captures only the path length ray state in global coordinates is

$$
\text{ }^{pathlength}T^{G} \;=\; \begin{bmatrix} 0 & 0 & 0 & 0 & 0 & 0 & 1 \end{bmatrix}
\qquad (8.7)
$$

As another example, the 2x7 output transform which takes the transverse aberration in ray state and casts it into detector coordinates *x* and *y* is

$$
\text{ }^{detector}T^{G} \;=\;
\begin{bmatrix} 0 & 0 & 0 & \hat{x} & 0 \\ 0 & 0 & 0 & \hat{y} & 0 \end{bmatrix}
\qquad (8.8)
$$

Both actuator and detector coordinates are specified by the user as part of the .in-file prescription.

With the input and output coordinate matrices included, the C-matrix takes the form:

$$
\begin{bmatrix} \hat{x}_{ray1} \\ \cdots \\ \hat{x}_{raym} \end{bmatrix}
\;=\;
\begin{bmatrix}
\text{ }^{out}T^{G}\left[\dfrac{\partial \hat{x}_n}{\partial \hat{x}_0}\,\text{ }^{G}T^{0} \;\; \dfrac{\partial \hat{x}_n}{\partial \hat{u}_1}\,\text{ }^{G}T^{1} \;\cdots\; \dfrac{\partial \hat{x}_n}{\partial \hat{u}_n}\,\text{ }^{G}T^{n}\right]_{ray1} \\[2ex]
\cdots \\[1ex]
\text{ }^{out}T^{G}\left[\dfrac{\partial \hat{x}_n}{\partial \hat{x}_0}\,\text{ }^{G}T^{0} \;\; \dfrac{\partial \hat{x}_n}{\partial \hat{u}_1}\,\text{ }^{G}T^{1} \;\cdots\; \dfrac{\partial \hat{x}_n}{\partial \hat{u}_n}\,\text{ }^{G}T^{n}\right]_{raynRay}
\end{bmatrix}
\begin{bmatrix} \hat{x}_0 \\ \hat{u}_1 \\ \cdots \\ \hat{u}_n \end{bmatrix}
\qquad (8.9)
$$

### Commands for Building Models

#### BUILD

The BUild command traces all of the rays in a beam from the source to the specified final element. It then computes the full C-matrix linear model, giving the sensitivity of each ray at the last element to perturbations of each element in the beam train. It has one argument, the number iElt of the last element in the propagation.

BUild performs the MACOS differential ray-trace calculations, and must be performed before exercising commands that manipulate linear models, such as LPE, LOP, LSP, LIN, or before EXPorting C-matrix information.

Note that the BUild command ignores the presence of reference surfaces and obscuring surfaces unless they are the last surface. The BUild command currently does not work if there are return surfaces in the beam train.

As an example, consider the Cassegrain.in (see Appendix A.1.2). This prescription has output coordinates set up to compute a wavefront model as discussed previously. Once the prescription is loaded using the OLD command, the C-matrix is easily built:

<!-- transcript: s08_build -->

#### DMBUILD

The DMBuild command is not used in this version.

### Commands for Changing Models

#### LPERTURB

The LPErturb command allows the user to enter an element perturbation for use with the linear model analysis commands. LPErturb takes different arguments depending on whether element input coordinates are to be used. The first argument iElt is the number of the element to be perturbed. Note that element 0 (the source) cannot be perturbed using LPErturb. The second argument is a yes-or-no response indicating whether input coordinates are to be used. If NO is selected, the next two prompts are for rotational and translational 3-vector perturbations. If YES is selected, a vector of the same length as the number of input coordinates is input in the input coordinate system.

After each LPErturb command, the current perturbation vector (in global coordinates) is printed out. The LPErturb command *adds* the perturbation data to the existing perturbation vector. Multiple LPErturb commands sum each perturbation with the previous perturbation vector. The perturbation vector can be reset to zero using the LREset command. Using Cassegrain.in:

<!-- transcript: s08_lperturb -->

#### LPREAD

The LPRead command reads element perturbations from a specified datafile. MACOS assumes that the name of the perturbation file is infile.pert. If this file does not exist, MACOS prompts for the perturbation file name.

The perturbation file must be in the following all numeric format. The perturbations must be in global coordinates. The order of the perturbation vectors in the .pert-file must be identical to the order of the elements (the first perturbation vector corresponds to the first element, the sixth perturbation vector corresponds with the sixth element). Each element must have a perturbation vector. If an element is unperturbed, then a string of zeroes should be entered as its perturbation vector. The perturbation vectors have 6 components: a 3-vector rotational perturbations and a 3-vector translational perturbations. The 6 components must be separated by spaces.

<!-- transcript: s08_lpread -->

#### LRESET

The LREset command zeros any previously entered linear perturbations. It should be used, when necessary, to reset the linear perturbation vector after LPErturb commands.

### Commands for Exercisizing Models

#### PARTIALS

The PARtials command prints individual optical sensitivity matrices for specified rays. The user is asked whether to print in input/output or global coordinates. For Cassegrain.in, the partials of the chief ray in global coordinates are:

<!-- transcript: s08_partials_chief -->

The specific partials take the form in Eqns. 8.3 and 8.4. The first matrix shows that source chief ray direction perturbations are magnified 3.036 times at the reference surface. The effect on ray direction at the reference surface of primary mirror tilts and decenters are magnified by 6.606 and 0.5671 times, respectively.

The path length partials show a factor-of-two sensitivity for this ray. Looking at the OPD partials only for a marginal ray, these numbers change:

<!-- transcript: s08_partials_ray2 -->

#### LOPD

The LOPd command generates an OPD map of the beam at the final surface of the linear model. It parallels the OPD command discussed in [Section ?](#ray-trace-commands). Using Cassegrain.in:

<!-- transcript: s08_lopd -->

#### LSPOT

The LSPot command generates a spot diagram of the beam at the final surface of the linear model. It parallels the SPOt command discussed in [Section ?](#ray-trace-commands). Using Cassegrain.in:

<!-- transcript: s08_lspot -->

#### LINTENSITY

The LINtensity command computes beam wavefront phase at the final surface of the linear model and then performs a single far-field propagation to the next element in the prescription. If the last element in the linear model is the last element in the prescription, the user is prompted for transmission distance and other information. LINtensity performs calculations that parallel the propagation commands discussed in Sections 6.3 and 7.2.1, respectively. Using Cassegrain.in:

<!-- transcript: s08_lint -->

#### LPIXILATE

The LPIxilate command computes beam wavefront phase at the final surface of the linear model and then performs a single far-field propagation to the next element in the prescription. If the last element in the linear model is the last element in the prescription, the user is prompted for transmission distance and other information. Finaly, the image is pixilated. The LPIxilate command parallels the PIxilate command discussed [Section ?](#pixilate). Using Cassegrain.in:

<!-- transcript: s08_lpix -->

### Commands for Exporting Models

The EXPort command provides a way to export the C-Matrix of the linear model and is discussed in [Section ?](#data-export).
