# Curvature–Transport Correspondence (CTC) – Research Archive

This repository hosts preprints and related manuscripts developed under the Curvature–Transport Correspondence (CTC) framework.

---

## 1. The Curvature–Transport Correspondence (CTC): From Quantum Effective Mass to Cosmological Dark Matter

**Author:** Ridwan Sakidja (2025)

**Download PDF:**  
👉 [Sakidja_Manuscript_CTC_Preprint.pdf](Sakidja_Manuscript_CTC_Preprint.pdf)

**DOI:**  
👉 https://doi.org/10.5281/zenodo.17805651

**Summary:**  
The Curvature–Transport Correspondence (CTC) proposes that field curvature across physical systems — from quantum effective mass to galactic rotation curves and cosmological density inference — arises from the divergence of underlying transport fluxes. This manuscript presents the core formulation, mathematical structure, and multi-scale implications.


## 2. The Curvature–Transport Correspondence (CTC) and Its Implications for Cosmic Structure Formation

**Author:** Ridwan Sakidja (2025)

**Download PDF:**

👉 [Sakidja_Manuscript_CTC_COSMO_Preprint.pdf](Sakidja_Manuscript_CTC_COSMO_Preprint.pdf)

**DOI:**  
👉 https://doi.org/10.5281/zenodo.17844548

**Summary:**
This manuscript extends the Curvature–Transport Correspondence to cosmic structure formation. In this framework, curvature arises from the divergence of a gravitational transport flux rather than from material density. Cosmic filaments emerge from spatial variations in divergence, black hole interiors avoid singularities through flux saturation, and the Higgs–Planck hierarchy reflects a transition in curvature-response stiffness. Together, these phenomena appear as manifestations of a single principle: curvature is the geometric response of a medium with finite transport capacity.

## 3. The Big Start: Cosmogenesis from a Finite Planck-Phase Boundary

**Author:** Ridwan Sakidja (2025)

**Download PDF:**

👉 [Sakidja_Manuscript_Big_Start.pdf](Sakidja_Manuscript_Big_Start.pdf)

**DOI:**  
👉 https://doi.org/10.5281/zenodo.17873646

**Summary:**
The standard Big Bang picture assumes that spacetime already exists at the moment of origin and is driven to a singular state of infinite curvature and infinite density. The Curvature Transport Correspondence (CTC) offers a different beginning in which the vacuum is a physical medium with stiffness that controls whether curvature or fields can exist at all. This leads to the Big Start, a finite radius Planck phase vacuum state where gravitational transport is saturated and no geometric degrees of freedom are present. This state is a pre geometric and zero entropy vacuum, directly realizing the insight of Sir Roger Penrose that the Universe must begin in an exceptionally ordered condition. As expansion relaxes vacuum stiffness, the layers of physics appear in sequence: curvature mobility at rG, transverse and quantum modes at rT, Higgs condensation and the appearance of mass at rH, and finally a classical FRW spacetime at rS. Cosmogenesis and gravitational collapse follow the same divergence law in opposite directions, and the Big Start replaces the classical singularity with a natural vacuum phase transition in which spacetime, gravity, quantum behavior, and mass appear only when the vacuum becomes soft enough to support them

## 4. Nonlocal Torque Coupling in Disk Galaxies: A Continuum Framework for Interpreting Rotation Curves

**Authors:** Ridwan Sakidja, Armitha Dutta and Justus Lau (2025)

**Download PDF:**  

👉 [SDL_PREPRINT_GALAXY_DYNAMICS_v1.pdf](SDL_PREPRINT_GALAXY_DYNAMICS_v1.pdf)

**DOI:**  
👉 https://doi.org/10.5281/zenodo.17959245

**Summary:**
Observed galactic rotation curves rise, bend, and remain flat at radii where a purely baryonic Newtonian potential predicts declining velocities. The conventional interpretation introduces non-luminous dark matter to supply the missing rotational support. Here we propose an alternative dynamical architecture based on rotational coherence transmission. In this framework the stellar disk is treated not as a set of mechanically isolated Keplerian annuli but as a continuum of energy-coupled shells capable of inheriting, transmitting, and eventually saturating rotational coherence from inner regions.
From this premise we derive three scale-free, falsifiable observables, namely the Curl Index, the Fractional Uplift, and the Surplus Index, quantify torque inheritance, rotational energy surplus, and radial coherence export. These indices can be computed directly from observed rotation curves without invoking additional mass. Applying the framework to the full SPARC database, we find that galaxies cluster into three dynamical regimes predicted by the model. The structured scaling of outer rotation curves is incompatible with strictly local baryonic dynamics but follows naturally from nonlocal torque transport. This particular work concerns only the dynamical inference drawn from galactic rotation curves. The framework is therefore offered as a galactic-scale dynamical alternative to the enclosed-mass interpretation.

## 5. How Symmetry Survives Disorder: An Expectation-Value Perspective

**Authors:** Ridwan Sakidja (2025)

**DOI:**  
👉 https://zenodo.org/records/18057598

**Summary:**
Symmetry is often treated as an exact microscopic property, even in systems with disorder and fluctuations. In realistic physical systems, however, individual configurations are generically distorted by noise or randomness. Here we examine how symmetry survives disorder by distinguishing between microscopic symmetry and symmetry realized in expectation values. Using numerical stress tests based on paraxial wave propagation, we study one-, two-, and three-dimensional disordered media with statistically symmetric but locally asymmetric refractive-index fluctuations. Individual realizations show strong, order-unity symmetry violation. In contrast, symmetry is robustly recovered in ensemble-averaged intensities, even in the presence of strong disorder. This behavior persists in genuinely three-dimensional random media with disorder varying along the propagation direction, demonstrating that symmetry recovery is not an artifact of reduced dimensionality. These results show that symmetry in disordered wave systems is a statistically stabilized property of observables rather than a microscopic constraint.

**Contents**

**Main Results Demonstrate**
1. Individual realizations exhibit order-unity symmetry violation
2. Ensemble-averaged observables recover symmetry robustly
3. Increasing dimensionality increases statistical stability
4. Disorder remains strong and active throughout propagation
5. Symmetry is a statistical property of observables, not a microscopic guarantee

The repository also contains four Jupyter notebooks, corresponding to increasing dimensionality and diagnostic depth:

**1. STRESS_TEST_1D.ipynb**

One-dimensional symmetry stress test
> Paraxial wave propagation with a single transverse coordinate,
> Statistically symmetric, correlated refractive-index disorder,
> No parity enforcement on disorder, grid, or propagation.

Computes:
> Microstate-level symmetry deviation,
> Distribution of symmetry deviation across realizations,
> Symmetry of ensemble-averaged intensity.

Purpose:
Establishes that symmetry is highly fragile at the microscopic level in 1D, but recovers only after ensemble averaging.

**2. STRESS_TEST_2D.ipynb**

Two-dimensional transverse disorder
> Full 2D transverse wave propagation,
> Same statistically symmetric disorder construction,
> No imposed spatial symmetrization.

Computes:
> Single-shot asymmetric intensity patterns,
> Ensemble-averaged intensity ⟨I(x,y)⟩,
> π-rotated comparisons,
> Symmetry deviation metrics,
> Spatial variance map Var[I(x,y)].

Purpose:
Demonstrates partial self-averaging in two dimensions and shows that symmetry recovery is spatially genuine, not a numerical artifact.

**3. STRESS_TEST_3D.ipynb**

Three-dimensional disordered medium
> Disorder varies in x, y, and along the propagation direction z,
> Longitudinal disorder correlations explicitly controlled.

Computes:
> Symmetry deviation as a function of propagation distance,
> Microstate versus ensemble behavior,
> Robustness of expectation-value symmetry in 3D.

Purpose:
Confirms that expectation-value symmetry survives in genuinely three-dimensional random media and is not an artifact of reduced dimensionality.

**4. STRESS_TEST_3D_B.ipynb**

Extended 3D diagnostics and parameter sweeps
> Additional longitudinal correlation-length sweeps,
> Ensemble statistics (mean and variance across realizations),
> Parametric stress tests beyond the main-text figures.

Purpose:
Provides supplementary diagnostics supporting the main conclusions and figures in the Supplementary Material.


## 6. The Bullet Cluster Revisited: Transport, Shocks, and Non-Equilibrium Structure

**Authors:** Ridwan Sakidja (2025)

**Download PDF:**  
👉 [The Bullet Cluster Revisited.pdf](The Bullet Cluster Revisited.pdf)

**DOI:**  
👉 https://doi.org/10.5281/zenodo.18092526

**Summary:**
The Bullet Cluster (1E 0657–558) is widely cited as decisive evidence for particle dark matter because weak-lensing mass peaks are offset from the X-ray–emitting gas and instead align with the collisionless galaxy component. This interpretation implicitly assumes that gravitational response must track local mass density even during a violent, dissipative merger. Here we present an alternative non-equilibrium interpretation within the Curvature–Transport Correspondence (CTC), in which curvature responds to the divergence of transport flux rather than to density alone. Using a simple and fully reproducible cartoon simulation, we show that merger-driven shocks naturally decouple density from coherent transport: the gas becomes dense and stalled, while collisionless components preserve ballistic motion. A nonlocal CTC lensing proxy constructed from ∇^2 Ψ=∇.J remains aligned with the collisionless component, reproducing the observed lensing–gas offset without invoking additional dark matter particles. This result clarifies the physical assumptions underlying the standard interpretation of the Bullet Cluster and highlights the role of non-equilibrium transport during cluster mergers.




