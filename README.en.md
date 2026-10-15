# Racing Molecule · 分子赛车

**[简体中文](README.md) · English**

Design a molecule and race it across 20 nm of gold.

Racing Molecule is an open-source game inspired by nanocar experiments. Build a molecular racer in the workshop, then watch it move, turn, and deform under thermal fluctuations, lattice corrugation, and an external driving field. Race trajectories are computed with classical molecular dynamics, so molecular structure affects performance through force-field interactions.

**[Play online](https://chc-cyber.github.io/racing-molecule/) · [Download releases](https://github.com/chc-cyber/racing-molecule/releases) · [Build and release guide](RELEASING.md)**

The project prioritizes playability while keeping its MD implementation open to inspection. The current model uses an effective driving field, a frozen gold surface, and an interfacial potential adapted for the game. Race results are not predictions of experimental nanocar speeds.

## Getting started

1. Enter the workshop and draw in 2D with Ketcher, edit in 3D with JSmol, or import a MOL file, a single-record SDF file, or a workshop project.
2. Choose a standard race or custom conditions and send the molecule to the track. The game checks the structure and jointly relaxes the molecule against the surface.
3. Start the race once preparation is complete. You can pause and resume, switch between ball-and-stick and space-filling views, and select 1× or 2× playback.
4. The race finishes when the center of mass has advanced a net 20 nm along the track. Results report average velocity in m/s, and the raw trajectory can be downloaded.

**Browser players do not need to install Ketcher, JSmol, LAMMPS, or Python.** The editors, computation module, and potential tables load with the website. MD runs in the browser on the player's device, with no simulation server required. The initial download is substantial; loading time depends on the network and device. Offline operation is not currently guaranteed.

Supported elements are H, B, C, N, O, F, Si, P, S, Cl, Br, and I. Race molecules may contain up to 256 atoms, including hydrogens, with a maximum span of 45 Å. Each chemical environment must still pass valence, force-field parameter, and geometry-relaxation checks.

## Default race settings

| Setting | Current default |
| --- | --- |
| Track | 20 nm × 20 nm |
| Heat-bath temperature | 150 K |
| Driving strength | 0.6, an effective internal game parameter |
| Gold surface | Au(111) proxy potential v2 |
| Lattice corrugation factor | 1.5 |
| MD timestep | 1 fs |
| Langevin damping time | 350 fs |
| Trajectory sampling interval | 50 fs |

Standard races use fixed rules. Custom races allow the exposed parameters to be adjusted before preparation.

Race timing starts from the prepared initial state. Completion requires a **net displacement** of 20 nm along the track relative to the initial center of mass; lateral detours do not count toward that distance. Results use simulation time, independently of browser frame rate or playback speed. The velocity conversion is:

$$
\bar v\,[\mathrm{m/s}]=1000\,\frac{\Delta x\,[\mathrm{nm}]}{\Delta t\,[\mathrm{ps}]}.
$$

At 1× playback, 100 ps of trajectory takes about 50 seconds to display; at 2×, about 25 seconds. Computation and playback run separately. Playback waits for buffered frames when necessary; changing playback speed does not change the timestep or race result.

## Computational model

### From simulation to animation

LAMMPS handles time integration and the heat bath. RDKit handles molecular processing and UFF force evaluation. The project's own force routines handle the surface, driving field, and boundaries. The forces are combined at each MD step and passed to LAMMPS; Three.js plays the resulting coordinates.

| Component | Implementation and role |
| --- | --- |
| Intramolecular interactions | RDKit UFF: bond stretching, angle bending, torsion, inversion, and nonbonded van der Waals interactions |
| C/H–Au interactions | A frozen periodic Au(111) surface based on SAIP, accelerated by proxy potential v2 by default |
| Other supported elements–Au | An isotropic 12–6 van der Waals approximation using mixed UFF parameters |
| Driving | A static nonuniform electric field acting on fixed partial charges and an effective polarization response |
| Boundaries | A conservative soft-wall potential limiting escape from the race area |
| Time evolution | LAMMPS velocity-Verlet integration with a Langevin heat bath |
| Rendering | Three.js displays actual MD coordinates; the gold atoms are a visual representation of the surface |

The browser version compiles LAMMPS, RDKit, and the project's force routines into WebAssembly and runs them in a Web Worker. The native version calls LAMMPS/RDKit from Python and uses compiled force kernels for acceleration. Both follow the same model conventions; see the [browser engine](wasm/src/engine.cpp) and [native engine](molecular_car/engine.py).

The conservative potential energy is:

$$
U=U_{\mathrm{UFF}}+U_{\mathrm{surface}}+U_{\mathrm{drive}}+U_{\mathrm{walls}},
\qquad \mathbf F_i^{\mathrm{cons}}=-\nabla_i U.
$$

LAMMPS uses `real` units: Å for length, fs for time, and kcal/mol for energy. Surface routines that work internally in eV convert their results when passing them to the engine.

### UFF, preparation, and molecular topology

UFF, the Universal Force Field, describes classical mechanical responses across a broad range of elements. The project uses [RDKit's UFF implementation](https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html) with fixed molecular connectivity. Bonds do not break or form, and chemical reactions do not occur during a race.

Preparation includes adding hydrogens, checking parameters, processing 3D conformations, and optimizing geometry. A 2D input can be embedded in 3D with ETKDGv3 and then relaxed with UFF; a 3D input is processed starting from its supplied conformation. Once the molecule is placed on the track, the driving field and heat bath are disabled for further relaxation of the combined intramolecular, surface, and boundary potential. The final norm of the force vector over all atoms must be no greater than **0.0001 kcal/mol/Å** before the race can begin. See [molecular processing](molecular_car/molecule.py) and [joint relaxation](molecular_car/relaxation.py).

This procedure seeks a usable local minimum. It does not guarantee the global minimum or constitute full thermal equilibration at the target temperature. Races currently ignore isotope differences and use standard elemental masses.

Gasteiger partial charges are assigned during preparation and used for the external-field response. The UFF internal energy does not receive an additional intramolecular Coulomb term based on these charges. Charges are not updated self-consistently, and charge transfer is not modeled.

### Ensemble and Langevin heat bath

The racer follows driven Langevin molecular dynamics, summarized by:

$$
m_i\ddot{\mathbf r}_i=-\nabla_i U-\frac{m_i}{\tau}\dot{\mathbf r}_i+\mathbf R_i(t).
$$

Here, $\tau$ is the damping time, and $\mathbf R_i$ is a random force consistent with the temperature and damping. Damping dissipates kinetic energy, while random forces introduce thermal fluctuations, allowing translation, rotation, and deformation alongside the directed motion. The Langevin bath approximates environmental dissipation and fluctuations; gold-substrate phonons are not simulated explicitly.

The implementation combines `fix nve` with `fix langevin`: the former provides velocity-Verlet integration, while the latter applies damping and random forces. **An integration command named `nve` does not mean the actual race is in the NVE ensemble.** With fixed particle number and box dimensions, no driving field, and sufficient equilibration, the Langevin bath targets the canonical NVT distribution. An actual race is externally driven and observed over a finite time, so it is a nonequilibrium process. See the LAMMPS documentation for [integration](https://docs.lammps.org/fix_nve.html) and [Langevin dynamics](https://docs.lammps.org/fix_langevin.html).

The heat bath acts on atomic velocities without forcing the total random force to vanish at every step, so the center of mass can also undergo random motion. The displayed internal temperature excludes center-of-mass translation. The setting of 150 K is a heat-bath target, not a requirement that every frame have exactly that temperature.

### Frozen Au(111) and SAIP

The reference surface contains six Au(111) layers. Their heights come from a gold slab previously relaxed with an Au EAM potential, while the in-plane lattice retains ideal ABC stacking. See the [surface reference data](molecular_car/data/au111_reference.json). Gold atoms remain fixed during a race; their motion is not integrated, and Au–Au internal energy is not repeatedly evaluated.

C/H–Au interactions are based on the parameters and functional form of SAIP, the semi-anisotropic interfacial potential. It includes dispersion attraction and short-range repulsion, with part of the directional dependence described by molecular local normals. An atom can therefore experience different surface interactions at different lattice positions or fragment orientations. The surface provides both vertical adsorption and lateral lattice corrugation.

The literature basis is a potential for interfaces between gold and graphitic systems: [Ouyang, Hod, and Guerra, 2021](https://doi.org/10.1021/acs.jctc.1c00622). LAMMPS documents the [SAIP/metal equations and parameters](https://docs.lammps.org/pair_saip_metal.html). Applying this model to arbitrary player-designed organic molecules is an extension made by this project; the original paper did not validate all such structures.

To accommodate freely designed molecules, the project adapts local-normal handling explicitly:

- Fixed chemical-bond connectivity selects normal-defining neighbors, preventing unrelated atoms brought close by folding from being mistaken for those neighbors.
- Where a local plane can be defined, its normal is computed from neighboring geometry and updated as the molecule moves. Forces arising from changes in that normal are propagated to the relevant neighbors.
- Terminal, four-coordinate, and designated linear environments that do not suit this definition use a fixed surface-axis direction. Near-degenerate geometries receive regularization and evaluation switching.

The surface extends periodically in-plane, with a smoothly tapered 16 Å cutoff for pair interactions. **Only surface query coordinates are made periodic; the racer remains unwrapped.** A long molecule spanning multiple unit cells retains its full coordinates for bonds, normals, center of mass, and external-field evaluation. Periodic image racers are not introduced.

For B, N, O, F, Si, P, S, Cl, Br, and I, the project does not claim equally applicable SAIP parameters. Instead, it geometrically mixes elemental and Au UFF van der Waals parameters to obtain isotropic 12–6 interactions, then sums and caches them over the frozen periodic lattice. These elements can participate in the game, but chemisorption, bonding to gold, and charge transfer are outside the model. See the [periodic surface implementation](molecular_car/periodic_au.py).

### How the racer is driven

Races use a nonuniform field fixed in space. Let $x$ be the forward track direction and $z$ the height above the surface. The potential and field are:

$$
\phi(x,z)=-(E_0+gx)z,\qquad
\mathbf E=-\nabla\phi=(gz,\,0,\,E_0+gx).
$$

The external-field energy is:

$$
U_{\mathrm{drive}}=\sum_i\left[q_i\phi(\mathbf r_i)-\frac12G\alpha_i|\mathbf E(\mathbf r_i)|^2\right].
$$

Here, $q_i$ is the fixed partial charge, $\alpha_i$ is an effective polarization-response coefficient assigned by element, and $G=500$ is a game response gain. For driving strength $s$ and field-growth parameter $R$, the implementation uses $E_0=0.5s/\sqrt R$ and $g=0.0015s\sqrt R$, with defaults $s=0.6$ and $R=1$. These quantities use internal units consistent with the energy expression. The displayed value of 0.6 must not be interpreted directly as 0.6 V/nm.

The field is constant in time and does not depend on a moving tip. Partial-charge responses and field-dependent polarization energy jointly affect motion; molecular shape, surface contact, and thermal fluctuations further determine the trajectory. There is no direct transverse electric-field component, but surface forces and random forces remain. See the [field implementation](molecular_car/fields.py) for equations and unit conventions.

## How proxy potential v2 works

Proxy potential v2 is a **numerical approximation and accelerated evaluation of the frozen periodic gold-surface potential energy**. It precomputes surface interactions for different positions and local orientations, then evaluates them quickly during a race. Racer coordinates, deformation, and dynamics are still computed step by step. It requires no machine-learning training and does not use prerecorded race trajectories.

### 1. Generate tables from the same frozen surface

The reference calculation sums Au lattice contributions within the cutoff. Repeating this has a cost even when the gold atoms are fixed, especially as the molecule's local normals change.

The implementation separates cacheable scalar baseline terms from orientation-dependent terms. Position-dependent baseline terms use periodic grids; the orientation-dependent carbon contribution uses the Fourier expansion and directional interpolation described below. Ordinary hydrogen interactions with fixed normals can be included in the cache. Other supported elements retain their own radial van der Waals caches. Thus, “proxy potential v2” does not mean that every element uses a single shared Fourier table.

### 2. Represent in-plane periodicity with a Fourier expansion

The orientation-dependent energy uses the in-plane periodicity of Au(111) in a finite expansion:

$$
U_{\mathrm{dir}}(x,y,z,\mathbf n)\approx
\sum_{\mathbf k\in K}\mu_{\mathbf k}\,\operatorname{Re}
\left[C_{\mathbf k}(z,\mathbf n)e^{i\mathbf k\cdot(x,y)}\right].
$$

Here, $\mathbf n$ is the local normal and $\mathbf k$ is a surface reciprocal-lattice vector. The current table retains **12 independent complex-coefficient modes, including the zero-frequency term**; multiplicities $\mu_{\mathbf k}$ account for conjugate modes. This expansion represents the spatial periodicity of a classical potential energy, not an electronic Bloch wavefunction.

The mode coefficients still depend on atomic height and normal direction. Height uses cubic interpolation; directional coordinates use bicubic Hermite interpolation, storing coefficient values and the required derivatives. The current height grid has 0.1 Å spacing and 95 layers, including interpolation padding. Each directional table uses 13 × 13 nodes.

### 3. Cover molecular rotation with three directional charts

Describing a normal only by $(n_x/n_z,n_y/n_z)$ becomes problematic when a fragment turns nearly vertical and the denominator approaches zero. Version 2 uses three overlapping directional charts:

$$
(n_x/n_z,n_y/n_z),\quad
(n_y/n_x,n_z/n_x),\quad
(n_z/n_y,n_x/n_y).
$$

Each chart interpolates within its valid region. In overlapping regions, normalized smooth weights combine the energies. Together, the charts cover the orientations of nondegenerate unit normals, reducing frequent, expensive reference evaluations for curved molecules such as C60 whose normals may lie nearly parallel to the surface.

Both coordinates of each chart range from −1.2 to 1.2, with node spacing of 0.2. Chart-switching derivatives are included when differentiating the energy, and local-normal chain-rule derivatives propagate back to the atoms defining the normal.

### 4. Return smoothly to reference evaluation at the limits

The directional proxy supports heights of 2.4–11.2 Å, with a 0.3 Å transition region at each boundary. Outside this range, or when normal geometry approaches degeneracy, the relevant directional terms return to frozen-lattice reference evaluation instead of extrapolating the table.

The transition is defined at the scalar-energy level. If $w$ is a position- and geometry-dependent weight:

$$
U=wU_p+(1-w)U_r,
$$

the force must include the derivative of that weight:

$$
\mathbf F=w\mathbf F_p+(1-w)\mathbf F_r-(U_p-U_r)\nabla w.
$$

The implementation retains this term and the relevant height, chart, and local-normal derivatives. Simply blending two forces while omitting derivatives of the switching weights would make force and energy inconsistent. “Reference” here means the project's frozen SAIP model, not an exact description of a real gold surface.

### 5. Apply the 1.5× lattice corrugation factor

The corrugation factor strengthens differences between lateral positions; it does not simply multiply the entire adsorption energy by 1.5. For an atom's surface energy $U_i$, the implementation first averages six translated positions at the same height and normal:

$$
\bar U_i=\frac16\sum_{a=0}^{1}\sum_{b=0}^{2}
U_i(x_i+aL_x/2,\,y_i+bL_y/3,\,z_i,\,\mathbf n_i),
$$

then uses:

$$
U_i^{(\lambda)}=U_i+(\lambda-1)h(z_i)(U_i-\bar U_i).
$$

$L_x$ and $L_y$ are the surface periods, and the default is $\lambda=1.5$. The height-protection function $h$ is zero below 2.0 Å, one above 2.6 Å, and smoothly connected in between, avoiding direct amplification of the corrugation correction at very short distances. The six-point average is this implementation's filtering convention, not an exact continuous average over the entire unit cell.

Translations change only phase in the Fourier representation. The accelerated C/H path can therefore combine the six translated contributions within mode evaluation, reducing repeated interpolation. Radial-cache paths for other elements retain the corresponding translated evaluations. The corrugation correction and its height derivative are both included in the forces.

### 6. Performance benefits and limits

The default path replaces many repeated Au-neighbor sums with a small number of mode sums, grid interpolations, and local-normal calculations. The browser evaluates these in double-precision C++/WASM, and potential tables are downloaded with lossless compression. The native version has corresponding compiled kernels.

This remains an approximation with finite modes and finite grids. The proxy reduces surface-force evaluation costs; UFF, molecular preparation, normal handling, and reference fallback still require computation, so different molecules do not receive the same speedup. Energy, forces, and race statistics should be compared under matched conditions; visual similarity alone does not establish equivalence.

Main implementation entry points:

- [Frozen periodic surface and reference forces](molecular_car/periodic_au.py)
- [v2 table generation and directional charts](molecular_car/fourier_atlas.py)
- [Native compiled force kernel](molecular_car/fourier_atlas_kernel.py)
- [Browser energy and force implementation](wasm/src/potentials.cpp)
- [Browser potential-table export](wasm/tools/export_fourier_atlas.py)

## Model boundaries and validation

The model includes intramolecular flexibility, surface lattice corrugation, local directional dependence, external-field response, and thermal fluctuations. It does not currently include dynamic Au atoms, phonon propagation, crystal defects, reactive bonding, charge transfer, self-consistent metal polarization, or electronic excitation by a scanning tunneling microscope. Heat-bath friction also does not explicitly depend on molecular contact area. Performance differences between structures should be interpreted within this game model.

Development checks cover energy and finite-difference forces, local-normal chain rules, periodic translations, interpolation and switching boundaries, and energy drift with both the driving field and heat bath disabled. Timestep comparisons and complete-race comparisons are also available. See the [development and MD validation notes](DEVELOPMENT_NOTES.md) for historical records and the [WASM documentation](wasm/README.md) for current build and validation entry points. These checks detect implementation errors and assess approximations; they do not replace experimental validation or guarantee frame-by-frame agreement of long stochastic trajectories across platforms.

This is an open-source preview. Local end-to-end flows have been checked in Chrome, Edge, and Firefox, and the public HTTPS site has passed 23 workflow and caching checks. See [release acceptance](docs/RELEASE_ACCEPTANCE.md). Physical phones and Mac/Safari still require validation. Reloading loses race trajectories that have not been downloaded; workshop projects should also be downloaded as backups.

## Development and native use

The same repository maintains both browser and native versions:

| Directory | Contents |
| --- | --- |
| `molecular_car/` | Native Python engine, potential models, and molecular preparation |
| `wasm/` | Browser C++ engine, build entry points, and potential-table export |
| `frontend/` | Track, workshop, and third-party editor resources |
| `web-ui/` | Game interface source and assets; `dist/` is part of the current interface entry point |
| `tests/` | Regression tests and numerical checks |
| `scripts/` | Build, acceptance, performance, and parameter-validation tools |

To build the browser version, start with the [build and release guide](RELEASING.md); pinned dependencies and commands are in the [WASM documentation](wasm/README.md). Editor modifications have separate [Ketcher build instructions](frontend/editor/ketcher-source/README.md) and [JSmol source notes](frontend/editor/jsmol-source/README.md).

The native Python version requires Windows, a compatible 64-bit CPython 3.11 installation, and the LAMMPS GUI distribution:

```powershell
.\scripts\setup.ps1 -PythonExe 'C:\path\to\python.exe' -LammpsHome 'C:\path\to\LAMMPS'
.\start-game.ps1
```

Open the local address shown by the launcher. Browser players do not need to perform these installation steps.

The repository retains build inputs and validation tools. Large trajectories, temporary scan outputs, local virtual environments, and toolchain caches are excluded from source control. Prebuilt website packages and corresponding third-party sources are available through [Releases](https://github.com/chc-cyber/racing-molecule/releases).

### Where is LAMMPS?

- **The `main` source branch** contains integration code, pinned versions, and build scripts, without copying in the complete upstream LAMMPS source tree or desktop installer. Browser builds download and verify the specified sources and compile them under the ignored `.local/wasm/` directory.
- **GitHub Pages and the browser release package** include a WASM module linking LAMMPS, RDKit, and the project's force routines. Players can run it by loading the website, without installing LAMMPS separately.
- **GitHub Releases** provide the corresponding LAMMPS source archive used for the build and related build materials. The native Python entry point uses a separately installed copy of LAMMPS.

## License

Project-owned code is licensed under [MIT](LICENSE). Third-party components and modifications retain their respective licenses. The linked WASM module containing LAMMPS is distributed under GPLv2; MIT does not cover the entire release package. Releases provide corresponding sources, build inputs, and checksums. See the [third-party notices](THIRD_PARTY.md) and [workshop third-party notices](frontend/editor/THIRD_PARTY.md).
