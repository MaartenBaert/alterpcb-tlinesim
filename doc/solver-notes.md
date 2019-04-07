Transmission Line Solver Notes
==============================

The solver uses the [finite element method](https://en.wikipedia.org/wiki/Finite_element_method) to solve the [macroscopic formulation of Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations#Macroscopic_formulation) in differential form. The [Galerkin method](https://en.wikipedia.org/wiki/Galerkin_method) of [mean weighted residuals](https://en.wikipedia.org/wiki/Method_of_mean_weighted_residuals) is used to discretize the problem. The solver does not actually calculate the electric and magnetic field directly, instead it uses the [potential field approach](https://en.wikipedia.org/wiki/Mathematical_descriptions_of_the_electromagnetic_field#Potential_field_approach). This greatly reduces the number of unknowns and also simplifies matrix generation, especially with regard to the handling of conductors and boundary conditions.

Limitations
-----------

### Modeling errors

Simulation results are only as accurate as the simulation model itself. Manufacturing variations will cause small differences in track width and thickness, solder mask thickness, substrate thickness, substrate permittivity and loss tangent, etc. In many cases the nominal value of certain parameters is not even known exactly. For example, PCB substrate material datasheets will typically specify just one value for the permittivity and loss tangent of the material, or possibly 2-4 values at different frequencies if you are lucky. In reality, most PCB substrate materials (like FR4) are anisotropic, which means that the permittivity and loss tangent is different in the horizontal and vertical direction. Furthermore, the permittivity and loss tangent depend on not only the frequency, but also the weave (the pattern of glass fiber strands) and resin content (which depends on the lamination process used by the PCB manufacturer for boards with more than 2 copper layers).

It's worth noting that FR4 is not really a material, it refers to a wide range of materials produced by different manufacturers with similar but not identical properties. For this reason you won't find FR4 in the material database. Instead you should find out which type of FR4 your PCB manufacturer uses, and design your transmission lines for this specific type.

The included material database contains material properties for various PCB substrate materials from a range of manufacturers. If you are using special RF substrates such as RO4350, the material properties will be reasonably accurate, because the material properties are well specified even at very high frequencies. However if you are using cheaper non-RF materials, the material properties may deviate significantly from the values in the database. There's a lot of guesswork involved in deriving these properties, for example the level of anisotropy is often derived from other (similar) materials because the actual anisotropy is not known. Similarly, if the material datasheet only specifies dielectric properties at a very low frequency (e.g. 1 MHz), the value at higher frequencies (e.g. 1 GHz) is extrapolated based on similar materials, and may be wildly inaccurate.

Another minor source of modeling errors is the track shape. Real tracks don't have a perfectly rectangular cross-section, the actual shape is usually more trapezoidal. This doesn't have a huge impact on the transmission line properties though, unless the track width is extremely small.

Due to manufacturing variations, it is not a good idea to design sensitive transmission lines with a track width or spacing too close to the minimum value allowed by the manufacturer. If the minimum track width is 100 µm, the actual track width might end up being 80 µm or 120 µm. The manufacturer doesn't guarantee that the width or spacing will be accurate, only that there will be no open or short circuits. So if you have sufficient space on your PCB, it's a good idea to pick a track width and spacing which is 2-4 times larger than the minimum value, at least if you care about accuracy.

### Discretization errors

Electic and magnetic potentials inside the cells are determined by linear interpolation of the potential values at the surrounding nodes. Smaller cells will improve the accuracy but also the processing time. The mesh generation algorithm tries to place smaller cells in areas where strong fields are expected in order to minimize both the error and the processing time.

### Skin effect approximation

Due to the [skin effect](https://en.wikipedia.org/wiki/Skin_effect), currents tend to flow in a very thin layer just below the surface of conductors. For example, the skin depth of copper is just 2.06 µm at a frequency of 1 GHz. Simulating these currents is challenging for a number of reasons. Since the skin depth is so small, very small cells are required to properly simulate the current distribution, which increases the processing time. It would be impractical to fill the entire conductor with such small cells, so instead only a thin layer of cells should be added (several times the skin depth). Since the skin depth is frequency-dependent, the mesh would need to be regenerated for each frequency, which slows down frequency sweeps. Furthermore, the addition of lossy (non-ideal) conductor cells results in a complex-valued, non-Hermitian matrix which must be solved with LU decomposition rather than Cholesky decomposition, which again makes the solver a lot slower.

Currently the solver avoids these problems by using ideal conductors for field calculations, and then adding the losses caused by the skin effect as a post-processing step. All currents flow on the surfaces of the conductors rather than inside the conductors, however the surfaces do have a finite resistance, which is calculated based on the conductivity and the skin depth. This approximation works very well when the skin depth is significantly smaller than the thickness of the conductors, which is generally the case at RF frequencies, but it breaks down at lower frequencies. In particular, the resistance of copper tracks may appear to be smaller than the DC resistance based on the cross-sectional area. In order to avoid these unrealistic results, the solver also calculates the DC resistance of each track and uses this as a lower bound for the AC resistance.

In the future, I might add a complete skin effect solver as an optional feature, but for most applications this should not be necessary.

### Quasi-TEM approximation

The solver assumes that the electric and magnetic field are always orthogonal to the direction of propagation, i.e. it assumes that the propagation mode will be a TEM mode. However in reality, this assumption is only true when the entire simulation domain is filled with a single dielectric with constant properties. This is the case for striplines (at least when the dielectric material is isotropic), but many other types of transmission lines like microstrip or coplanar waveguides have their electromagnetic fields located partially inside the substrate, and partially in the air above the substrate. As a result, the field will not be perfectly orthogonal to the direction of propagation. At low frequencies, i.e. when the cross-section of the transmission line is significantly smaller than the wavelength, the resulting mode will still be very similar to a TEM mode, so it is often called a quasi-TEM mode.

The solver implicitly assumes that all fields are orthogonal to the direction of propagation, even for quasi-TEM modes. At very high frequencies, where the transmission line may start to behave more like a waveguide, this assumption is no longer accurate. This affects some types of transmission lines more than others: transmission lines which have their field confined in a small area (like coplanar waveguides) are far less affected than transmission lines which have a more spread out field (like microstrips). Differential transmission lines are generally less affected than single-ended transmission lines.

Currently the solver makes no attempt to determine whether the quasi-TEM assumption is still accurate at the required frequency. I am planning to add something like this in the future. But for now, you should manually verify that the area around the transmission line which contains most of the electromagnetic energy is significantly smaller than the wavelength.

Boundary condition handling
---------------------------

Because of the way the potential equations are set up, Dirichlet boundary conditions correspond to perfect electric conductors, and Neumann boundary conditions correspond to perfect magnetic conductors.

Matrix building happens cell-by-cell rather than node-by-node. Within each cell, the equations for each pair of nodes are added to the matrix. This is not terribly efficient because there will be a lot of duplicate calculations, but this approach makes it very easy to deal with varying cell sizes (or even mixed cell types), discontinuities caused by changes in material properties, and boundary conditions. Cells that are not part of the simulation domain (i.e. cells outside the boundaries or inside conductors) are simply ignored. These cells will implicitly become a perfect electric conductor if the potential values of the surrounding nodes are fixed, or a perfect magnetic conductor if the potential values of the surrounding nodes are free variables. Floating conductors work just like perfect electric conductors, however their potential is a free variable rather than a fixed value. The corresponding equation enforces zero total charge/current within the floating conductor. When the same free variable is assigned to all boundary nodes of the conductor, this equation is constructed implicitly without requiring any extra code.

Eigenmode analysis
------------------

For many types of transmission lines, signals can propagate in more than one 'mode'. The most well-known example is the 'common mode' and 'differential mode' of a differential transmission line. This concept can be extended to more complex transmission lines as well. Since any linear combination of two or more modes is itself also a valid mode, there are an infinite number of possible propagation modes. In order to determine the properties of these propagation modes, we can determine the so-called 'eigenmodes' of the transmission line. Eigenmodes are propagation modes which do not interact with each other - they can be considered 'pure' modes. All other modes can be expressed as linear combinations of the eigenmodes of the transmission line.

In the case of a symmetric differential transmission line, the eigenmodes correspond to 'common mode' and 'differential mode'. In general, when we consider only TEM and quasi-TEM modes, the number of eigenmodes is equal to the number of conductors minus one. For example, a transmission line consisting of 4 tracks above a ground plane has 5 conductors, and will have 4 eigenmodes which will likely look something like this:

	Mode 1: ++++
	Mode 2: ++--
	Mode 3: +-+-
	Mode 4: +--+

Each eigenmode has a corresponding propagation constant that determines how fast the signal propagates and how quickly the amplitude decreases due to losses. Other modes which are actually linear combinations of eigenmodes with different propagation constants will disperse over time, since the different underlying eigenmodes propagate at different speeds. This is why we generally want to use only eigenmodes when we use transmission lines to transmit data.

### Calculating eigenmodes

Determining the eigenmodes of a transmission line boils down to solving a large sparse eigenvalue problem, where the eigenvalues correspond to the propagation constants and the eigenvectors correspond to the electromagnetic fields. Normally this would be a very computationally intensive process, however the solver is able to decrease the processing time drastically by taking advantage of the fact that we are only interested in quasi-TEM modes. This reduces the problem to the following steps:
- Calculate the electric and magnetic potential (requires solving two large sparse linear systems).
- Derive charge and current distributions from the electric and magnetic potentials (requires solving a smaller sparse linear system for surface currents).
- Calculate dielectric losses based on the electric potential.
- Calculate resistive losses based on the surface current distribution.
- Derive the inductance (L), capacitance (C), resistance (R) and conductance (G) matrices that are used by the telegrapher's equations. These are square matrices with sizes corresponding to the number of eigenmodes.
- Determine the eigenmodes and corresponding propagation constants by calculating the eigenvalue decomposition of (R + j * omega * L) * (G + j * omega * C). The propagation constants correspond to the square root of the eigenvalues.
The advantage of this approach is that the eigenvalue problem that needs to be solved is tiny, which makes this method much easier. The downside is that only quasi-TEM modes can be found this way.

### Transmission line coupling

Eigenmodes can be used to analyze the coupling between two transmission lines. When two transmission lines are in close proximity, they should be treated as a single transmission line with twice as many eigenmodes. The eigenmodes of the combined transmission line are called supermodes. When the space between the transmission lines is large, the supermodes will be nearly indistinguishable from the original eigenmodes. However when they are brought closer together, the interaction becomes stronger, and the supermodes may look very different from the original eigenmodes.

By decomposing each original mode into a linear combination of supermodes, it is possible to analyze the propagation of each supermode separately. When the supermodes reach the other end of the transmission line, they are converted back to the original modes. Since different supermodes will likely have different propagation constants, there may be a time difference (or phase shift in the frequency domain) between the supermodes. As a result, the reconstruction of the original modes won't be perfect and some energy will 'leak' into other modes. This is known as mode coupling. An important conclusion is that when all supermodes have the same propagation constant, there will be no mode coupling. This is the case for pure TEM supermodes, which occur when the entire simulation domain has a constant permittivity, which is the case for striplines in an isotropic dielectric material.

It may sound very counterintuitive that two transmission lines which are placed very close together don't couple at all simply because the propagation constants of the supermodes are the same, but this is actually true. There is actually an intuitive explanation for this phenomenon. When two tracks are in close proximity, there will be both capacitive coupling (mutual capacitance) and inductive coupling (mutual inductance). The capacitive coupling induces a voltage in the other track, however the inductive coupling induces a current in the opposite direction in the other track. If the ratio matches the characteristic impedance of the track, the two effects will cancel out, and there will be no coupling. This happens to be the case for pure TEM transmission lines.

Aside from mode coupling, there is also another type of coupling known as butt coupling. It occurs when a transmission line abruptly changes to a different type, thus forcing a sudden change in propagation mode. This occurs for example near the start and end of a coupled section of a transmission line:

	            butt coupling           butt coupling
	                  ┊                       ┊
	                  ┊ ◁┈┈ mode coupling ┈┈▷ ┊
	                  ▽                       ▽
	┎───────────────────────────────────────────────────────────┒
	┖───────────────────────────────────────────────────────────┚
	                  ╭───────────────────────╮
	                  │ ╭───────────────────╮ │
	                  │ │                   │ │
	                  ┕━┙                   ┕━┙

Butt coupling is outside the scope of this solver, because is a lot harder to analyze. It is essentially a 3D problem whereas mode coupling is a 2D problem. Luckily the effect of butt coupling is often quite small, at least at low frequencies, so in many cases it can be ignored. Butt coupling becomes more significant at higher frequencies, or when the coupled section is short, or when there is no significant mode coupling (e.g. between striplines).

Mode coupling always occurs in the same direction of propagation, i.e. it will never induce a signal that propagates in the opposite direction, neither in the original track nor in the coupled tracks. Butt coupling on the other hand will usually induce signals that propagate in both directions, however the forward-propagating and backward-propagating signals may be very different due to constructive and destructive interference. The well-known quarter-wavelength directional coupler is based on this principle: the coupled signal which travels in the opposite direction is the result of butt coupling, not mode coupling. Since the two transitions are opposites, they produce signals with opposite polarity, which causes the forward-propagating signal to cancel out. However since the transitions are spaced by a quarter wavelength, the backward-propagating signals have a 180-degree phase shift which causes them to add up constructively. This is a good example of an application where butt coupling clearly can't be ignored. For most applications, you should probably focus on mode coupling when signals propagate in the same direction, and focus on butt coupling when signals propagate in opposite directions.

It's also worth noting that both types of coupling can become a lot worse when a transmission line isn't properly terminated (matched). Without good matching, the signals will reflect back and forth multiple times, which gives them multiple chances to couple with the other transmission lines, possibly creating constructive interference. This can be avoided by terminating at least one side of the transmission line, or ideally both sides. If you are dealing with high-speed signals and receiver-side termination isn't practical (e.g. because of high current consumption for regular CMOS logic levels), adding some series resistors as transmitter-side termination is a good idea.

### Characteristic impedance

Propagation modes in a quasi-TEM transmission line have a voltage component and a current component. The voltage to current ratio is known as the characteristic impedance of the propagation mode. This definition works fine for simple transmission lines with only two terminals, or symmetric transmission lines with three terminals (differential pairs). Unfortunately it doesn't work in general for more complicated transmission lines with multiple propagation modes, because there is no guarantee that the voltage and current across the terminals will be distributed in the same way, not even for eigenmodes. For example, applying the voltages (1, 0, 0, 0) may result in currents (0.1, -0.02, -0.03, -0.05). It doesn't really make sense to specify a voltage to current ratio when this ratio is different for every pair of terminals. What we can do instead is determine a characteristic impedance matrix, which describes the relation between the voltage and current 'vector' for all terminals simultaneously.

Although the characteristic impedance matrix describes the behavior of the transmission line perfectly, it is not very useful for the designer. It would be more convenient to extract an approximate characteristic impedance for each mode and describe the coupling by other means (e.g. S-parameters). We could extract the diagonal elements of the characteristic impedance matrix and ignore the rest, but this usually results in an overestimation, because it corresponds to leaving the unused tracks floating, which reduces the capacitance and increases the inductance. We could instead use the diagonal elements of the characteristic admittance matrix (the inverse of the characteristic impedance matrix), but this results in an underestimation, because it corresponds to shorting the unused tracks to ground, which increases the capacitance and reduces the inductance. The solver uses a simple compromise: it calculates the impedance and admittance matrices, extracts the diagonal elements, and uses the geometric mean of the two impedance values obtained this way. Usually the resulting values are very close to optimal in terms of minimizing the amplitude of reflections caused by impedance mismatch, so they should be adequate for practical transmission line design.

### Propagation constants

Propagation constants of eigenmodes can be used to determine the speed of propagation as well as the attenuation of the signal. However the designer may want to know these properties even for modes that aren't true eigenmodes. In this case the solver uses the Rayleigh quotient to determine an approximate propagation constant.

TODO:
- propagation constant (coupled)
- coupling S-parameters?
- iterative RF mode solver (Rayleigh quotient iteration)

