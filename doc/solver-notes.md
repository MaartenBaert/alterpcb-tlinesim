Transmission Line Solver Notes
==============================

The solver uses the [finite element method](https://en.wikipedia.org/wiki/Finite_element_method) to solve the [macroscopic formulation of Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations#Macroscopic_formulation) in differential form. The solver does not actually calculate the electric and magnetic field directly, instead it uses the [potential field approach](https://en.wikipedia.org/wiki/Mathematical_descriptions_of_the_electromagnetic_field#Potential_field_approach). This greatly reduces the number of unknowns and also simplifies matrix generation, especially with regard to the handling of conductors and boundary conditions.

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

Currently the solver avoids these problems by using ideal conductors for field calculations, and then adding the losses caused by the skin effect as a post-processing step. All currents flow on the surfaces of the conductors rather than inside the conductors, however the surfaces do have a finite resistance, which is calculated based on the conductivity and the skin depth. This approximation works very well when the skin depth is significantly smaller than the thickness of the conductors, which is generally the case at RF frequencies, but it breaks down at lower frequencies. In particular, the resistance of copper tracks may appear to be smaller than the DC resistance based on the cross-sectional area. I am planning to address this by also calculating the DC resistance, and adjusting the AC resistance where necessary.

In the future, I might add a complete skin effect solver as an optional feature, but for most applications this should not be necessary.

### Quasi-TEM approximation

The solver assumes that the electric and magnetic field are always orthogonal to the direction of propagation, also known as a TEM mode. However in reality, this assumption is true only when the entire simulation domain is filled with a single dielectric with constant properties. This is the case for striplines, but many other types of transmission lines like microstrip or coplanar waveguides have their electromagnetic fields located partially inside the substrate, and partially in the air above the substrate. As a result, the field will not be perfectly orthogonal to the direction of propagation. At low frequencies, i.e. when the cross-section of the transmission line is significantly smaller than the wavelength, the resulting mode will still be very similar to a TEM mode, so it is often called a quasi-TEM mode.

The solver implicitly assumes that all fields are orthogonal to the direction of propagation, even for quasi-TEM modes. At very high frequencies, where the transmission line may start to behave more like a waveguide, this assumption is no longer accurate. This affects some types of transmission lines more than others: transmission lines which have their field confined in a small area (like coplanar waveguides) are far less affected than transmission lines which have a more spread out field (like microstrip). Differential transmission lines are generally less affected than single-ended transmission lines.

Currently the solver makes no attempt to determine whether the quasi-TEM assumption is still accurate at the required frequency. I am planning to add something like this in the future. But for now, you should manually verify that the area around the transmission line which contains most of the electromagnetic energy is significantly smaller than the wavelength.

Boundary condition handling
---------------------------

Because of the way the potential equations are set up, Dirichlet boundary conditions correspond to perfect electric conductors, and Neumann boundary conditions correspond to perfect magnetic conductors.

Matrix generations happens cell-by-cell rather than node-by-node. Within each cell, the equations for each pair of nodes are added to the matrix. This is not terribly efficient because there will be a lot of duplicate calculations, but this approach makes it very easy to deal with varying cell sizes (or even mixed cell types), discontinuities caused by changes in material properties, and boundary conditions. Cells that are not part of the simulation domain (i.e. cells outside the boundaries or inside conductors) are simply ignored. These cells will implicitly become a perfect electric conductor if the potential values of the surrounding nodes are fixed, or a perfect magnetic conductor if the potential values of the surrounding nodes are free variables. Floating conductors work just like perfect electric conductors, however their potential is a free variable rather than a fixed value. The corresponding equation enforces zero total charge/current within the floating conductor. When the same free variable is assigned to all boundary nodes of the conductor, this equation is constructed implicitly without requiring any extra code.

Eigenmode decomposition
-----------------------

(TODO)
