Instructions to generate data and image for Figure 5.

1. Compile `meta3.c` with custom libraries: `randomGenerator.c`, `ran4.c` and `colors.c`.
2. Run `meta3.c` to generate data to the standard output. Save data to a `.dat` file, e.g. `TYPE4_NN10_N20_(1).dat`.
2. Recompile and run with 10 iterations for each `N` value `(20, 40, 60, 80, 100)`
   by changing the value of the macro `#define N` in `meta3.c`.
3. Each simulation's generated data should be saved in separate .dat files, e.g. `TYPE4_NN10_N20_(2).dat`, `TYPE4_NN10_N20_(3).dat`, etc.
   Pre-generated data files are provided as `TYPE4_NN10_N*_(*).dat`.
4. Run the *Mathematica* code `plotter.nb` (optimized for *Mathematica* 11.0.0.0). It reads the generated data files and creates the figure.
   Pre-generated image is provided as `combined.tiff`.

(Reference name: `...NRH\Storkey\meta\TYPE4_NN10_Nchg_P10_me`)

By István Zachar
2016
