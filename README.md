# PhD Dissertation - Jeremy L Thompson

## Organization

The raw LaTeX source for this dissertation can be found in the `thesis` and `presentation` folders.
The images and charts can be found in the `img` folder.
The data and scripts for all charts and analysis in this dissertation can be found in the `data`, `scripts`, and `jupyter` folders.

## Reproducibility

The software used in this dissertation is all open source and freely available.

The Local Fourier Analysis can be reproduced with the Julia package [LFAToolkit.jl](https://www.github.com/jeremylt/LFAToolkit.jl). LFAToolkit.jl is a toolkit for analyzing the performance of preconditioners for arbitrary, user provided weak forms of second order partial differential equations.

The high-order matrix-free finite element operators and Neo-Hookean hyperelasticity mini-app can be found in the [libCEED](https://www.github.com/CEED/libCEED) GitHub repository.
libCEED is a low-level API library for the efficient high-order discretization methods developed by the ECP co-design Center for Efficient Exascale Discretizations (CEED). While the focus is on high-order finite elements, the approach used in libCEED is mostly algebraic and thus applicable to other discretizations in factored form.

Finally, the linear and non-linear solver and preconditioning infrastructure can be found in the scientific toolkit [PETSc](https://www.mcs.anl.gov/petsc/).
PETSc is a suite of data structures and routines for the scalable, parallel solution of scientific applications modeled by partial differential equations. It supports MPI, and GPUs through CUDA, HIP, or OpenCL, as well as hybrid MPI-GPU parallelism.
