# Metodi Numerici Avanzati per Equazioni alle Derivate Parziali

Codes for the examination of the **Metodi Numerici Avanzati per Equazioni alle Derivate Parziali**[^1] course at **UniMiB**.

[^1]: Advanced Numerical Methods for Partial Differential Equations.

Implementation of the Crouzeix-Raviart Element to solve the Stokes problem with homogeneous Dirichlet boundary data.

## Contents

- `/src/*`
    - `/src/solver.m` Solver implementation for the Stokes problem on the mesh.
    - `/src/basis.m` Basis functions for velocity.
    - `/src/gradbasis.m` Computation of gradients of basis functions for velocity.
    - `/src/pressurebasis.m` Basis functions for pressure.
    - `/src/loading.m` Source definition.
    - `/src/exact.m` Implementation of the exact solution for the specified source.
    - `/src/quadrature.m` Quadrature nodes for the reference element.
    - `/src/errorTrend.m` Generation of plots to visualize error trends.
    - `/src/meshinfo.m` Sizes and DOFs info.
    - `/src/reload.m` Wrapper for the errorTrend function.
    - `/src/quadmeshes.mat` Predefined mesh.
    