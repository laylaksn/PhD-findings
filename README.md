# PhD-findings
Thesis and code to be uploaded soon...

# Summary

This thesis is concerned with the efficient approximation of solutions to high-dimensional
partial differential equations (PDEs) using the Proper Generalised Decomposition (PGD).
The PGD represents unknown fields as a series of low-dimensional problems, which can
overcome the computational challenge associated with the ‘curse of dimensionality’. However, convergence of PGD algorithms is not guaranteed, particularly when using Galerkin
formulations.
To address these limitations, PGD algorithms based on rigorously derived Least-Squares
formulations are developed based on optimal first-order reformulations of the problems. The
theoretical framework is based on ADN theory, which is used to address practical issues
and convergence considerations. A novel approach to solving elliptic PDEs is presented
based on a PGD which solves the minimisation of a discrete rather than a continuous Least-Squares functional. The spectral element method which combines the flexibility of finite
element methods with the accuracy of spectral methods is used as the foundation of the
discretisation.
The application of the Least-Squares PGD algorithm to the Poisson, convection-diffusion
and the Stokes problems demonstrates the efficiency, robustness and precision of this novel
approach. Notably, the algorithm maintains high accuracy for convection-dominated problems, which are prone to numerical instability.
The optimal Least-Squares PGD is extended to incorporate a mimetic framework for
solving the Poisson problem. This ensures that the algebraic equations obtained through
discretisation are compatible and physically consistent, providing accurate numerical representations of the physical solution. Results validate, for the first time, the successful
implementation of the mimetic approach within a Least-Squares PGD algorithm. The performance and properties of the Least-Squares PGD algorithm on the series of benchmark
problems studied in this thesis demonstrate the efficiency and accuracy of this approach for
solving PDEs and provide a solid foundation for further development of this technique.
