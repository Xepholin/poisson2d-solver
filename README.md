# Poisson2D Solver

Création d'une matrice creuse via un stencil, avec stockage en CSR.

Stencil:
```math
    \begin{pmatrix}
       & -1 &    \\
    -1 &  4 & -1 \\
       & -1 &  
    \end{pmatrix}
```

Puis résolution du système avec les méthodes de Jacobi et de Gauss-Seidel.

# Usage
```sh
$ make
$ ./poisson N MAX_ITER THRESHOLD
```
