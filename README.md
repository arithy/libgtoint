# LibGtoint #

## Overview ##

**LibGtoint** is an analytical GTO (Gaussian type orbital) integral library for C and Fortran.

> *Firstly note that this library is not fast and currently it is just an implementation for study. It would be much appreciated if you could share your expertise on more efficient implementation by [Issues](https://github.com/arithy/libgtoint/issues), [Pull Requests](https://github.com/arithy/libgtoint/pulls), or [Wiki](https://arithy.org/libgtoint/wiki/).*

Its main features are as follows:
- Supports Cartesian GTOs and spherical GTOs,
- Supports GTOs with any angular momentum quantum numbers (s, p, d, f, ...),
- Supports analytical derivatives with any orders [*1], and
- Implemented compactly.

It supports the following integrals:
- overlap integrals,
- kinetic energy integrals,
- nuclear attraction integrals,
- electron repulsion integrals (i.e. two-electron integrals),
- multipole moment integrals, and
- scalar ECP (effective core potential) integrals.

Currently, LibGtoint uses only Obara-Saika scheme.

[*1]: As for ECP integrals, there is a limitation. For the details, see *Notice* in the section *API*.

## Installation ##

To install LibGtoint, [CMake](https://cmake.org/) 3.14 or higher is required to be installed in your system.

### For Unix-like OS ###

If you use Unix or Unix-like OS such as Linux and macOS, you can install LibGtoint by executing the following commands:

```sh
mkdir build
cd build
cmake ..
make
make test
sudo make install
```

To use specific C and Fortran compilers, use the `cmake` options `-DCMAKE_C_COMPILER=`*C-Compiler*  and `-DCMAKE_Fortran_COMPILER=`*Fortran-Compiler*  respectively.

If no need of the Fortran interface, use the `cmake` option `-DFortran=OFF`.

By default, a static library is built. If you need a shared object version of this library, use the `cmake` option `-DBUILD_SHARED_LIBS=1`.

The default installation directory is `/usr/local`. If you want to change it, use the `cmake` option `-DCMAKE_INSTALL_PREFIX=`*Installation-Path*.

### For Windows ###

If you use Windows, you have two options shown in the sections below.

#### Using Visual Studio ####

To install LibGtoint which is built by using Microsoft's genuine build tools, [Build Tools for Visual Studio](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019) has to be installed in your system.
After you have installed it, you can install LibGtoint by executing the following commands using 'Developer Command Prompt for VS 2019' or 'Developer PowerShell for VS 2019':

```sh
mkdir build
cd build
cmake -DFortran=OFF -DCMAKE_INSTALL_PREFIX=%USERPROFILE%\libgtoint ..
MSBuild ALL_BUILD.vcxproj -p:Configuration=Release
MSBuild RUN_TESTS.vcxproj -p:Configuration=Release
MSBuild INSTALL.vcxproj -p:Configuration=Release
```

By default, a static library is built. If you need a dynamic link library, use the `cmake` option `-DBUILD_SHARED_LIBS=1`.

In the above commands, the installation directory is set to that directly under your account directory, ex, `C:\Users\`*Your-User-Name*`\libgtoint`. If you want to change it, modify the `cmake` option `-DCMAKE_INSTALL_PREFIX=`*Installation-Path*.

**Known Issues**:
- As of `cmake` 3.20.1, `cmake` fails to generate a build configuration for compiling Fortran sources even if Intel Fortran compiler for Windows exists. So, you cannot help but specify the `cmake` option `-DFortran=OFF`.
- As of `cmake` 3.20.1, `cmake` does not see C compilers other than Microsoft Visual C++ compiler (MSVC) even if Intel C/C++ compiler for Windows exists. So, the `cmake` option `-DCMAKE_C_COMPILER=icl` has no effect.

#### Using MinGW-w64 ####

To install LibGtoint which is built by MinGW tool chain, you can install LibGtoint by executing the following commands:

```sh
mkdir build
cd build
cmake -G "MSYS Makefiles" -DCMAKE_INSTALL_PREFIX=/usr/local ..
make
make test
make install
```

To use specific C and Fortran compilers, use the `cmake` options `-DCMAKE_C_COMPILER=`*C-Compiler*  and `-DCMAKE_Fortran_COMPILER=`*Fortran-Compiler*  respectively.

If no need of the Fortran interface, use the `cmake` option `-DFortran=OFF`.

In the above commands, the installation directory is set to `/usr/local`. If you want to change it, modify the `cmake` option `-DCMAKE_INSTALL_PREFIX=`*Installation-Path*.

**Known Issue**:
- As of GNU Binutils 2.36.1, a DLL version of LibGtoint causes abnormal termination. So, do not specify the `cmake` option `-DBUILD_SHARED_LIBS=1`.

## API ##

### Preparation ###

To use the API in a C program, the header file `gtoint.h` is needed.

```c
#include <gtoint.h>
```

To use the API in a Fortran program, the module `gtoint` is needed.

```f90
use gtoint
```

### Basic types ###

In the API for C, two data types are defined.
- The data type `gtoint_double3_t` has three `double` type member variables named `x`, `y`, and `z`.
- The data type `gtoint_int3_t` has three `int` type member variables named `x`, `y`, and `z`.

### Error codes ###

As for C, most of API functions return error codes.
As for Fortran, most of API routines pass error codes via the argument `err`.

The error codes are shown below.
- `GTOINT_ERROR_OK`: successfully done with no error
- `GTOINT_ERROR_ARGUMENT`: invalid argument
- `GTOINT_ERROR_MEMORY`: out of memory
- `GTOINT_ERROR_UNSUPPORTED`: unsupported functionality
- `GTOINT_ERROR_INTERNAL`: internal error

### Creation of an integrator ###

To construct basis functions and ECPs, and to compute several kinds of integrals, an *integrator* is required.
It can be created by using the following function or routine.

- C
    ```c
    gtoint_error_t gtoint_integrator_create(gtoint_integrator_t *itg);
    ```
- Fortran
    ```f90
    subroutine gtoint_integrator_create(itg, err)
       type(C_PTR), intent(out) :: itg
       integer(C_INT), intent(out) :: err
    end subroutine
    ```

**Arguments**
- `itg`: the integrator to be created
- `err`: the error code (only Fortran)

**Return Value**
- the error code (only C)

### Disposal of an integrator ###

The integrator must be disposed by using the following function or routine when the integrator is no more needed.

- C
    ```c
    void gtoint_integrator_destroy(gtoint_integrator_t itg);
    ```
- Fortran
    ```f90
    subroutine gtoint_integrator_destroy(itg)
       type(C_PTR), intent(in) :: itg
    end subroutine
    ```

**Argument**
- `itg`: the integrator to be disposed; nothing is done if `GTOINT_NULL` (in C) or `c_null_ptr` (in Fortran)

### Duplicate of an integrator ###

The integrator can be duplicated by using the following function or routine.

- C
    ```c
    gtoint_error_t gtoint_integrator_copy(gtoint_integrator_t *itg, gtoint_integrator_t src);
    ```
- Fortran
    ```f90
    subroutine gtoint_integrator_copy(itg, src, err)
       type(C_PTR), intent(out) :: itg
       type(C_PTR), intent(in) :: src
       integer(C_INT), intent(out) :: err
    end subroutine
    ```

**Argument**
- `itg`: the integrator to be created by copying
- `src`: the integrator to be copied
- `err`: the error code (only Fortran)

**Return Value**
- the error code (only C)

### Cleanup of work memory ###

The work memory automatically grown in the integrator can be deallocated by using the following function or routine.

- C
    ```c
    void gtoint_integrator_cleanup_memory(gtoint_integrator_t itg);
    ```
- Fortran
    ```f90
    subroutine gtoint_integrator_cleanup_memory(itg)
       type(C_PTR), intent(in) :: itg
    end subroutine
    ```

**Argument**
- `itg`: the integrator whose work memory is to be deallocated

### Change of the integral error tolerance ###

The error tolerance of integrals can be changed by using the following function or routine.
The error tolerance is used for nuclear attraction integrals, electron repulsion integrals, and ECP integrals.
The default value is `1e-10`.

- C
    ```c
    void gtoint_integrator_set_error_tolerance(gtoint_integrator_t itg, double tol);
    ```
- Fortran
    ```f90
    subroutine gtoint_integrator_set_error_tolerance(itg, tol)
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: tol
    end subroutine
    ```

**Arguments**
- `itg`: the integrator whose error tolerance is to be changed
- `tol`: the new error tolerance

### Retrieval of the integral error tolerance ###

The error tolerance of integrals can be retrieved by using the following function or routine.
The error tolerance is used for nuclear attraction integrals, electron repulsion integrals, and ECP integrals.

- C
    ```c
    double gtoint_integrator_get_error_tolerance(gtoint_integrator_t itg);
    ```
- Fortran
    ```f90
    subroutine gtoint_integrator_get_error_tolerance(itg, tol)
       type(C_PTR), intent(in) :: itg
       real(8), intent(out) :: tol
    end subroutine
    ```

**Arguments**
- `itg`: the integrator whose error tolerance is to be retrieved
- `tol`: the current error tolerance (only Fortran)

**Return Value**
- the current error tolerance (only C)

### Change of the integral cutoff ###

The cutoff of integrals can be changed by using the following function or routine.
The cutoff is used for ECP integrals.
The default value is `1e-15`.

- C
    ```c
    void gtoint_integrator_set_cutoff(gtoint_integrator_t itg, double cut);
    ```
- Fortran
    ```f90
    subroutine gtoint_integrator_set_cutoff(itg, cut)
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: cut
    end subroutine
    ```

**Arguments**
- `itg`: the integrator whose cutoff is to be changed
- `cut`: the new cutoff

### Retrieval of the integral cutoff ###

The cutoff of integrals can be retrieved by using the following function or routine.
The cutoff is used for ECP integrals.

- C
    ```c
    double gtoint_integrator_get_cutoff(gtoint_integrator_t itg);
    ```
- Fortran
    ```f90
    subroutine gtoint_integrator_get_cutoff(itg, cut)
       type(C_PTR), intent(in) :: itg
       real(8), intent(out) :: cut
    end subroutine
    ```

**Arguments**
- `itg`: the integrator whose cutoff is to be retrieved
- `cut`: the current cutoff (only Fortran)

**Return Value**
- the current cutoff (only C)

### Construction of basis functions ###

The basis functions are defined using the following function or routine.
The basis functions are automatically normalized.

- C
    ```c
    gtoint_error_t gtoint_basis_shell_create(
        gtoint_basis_shell_t *bas, gtoint_integrator_t itg,
        double sf, int na, const int *a, int ng, const double *g, const double *c
    );
    ```
- Fortran
    ```f90
    subroutine gtoint_basis_shell_create(bas, itg, sf, na, a, ng, g, c, err)
       type(C_PTR), intent(out) :: bas
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: sf
       integer, intent(in) :: na
       integer(C_INT), intent(in) :: a(*)
       integer, intent(in) :: ng
       real(8), intent(in) :: g(*)
       real(8), intent(in) :: c(*)
       integer(C_INT), intent(out) :: err
    end subroutine
    ```

**Arguments**
- `bas`: the basis shell to be created
- `itg`: the integrator
- `sf`: the scaling factor (the squared value is multiplied to GTO exponents)
- `na`: the number of the shells
- `a`: the array of the angular momentum quantum numbers of the respect shells
    - the array dimension is `na`
    - to specify a spherical GTO, perform bitwise-not of the number using the operator `~` in C or the intrinsic function `not` in Fortran
- `ng`:the number of the contraction
- `g`: the array of the GTO exponents
    - the array dimension is `ng`
- `c`: the array of the contract coefficients of the normalized GTOs
    - the array dimensions are `[na][ng]` in C-style notation or `(ng,na)` in Fortran-style notation
- `err`: the error code (only Fortran)

**Return Value**
- the error code (only C)

### Destruction of basis functions ###

The basis functions must be destroyed by using the following function or routine when they are no more needed.

- C
    ```c
    void gtoint_basis_shell_destroy(gtoint_basis_shell_t bas);
    ```
- Fortran
    ```f90
    subroutine gtoint_basis_shell_destroy(bas)
       type(C_PTR), intent(in) :: bas
    end subroutine
    ```

**Argument**
- `bas`: the basis shell to be destroyed; nothing is done if `GTOINT_NULL` (in C) or `c_null_ptr` (in Fortran)

### Copy of basis functions ###

The basis functions can be copied using the following function or routine.

- C
    ```c
    gtoint_error_t gtoint_basis_shell_copy(gtoint_basis_shell_t *bas, gtoint_basis_shell_t src);
    ```
- Fortran
    ```f90
    subroutine gtoint_basis_shell_copy(bas, src, err)
       type(C_PTR), intent(out) :: bas
       type(C_PTR), intent(in) :: src
       integer(C_INT), intent(out) :: err
    end subroutine
    ```

**Arguments**
- `bas`: the basis shell to be created by copying
- `src`: the basis shell to be copied
- `err`: the error code (only Fortran)

**Return Value**
- the error code (only C)

### Count of basis functions ###

The number of the basis functions can be retrieved by using the following function or routine.

- C
    ```c
    int gtoint_basis_shell_get_count(gtoint_basis_shell_t bas);
    ```
- Fortran
    ```f90
    subroutine gtoint_basis_shell_get_count(bas, nb)
       type(C_PTR), intent(in) :: bas
       integer, intent(out) :: nb
    end subroutine
    ```

**Arguments**
- `bas`: the basis shell whose number of the basis functions is to be retrieved
- `nb`: the number of the basis functions (only Fortran)

**Return Value**
- the number of the basis functions (only C)

### Construction of scalar ECP ###

The scalar ECP is defined using the following function or routine.

- C
    ```c
    gtoint_error_t gtoint_ecp_shell_create(
        gtoint_ecp_shell_t *ecp, gtoint_integrator_t itg,
        int a, int r, int ng, const double *g, const double *c
    );
    ```
- Fortran
    ```f90
    subroutine gtoint_ecp_shell_create(ecp, itg, a, r, ng, g, c, err)
       type(C_PTR), intent(out) :: ecp
       type(C_PTR), intent(in) :: itg
       integer(C_INT), intent(in) :: a
       integer, intent(in) :: r
       integer, intent(in) :: ng
       real(8), intent(in) :: g(*)
       real(8), intent(in) :: c(*)
       integer(C_INT), intent(out) :: err
    end subroutine
    ```

**Arguments**
- `ecp`: the ECP shell to be created
- `itg`: the integrator
- `a`: the angular momentum quantum number; -1 means 'UL' potential
- `r`: the power restricted to the values 0, 1, and 2
- `ng`:the number of the contraction
- `g`: the array of the exponents
    - the array dimension is `ng`
- `c`: the array of the contract coefficients
    - the array dimension is `ng`
- `err`: the error code (only Fortran)

**Return Value**
- the error code (only C)

### Destruction of scalar ECP ###

The ECP must be destroyed by using the following function or routine when it is no more needed.

- C
    ```c
    void gtoint_ecp_shell_destroy(gtoint_ecp_shell_t ecp);
    ```
- Fortran
    ```f90
    subroutine gtoint_ecp_shell_destroy(ecp)
       type(C_PTR), intent(in) :: ecp
    end subroutine
    ```

**Argument**
- `ecp`: the ECP shell to be destroyed; nothing is done if `GTOINT_NULL` (in C) or `c_null_ptr` (in Fortran)

### Copy of scalar ECP ###

The scalar ECP can be copied using the following function or routine.

- C
    ```c
    gtoint_error_t gtoint_ecp_shell_copy(gtoint_ecp_shell_t *ecp, gtoint_ecp_shell_t src);
    ```
- Fortran
    ```f90
    subroutine gtoint_ecp_shell_copy(ecp, src, err)
       type(C_PTR), intent(out) :: ecp
       type(C_PTR), intent(in) :: src
       integer(C_INT), intent(out) :: err
    end subroutine
    ```

**Arguments**
- `ecp`: the ECP shell to be created by copying
- `src`: the ECP shell to be copied
- `err`: the error code (only Fortran)

**Return Value**
- the error code (only C)

### Computation of integrals ###

The several kinds of integrals can be computed by using the following functions or routines.

- C
    ```c
    gtoint_error_t gtoint_compute_overlap_integrals(
        gtoint_integrator_t itg,
        const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
        const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
        int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
        double *out
    );

    gtoint_error_t gtoint_compute_kinetic_energy_integrals(
        gtoint_integrator_t itg,
        const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
        const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
        int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
        double *out
    );

    gtoint_error_t gtoint_compute_nuclear_attraction_integrals(
        gtoint_integrator_t itg,
        const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
        const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
        const gtoint_double3_t *pc,
        int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *dc,
        double *out
    );

    gtoint_error_t gtoint_compute_multipole_moment_integrals(
        gtoint_integrator_t itg,
        const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
        const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
        const gtoint_double3_t *pm, int nam, const gtoint_int3_t *am,
        int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1,
        double *out
    );

    gtoint_error_t gtoint_compute_electron_repulsion_integrals(
        gtoint_integrator_t itg,
        const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
        const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
        const gtoint_double3_t *p2, gtoint_basis_shell_t bas2,
        const gtoint_double3_t *p3, gtoint_basis_shell_t bas3,
        int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *d2, const gtoint_int3_t *d3,
        double *out
    );

    gtoint_error_t gtoint_compute_ecp_integrals(
        gtoint_integrator_t itg,
        const gtoint_double3_t *p0, gtoint_basis_shell_t bas0,
        const gtoint_double3_t *p1, gtoint_basis_shell_t bas1,
        const gtoint_double3_t *pc, gtoint_ecp_shell_t ecp,
        int nd, const gtoint_int3_t *d0, const gtoint_int3_t *d1, const gtoint_int3_t *dc,
        double *out
    );
    ```
- Fortran
    ```f90
    subroutine gtoint_compute_overlap_integrals &
          & (itg, p0, bas0, p1, bas1, nd, d0, d1, out, err)
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: p0(3)
       type(C_PTR), intent(in) :: bas0
       real(8), intent(in) :: p1(3)
       type(C_PTR), intent(in) :: bas1
       integer, intent(in) :: nd
       integer(C_INT), intent(in) :: d0(3,nd)
       integer(C_INT), intent(in) :: d1(3,nd)
       real(8), intent(out) :: out(*)
       integer(C_INT), intent(out) :: err
    end subroutine

    subroutine gtoint_compute_kinetic_energy_integrals &
          & (itg, p0, bas0, p1, bas1, nd, d0, d1, out, err)
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: p0(3)
       type(C_PTR), intent(in) :: bas0
       real(8), intent(in) :: p1(3)
       type(C_PTR), intent(in) :: bas1
       integer, intent(in) :: nd
       integer(C_INT), intent(in) :: d0(3,nd)
       integer(C_INT), intent(in) :: d1(3,nd)
       real(8), intent(out) :: out(*)
       integer(C_INT), intent(out) :: err
    end subroutine

    subroutine gtoint_compute_nuclear_attraction_integrals &
          & (itg, p0, bas0, p1, bas1, pc, nd, d0, d1, dc, out, err)
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: p0(3)
       type(C_PTR), intent(in) :: bas0
       real(8), intent(in) :: p1(3)
       type(C_PTR), intent(in) :: bas1
       real(8), intent(in) :: pc(3)
       integer, intent(in) :: nd
       integer(C_INT), intent(in) :: d0(3,nd)
       integer(C_INT), intent(in) :: d1(3,nd)
       integer(C_INT), intent(in) :: dc(3,nd)
       real(8), intent(out) :: out(*)
       integer(C_INT), intent(out) :: err
    end subroutine

    subroutine gtoint_compute_multipole_moment_integrals &
          & (itg, p0, bas0, p1, bas1, pm, nam, am, nd, d0, d1, out, err)
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: p0(3)
       type(C_PTR), intent(in) :: bas0
       real(8), intent(in) :: p1(3)
       type(C_PTR), intent(in) :: bas1
       real(8), intent(in) :: pm(3)
       integer, intent(in) :: nam
       integer(C_INT), intent(in) :: am(3,nam)
       integer, intent(in) :: nd
       integer(C_INT), intent(in) :: d0(3,nd)
       integer(C_INT), intent(in) :: d1(3,nd)
       real(8), intent(out) :: out(*)
       integer(C_INT), intent(out) :: err
    end subroutine

    subroutine gtoint_compute_electron_repulsion_integrals &
          & (itg, p0, bas0, p1, bas1, p2, bas2, p3, bas3, nd, d0, d1, d2, d3, out, err)
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: p0(3)
       type(C_PTR), intent(in) :: bas0
       real(8), intent(in) :: p1(3)
       type(C_PTR), intent(in) :: bas1
       real(8), intent(in) :: p2(3)
       type(C_PTR), intent(in) :: bas2
       real(8), intent(in) :: p3(3)
       type(C_PTR), intent(in) :: bas3
       integer, intent(in) :: nd
       integer(C_INT), intent(in) :: d0(3,nd)
       integer(C_INT), intent(in) :: d1(3,nd)
       integer(C_INT), intent(in) :: d2(3,nd)
       integer(C_INT), intent(in) :: d3(3,nd)
       real(8), intent(out) :: out(*)
       integer(C_INT), intent(out) :: err
    end subroutine

    subroutine gtoint_compute_ecp_integrals &
          & (itg, p0, bas0, p1, bas1, pc, ecp, nd, d0, d1, dc, out, err)
       type(C_PTR), intent(in) :: itg
       real(8), intent(in) :: p0(3)
       type(C_PTR), intent(in) :: bas0
       real(8), intent(in) :: p1(3)
       type(C_PTR), intent(in) :: bas1
       real(8), intent(in) :: pc(3)
       type(C_PTR), intent(in) :: ecp
       integer, intent(in) :: nd
       integer(C_INT), intent(in) :: d0(3,nd)
       integer(C_INT), intent(in) :: d1(3,nd)
       integer(C_INT), intent(in) :: dc(3,nd)
       real(8), intent(out) :: out(*)
       integer(C_INT), intent(out) :: err
    end subroutine
    ```

**Arguments**
- `itg`: the integrator
- `bas0`, `bas1`, `bas2`, `bas3`: the basis shells
- `ecp`: the ECP shell
- `p0`, `p1`, `p2`, `p3`: the coordinates of the basis shells centers (in Bohr)
- `pc`: the coordinates of the charge center or the ECP center (in Bohr)
- `pm`: the coordinates of the multipole center
- `nd`: the number of the derivatives
- `d0`, `d1`, `d2`, `d3`: the array of the coordinate derivative orders of the basis shells centers
- `dc`: the array of the coordinate derivative orders of the charge center or the ECP center
- `nam`: the number of the multipole moments
- `am`: the array of the multipole moments
- `out`: the array of the integral values to be computed
    - as for `gtoint_compute_electron_repulsion_integrals()`, the array dimensions are `[nd][n3][n2][n1][n0]` in C-style notation or `(n0,n1,n2,n3,nd)` in Fortran-style notation
    - as for `gtoint_compute_multipole_moment_integrals()`, the array dimensions are `[nd][nam][n1][n0]` in C-style notation or `(n0,n1,nam,nd)` in Fortran-style notation
    - as for the others, the array dimensions are `[nd][n1][n0]` in C-style notation or `(n0,n1,nd)` in Fortran-style notation
    - here, `n0` is the number of the basis functions in the basis shell `bas0`, `n1` is the number of the basis functions in the basis shell `bas1`, and so on.
- `err`: the error code (only Fortran)

**Return Value**
- the error code (only C)

**Notice**
- The Cartesian basis functions of the same angular momentum quantum number are ordered alphabetically according to their name (ex. dxx, dxy, dxz, dyy, dyz, dzz).
  The spherical basis functions of the same angular momentum quantum number are ordered numerically according to their numbering (ex. d-1, d-2, d0, d+1, d+2).
- As for `gtoint_compute_nuclear_attraction_integrals()`, each integral value have to be multiplied by the negated charge of the nucleus at the position `pc`.
- As for `gtoint_compute_ecp_integrals()`, there is a limitation that the sum of the derivative orders of each coordinate component (x, y, or z) among 3 centers must be less than or equal to 2.
  The error code `GTOINT_ERROR_UNSUPPORTED` is returned if violating the limitation.

## Example ##

If we have the following basis set (written in Gaussian format):
```
Fe     0
SP   4   1.00
     70.28000                -.002611               -.007940
      6.061000               -.692435               -.290151
      4.134000               0.362530               0.591028
      1.421000               1.140645               0.719448
SP   2   1.00
      1.978000               -.098172               -.033731
      0.121300               1.026957               1.004462
SP   1   1.00
      0.512100               1.000000               1.000000
SP   1   1.00
      0.041000               1.000000               1.000000
D    4   1.00
     47.10000                0.026608
     13.12000                0.152010
      4.478000               0.413827
      1.581000               0.605542
D    1   1.00
      0.510000               1.000000
D    1   1.00
      0.138200               1.000000
****

FE     0
FE-ECP     2     10
d potential
  1
1     17.6191000             -3.8942300
s-d potential
  3
0      2.3315000              3.8170600
2      6.0836500            172.0534900
2      5.2665900           -144.7005600
p-d potential
  2
0     49.4045600              4.1373700
2     11.4018300             82.7369600
```
we can construct the basis shells and the ECP shells by the code shown below:
- C
    ```c
    gtoint_error_t e = GTOINT_ERROR_OK;
    gtoint_integrator_t itg = 0;
    gtoint_basis_shell_t bas[7] = { 0 };
    gtoint_ecp_shell_t ecp[5] = { 0 };

    /* Creation of Integrator */
    if ((e = gtoint_integrator_create(&itg)) != GTOINT_ERROR_OK) {
        printf("ERROR: gtoint_integrator_create\n");
        exit(1);
    }

    /* Construction of Basis Shells */
    {
        static const int a[2] = { 0, 1 }; /* s, p */
        static const double g[4] = {
            70.28000, 6.061000, 4.134000, 1.421000
        };
        static const double c[8] = {
            -0.002611, -0.692435, 0.362530, 1.140645,
            -0.007940, -0.290151, 0.591028, 0.719448
        };
        if ((e = gtoint_basis_shell_create(&(bas[0]), itg, 1.0, 2, a, 4, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_basis_shell_create\n");
            exit(1);
        }
    }
    {
        static const int a[2] = { 0, 1 }; /* s, p */
        static const double g[2] = {
            1.978000, 0.121300
        };
        static const double c[4] = {
            -0.098172, 1.026957,
            -0.033731, 1.004462
        };
        if ((e = gtoint_basis_shell_create(&(bas[1]), itg, 1.0, 2, a, 2, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_basis_shell_create\n");
            exit(1);
        }
    }
    {
        static const int a[2] = { 0, 1 }; /* s, p */
        static const double g[1] = {
            0.512100
        };
        static const double c[2] = {
            1.000000,
            1.000000
        };
        if ((e = gtoint_basis_shell_create(&(bas[2]), itg, 1.0, 2, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_basis_shell_create\n");
            exit(1);
        }
    }
    {
        static const int a[2] = { 0, 1 }; /* s, p */
        static const double g[1] = {
            0.041000
        };
        static const double c[2] = {
            1.000000,
            1.000000
        };
        if ((e = gtoint_basis_shell_create(&(bas[3]), itg, 1.0, 2, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_basis_shell_create\n");
            exit(1);
        }
    }
    {
        static const int a[1] = { ~2 }; /* d (spherical) */
        static const double g[4] = {
            47.10000, 13.12000, 4.478000, 1.581000
        };
        static const double c[4] = {
            0.026608, 0.152010, 0.413827, 0.605542
        };
        if ((e = gtoint_basis_shell_create(&(bas[4]), itg, 1.0, 1, a, 4, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_basis_shell_create\n");
            exit(1);
        }
    }
    {
        static const int a[1] = { ~2 }; /* d (spherical) */
        static const double g[1] = {
            0.510000
        };
        static const double c[1] = {
            1.000000
        };
        if ((e = gtoint_basis_shell_create(&(bas[5]), itg, 1.0, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_basis_shell_create\n");
            exit(1);
        }
    }
    {
        static const int a[1] = { ~2 }; /* d (spherical) */
        static const double g[1] = {
            0.138200
        };
        static const double c[1] = {
            1.000000
        };
        if ((e = gtoint_basis_shell_create(&(bas[6]), itg, 1.0, 1, a, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_basis_shell_create\n");
            exit(1);
        }
    }

    /* Construction of ECP Shells */
    {
        static const double g[1] = {
            17.6191000
        };
        static const double c[1] = {
            -3.8942300
        };
        if ((e = gtoint_ecp_shell_create(&(ecp[0]), itg, -1, 1, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_ecp_shell_create\n");
            exit(1);
        }
    }
    {
        static const double g[1] = {
            2.3315000
        };
        static const double c[1] = {
            3.8170600
        };
        if ((e = gtoint_ecp_shell_create(&(ecp[1]), itg, 0, 0, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_ecp_shell_create\n");
            exit(1);
        }
    }
    {
        static const double g[2] = {
            6.0836500, 5.2665900
        };
        static const double c[2] = {
            172.0534900, -144.7005600
        };
        if ((e = gtoint_ecp_shell_create(&(ecp[2]), itg, 0, 2, 2, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_ecp_shell_create\n");
            exit(1);
        }
    }
    {
        static const double g[1] = {
            49.4045600
        };
        static const double c[1] = {
            4.1373700
        };
        if ((e = gtoint_ecp_shell_create(&(ecp[3]), itg, 1, 0, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_ecp_shell_create\n");
            exit(1);
        }
    }
    {
        static const double g[1] = {
            11.4018300
        };
        static const double c[1] = {
            82.7369600
        };
        if ((e = gtoint_ecp_shell_create(&(ecp[4]), itg, 1, 2, 1, g, c)) != GTOINT_ERROR_OK) {
            printf("ERROR: gtoint_ecp_shell_create\n");
            exit(1);
        }
    }
    ```
- Fortran
    ```f90
    use gtoint
    implicit none
    type(C_PTR) :: itg, bas(7), ecp(5)
    integer(C_INT) :: a(2)
    real(8) :: g(4), c(8)
    integer :: err

    itg = c_null_ptr
    bas = c_null_ptr
    ecp = c_null_ptr

    ! Creation of Integrator

    call gtoint_integrator_create(itg, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_integrator_create"
       stop 1
    end if

    ! Construction of Basis Shells

    a(1) = 0 ! s
    a(2) = 1 ! p
    g(1) = 70.28000
    g(2) = 6.061000
    g(3) = 4.134000
    g(4) = 1.421000
    c(1) = -0.002611
    c(2) = -0.692435
    c(3) =  0.362530
    c(4) =  1.140645
    c(5) = -0.007940
    c(6) = -0.290151
    c(7) =  0.591028
    c(8) =  0.719448
    call gtoint_basis_shell_create(bas(1), itg, 1.0d0, 2, a, 4, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_basis_shell_create"
       stop 1
    end if

    a(1) = 0 ! s
    a(2) = 1 ! p
    g(1) = 1.978000
    g(2) = 0.121300
    c(1) = -0.098172
    c(2) =  1.026957
    c(3) = -0.033731
    c(4) =  1.004462
    call gtoint_basis_shell_create(bas(2), itg, 1.0d0, 2, a, 2, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_basis_shell_create"
       stop 1
    end if

    a(1) = 0 ! s
    a(2) = 1 ! p
    g(1) = 0.512100
    c(1) = 1.000000
    c(2) = 1.000000
    call gtoint_basis_shell_create(bas(3), itg, 1.0d0, 2, a, 1, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_basis_shell_create"
       stop 1
    end if

    a(1) = 0 ! s
    a(2) = 1 ! p
    g(1) = 0.041000
    c(1) = 1.000000
    c(2) = 1.000000
    call gtoint_basis_shell_create(bas(4), itg, 1.0d0, 2, a, 1, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_basis_shell_create"
       stop 1
    end if

    a(1) = not(2) ! d (spherical)
    g(1) = 47.10000
    g(2) = 13.12000
    g(3) = 4.478000
    g(4) = 1.581000
    c(1) = 0.026608
    c(2) = 0.152010
    c(3) = 0.413827
    c(4) = 0.605542
    call gtoint_basis_shell_create(bas(5), itg, 1.0d0, 1, a, 4, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_basis_shell_create"
       stop 1
    end if

    a(1) = not(2) ! d (spherical)
    g(1) = 0.510000
    c(1) = 1.000000
    call gtoint_basis_shell_create(bas(6), itg, 1.0d0, 1, a, 1, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_basis_shell_create"
       stop 1
    end if

    a(1) = not(2) ! d (spherical)
    g(1) = 0.138200
    c(1) = 1.000000
    call gtoint_basis_shell_create(bas(7), itg, 1.0d0, 1, a, 1, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_basis_shell_create"
       stop 1
    end if

    ! Construction of ECP Shells

    g(1) = 17.6191000
    c(1) = -3.8942300
    call gtoint_ecp_shell_create(ecp(1), itg, -1, 1, 1, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_ecp_shell_create"
       stop 1
    end if

    g(1) = 2.3315000
    c(1) = 3.8170600
    call gtoint_ecp_shell_create(ecp(2), itg, 0, 0, 1, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_ecp_shell_create"
       stop 1
    end if

    g(1) = 6.0836500
    g(2) = 5.2665900
    c(1) =  172.0534900
    c(2) = -144.7005600
    call gtoint_ecp_shell_create(ecp(3), itg, 0, 2, 2, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_ecp_shell_create"
       stop 1
    end if

    g(1) = 49.4045600
    c(1) =  4.1373700
    call gtoint_ecp_shell_create(ecp(4), itg, 1, 0, 1, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_ecp_shell_create"
       stop 1
    end if

    g(1) = 11.4018300
    c(1) = 82.7369600
    call gtoint_ecp_shell_create(ecp(5), itg, 1, 2, 1, g, c, err)
    if (err /= GTOINT_ERROR_OK) then
       print *, "ERROR: gtoint_ecp_shell_create"
       stop 1
    end if
    ```

Then, we can compute several types of integrals among the basis shells.

We must destroy the basis shells and the ECP shells after their use by the code shown below:
- C
    ```c
    gtoint_ecp_shell_destroy(ecp[4]);
    gtoint_ecp_shell_destroy(ecp[3]);
    gtoint_ecp_shell_destroy(ecp[2]);
    gtoint_ecp_shell_destroy(ecp[1]);
    gtoint_ecp_shell_destroy(ecp[0]);
    gtoint_basis_shell_destroy(bas[6]);
    gtoint_basis_shell_destroy(bas[5]);
    gtoint_basis_shell_destroy(bas[4]);
    gtoint_basis_shell_destroy(bas[3]);
    gtoint_basis_shell_destroy(bas[2]);
    gtoint_basis_shell_destroy(bas[1]);
    gtoint_basis_shell_destroy(bas[0]);
    gtoint_integrator_destroy(itg);
    ```
- Fortran
    ```f90
    call gtoint_ecp_shell_destroy(ecp(5));
    call gtoint_ecp_shell_destroy(ecp(4));
    call gtoint_ecp_shell_destroy(ecp(3));
    call gtoint_ecp_shell_destroy(ecp(2));
    call gtoint_ecp_shell_destroy(ecp(1));
    call gtoint_basis_shell_destroy(bas(7));
    call gtoint_basis_shell_destroy(bas(6));
    call gtoint_basis_shell_destroy(bas(5));
    call gtoint_basis_shell_destroy(bas(4));
    call gtoint_basis_shell_destroy(bas(3));
    call gtoint_basis_shell_destroy(bas(2));
    call gtoint_basis_shell_destroy(bas(1));
    call gtoint_integrator_destroy(itg);
    ```
