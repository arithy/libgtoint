! LibGtoint: an analytical GTO integral library for C and Fortran.
!
! Copyright (c) 2020-2021 Arihiro Yoshida. All rights reserved.
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.

module gtoint
   use iso_c_binding

   implicit none

   integer(C_INT), parameter :: GTOINT_ERROR_OK          = 0
   integer(C_INT), parameter :: GTOINT_ERROR_ARGUMENT    = 1
   integer(C_INT), parameter :: GTOINT_ERROR_MEMORY      = 2
   integer(C_INT), parameter :: GTOINT_ERROR_UNSUPPORTED = -2
   integer(C_INT), parameter :: GTOINT_ERROR_INTERNAL    = -1

contains

   subroutine gtoint_integrator_create(itg, err)
      type(C_PTR), intent(out) :: itg
      integer(C_INT), intent(out) :: err
      interface
         function c_gtoint_integrator_create(itg) &
               & bind(c, name = 'gtoint_integrator_create')
            import C_PTR, C_INT
            type(C_PTR), intent(out) :: itg
            integer(C_INT) :: c_gtoint_integrator_create
         end function
      end interface
      err = c_gtoint_integrator_create(itg)
   end subroutine

   subroutine gtoint_integrator_destroy(itg)
      type(C_PTR), intent(in) :: itg
      interface
         subroutine c_gtoint_integrator_destroy(itg) &
               & bind(c, name = 'gtoint_integrator_destroy')
            import C_PTR
            type(C_PTR), value :: itg
         end subroutine
      end interface
      call c_gtoint_integrator_destroy(itg)
   end subroutine

   subroutine gtoint_integrator_cleanup_memory(itg)
      type(C_PTR), intent(in) :: itg
      interface
         subroutine c_gtoint_integrator_cleanup_memory(itg) &
               & bind(c, name = 'gtoint_integrator_cleanup_memory')
            import C_PTR
            type(C_PTR), value :: itg
         end subroutine
      end interface
      call c_gtoint_integrator_cleanup_memory(itg)
   end subroutine

   subroutine gtoint_integrator_set_error_tolerance(itg, tol)
      type(C_PTR), intent(in) :: itg
      real(8), intent(in) :: tol
      interface
         subroutine c_gtoint_integrator_set_error_tolerance(itg, tol) &
               & bind(c, name = 'gtoint_integrator_set_error_tolerance')
            import C_PTR, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE), value :: tol
         end subroutine
      end interface
      call c_gtoint_integrator_set_error_tolerance(itg, tol)
   end subroutine

   subroutine gtoint_integrator_get_error_tolerance(itg, tol)
      type(C_PTR), intent(in) :: itg
      real(8), intent(out) :: tol
      interface
         function c_gtoint_integrator_get_error_tolerance(itg) &
               & bind(c, name = 'gtoint_integrator_get_error_tolerance')
            import C_PTR, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE) :: c_gtoint_integrator_get_error_tolerance
         end function
      end interface
      tol = c_gtoint_integrator_get_error_tolerance(itg)
   end subroutine

   subroutine gtoint_integrator_set_cutoff(itg, cut)
      type(C_PTR), intent(in) :: itg
      real(8), intent(in) :: cut
      interface
         subroutine c_gtoint_integrator_set_cutoff(itg, cut) &
               & bind(c, name = 'gtoint_integrator_set_cutoff')
            import C_PTR, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE), value :: cut
         end subroutine
      end interface
      call c_gtoint_integrator_set_cutoff(itg, cut)
   end subroutine

   subroutine gtoint_integrator_get_cutoff(itg, cut)
      type(C_PTR), intent(in) :: itg
      real(8), intent(out) :: cut
      interface
         function c_gtoint_integrator_get_cutoff(itg) &
               & bind(c, name = 'gtoint_integrator_get_cutoff')
            import C_PTR, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE) :: c_gtoint_integrator_get_cutoff
         end function
      end interface
      cut = c_gtoint_integrator_get_cutoff(itg)
   end subroutine

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
      interface
         function c_gtoint_basis_shell_create(bas, itg, sf, na, a, ng, g, c) &
               & bind(c, name = 'gtoint_basis_shell_create')
            import C_PTR, C_INT, C_DOUBLE
            type(C_PTR), intent(out) :: bas
            type(C_PTR), value :: itg
            real(C_DOUBLE), value :: sf
            integer(C_INT), value :: na
            integer(C_INT), intent(in) :: a(*)
            integer(C_INT), value :: ng
            real(C_DOUBLE), intent(in) :: g(*)
            real(C_DOUBLE), intent(in) :: c(*)
            integer(C_INT) :: c_gtoint_basis_shell_create
         end function
      end interface
      err = c_gtoint_basis_shell_create(bas, itg, sf, na, a, ng, g, c)
   end subroutine

   subroutine gtoint_basis_shell_destroy(bas)
      type(C_PTR), intent(in) :: bas
      interface
         subroutine c_gtoint_basis_shell_destroy(bas) &
               & bind(c, name = 'gtoint_basis_shell_destroy')
            import C_PTR
            type(C_PTR), value :: bas
         end subroutine
      end interface
      call c_gtoint_basis_shell_destroy(bas)
   end subroutine

   subroutine gtoint_basis_shell_get_count(bas, nb)
      type(C_PTR), intent(in) :: bas
      integer, intent(out) :: nb
      interface
         function c_gtoint_basis_shell_get_count(bas) &
               & bind(c, name = 'gtoint_basis_shell_get_count')
            import C_PTR, C_INT
            type(C_PTR), value :: bas
            integer(C_INT) :: c_gtoint_basis_shell_get_count
         end function
      end interface
      nb = c_gtoint_basis_shell_get_count(bas)
   end subroutine

   subroutine gtoint_ecp_shell_create(ecp, itg, a, r, ng, g, c, err)
      type(C_PTR), intent(out) :: ecp
      type(C_PTR), intent(in) :: itg
      integer(C_INT), intent(in) :: a
      integer, intent(in) :: r
      integer, intent(in) :: ng
      real(8), intent(in) :: g(*)
      real(8), intent(in) :: c(*)
      integer(C_INT), intent(out) :: err
      interface
         function c_gtoint_ecp_shell_create(ecp, itg, a, r, ng, g, c) &
               & bind(c, name = 'gtoint_ecp_shell_create')
            import C_PTR, C_INT, C_DOUBLE
            type(C_PTR), intent(out) :: ecp
            type(C_PTR), value :: itg
            integer(C_INT), value :: a
            integer(C_INT), value :: r
            integer(C_INT), value :: ng
            real(C_DOUBLE), intent(in) :: g(*)
            real(C_DOUBLE), intent(in) :: c(*)
            integer(C_INT) :: c_gtoint_ecp_shell_create
         end function
      end interface
      err = c_gtoint_ecp_shell_create(ecp, itg, a, r, ng, g, c)
   end subroutine

   subroutine gtoint_ecp_shell_destroy(ecp)
      type(C_PTR), intent(in) :: ecp
      interface
         subroutine c_gtoint_ecp_shell_destroy(ecp) &
               & bind(c, name = 'gtoint_ecp_shell_destroy')
            import C_PTR
            type(C_PTR), value :: ecp
         end subroutine
      end interface
      call c_gtoint_ecp_shell_destroy(ecp)
   end subroutine

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
      interface
         function c_gtoint_compute_overlap_integrals &
               & (itg, p0, bas0, p1, bas1, nd, d0, d1, out) &
               & bind(c, name = 'gtoint_compute_overlap_integrals')
            import C_PTR, C_INT, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE), intent(in) :: p0(*)
            type(C_PTR), value :: bas0
            real(C_DOUBLE), intent(in) :: p1(*)
            type(C_PTR), value :: bas1
            integer(C_INT), value :: nd
            integer(C_INT), intent(in) :: d0(*)
            integer(C_INT), intent(in) :: d1(*)
            real(C_DOUBLE), intent(out) :: out(*)
            integer(C_INT) :: c_gtoint_compute_overlap_integrals
         end function
      end interface
      err = c_gtoint_compute_overlap_integrals &
            & (itg, p0, bas0, p1, bas1, nd, d0, d1, out)
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
      interface
         function c_gtoint_compute_kinetic_energy_integrals &
               & (itg, p0, bas0, p1, bas1, nd, d0, d1, out) &
               & bind(c, name = 'gtoint_compute_kinetic_energy_integrals')
            import C_PTR, C_INT, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE), intent(in) :: p0(*)
            type(C_PTR), value :: bas0
            real(C_DOUBLE), intent(in) :: p1(*)
            type(C_PTR), value :: bas1
            integer(C_INT), value :: nd
            integer(C_INT), intent(in) :: d0(*)
            integer(C_INT), intent(in) :: d1(*)
            real(C_DOUBLE), intent(out) :: out(*)
            integer(C_INT) :: c_gtoint_compute_kinetic_energy_integrals
         end function
      end interface
      err = c_gtoint_compute_kinetic_energy_integrals &
            & (itg, p0, bas0, p1, bas1, nd, d0, d1, out)
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
      interface
         function c_gtoint_compute_nuclear_attraction_integrals &
               & (itg, p0, bas0, p1, bas1, pc, nd, d0, d1, dc, out) &
               & bind(c, name = 'gtoint_compute_nuclear_attraction_integrals')
            import C_PTR, C_INT, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE), intent(in) :: p0(*)
            type(C_PTR), value :: bas0
            real(C_DOUBLE), intent(in) :: p1(*)
            type(C_PTR), value :: bas1
            real(C_DOUBLE), intent(in) :: pc(*)
            integer(C_INT), value :: nd
            integer(C_INT), intent(in) :: d0(*)
            integer(C_INT), intent(in) :: d1(*)
            integer(C_INT), intent(in) :: dc(*)
            real(C_DOUBLE), intent(out) :: out(*)
            integer(C_INT) :: c_gtoint_compute_nuclear_attraction_integrals
         end function
      end interface
      err = c_gtoint_compute_nuclear_attraction_integrals &
            & (itg, p0, bas0, p1, bas1, pc, nd, d0, d1, dc, out)
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
      interface
         function c_gtoint_compute_multipole_moment_integrals &
               & (itg, p0, bas0, p1, bas1, pm, nam, am, nd, d0, d1, out) &
               & bind(c, name = 'gtoint_compute_multipole_moment_integrals')
            import C_PTR, C_INT, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE), intent(in) :: p0(*)
            type(C_PTR), value :: bas0
            real(C_DOUBLE), intent(in) :: p1(*)
            type(C_PTR), value :: bas1
            real(C_DOUBLE), intent(in) :: pm(*)
            integer(C_INT), value :: nam
            integer(C_INT), intent(in) :: am(*)
            integer(C_INT), value :: nd
            integer(C_INT), intent(in) :: d0(*)
            integer(C_INT), intent(in) :: d1(*)
            real(C_DOUBLE), intent(out) :: out(*)
            integer(C_INT) :: c_gtoint_compute_multipole_moment_integrals
         end function
      end interface
      err = c_gtoint_compute_multipole_moment_integrals &
            & (itg, p0, bas0, p1, bas1, pm, nam, am, nd, d0, d1, out)
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
      interface
         function c_gtoint_compute_electron_repulsion_integrals &
               & (itg, p0, bas0, p1, bas1, p2, bas2, p3, bas3, nd, d0, d1, d2, d3, out) &
               & bind(c, name = 'gtoint_compute_electron_repulsion_integrals')
            import C_PTR, C_INT, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE), intent(in) :: p0(*)
            type(C_PTR), value :: bas0
            real(C_DOUBLE), intent(in) :: p1(*)
            type(C_PTR), value :: bas1
            real(C_DOUBLE), intent(in) :: p2(*)
            type(C_PTR), value :: bas2
            real(C_DOUBLE), intent(in) :: p3(*)
            type(C_PTR), value :: bas3
            integer(C_INT), value :: nd
            integer(C_INT), intent(in) :: d0(*)
            integer(C_INT), intent(in) :: d1(*)
            integer(C_INT), intent(in) :: d2(*)
            integer(C_INT), intent(in) :: d3(*)
            real(C_DOUBLE), intent(out) :: out(*)
            integer(C_INT) :: c_gtoint_compute_electron_repulsion_integrals
         end function
      end interface
      err = c_gtoint_compute_electron_repulsion_integrals &
            & (itg, p0, bas0, p1, bas1, p2, bas2, p3, bas3, nd, d0, d1, d2, d3, out)
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
      interface
         function c_gtoint_compute_ecp_integrals &
               & (itg, p0, bas0, p1, bas1, pc, ecp, nd, d0, d1, dc, out) &
               & bind(c, name = 'gtoint_compute_ecp_integrals')
            import C_PTR, C_INT, C_DOUBLE
            type(C_PTR), value :: itg
            real(C_DOUBLE), intent(in) :: p0(*)
            type(C_PTR), value :: bas0
            real(C_DOUBLE), intent(in) :: p1(*)
            type(C_PTR), value :: bas1
            real(C_DOUBLE), intent(in) :: pc(*)
            type(C_PTR), value :: ecp
            integer(C_INT), value :: nd
            integer(C_INT), intent(in) :: d0(*)
            integer(C_INT), intent(in) :: d1(*)
            integer(C_INT), intent(in) :: dc(*)
            real(C_DOUBLE), intent(out) :: out(*)
            integer(C_INT) :: c_gtoint_compute_ecp_integrals
         end function
      end interface
      err = c_gtoint_compute_ecp_integrals &
            & (itg, p0, bas0, p1, bas1, pc, ecp, nd, d0, d1, dc, out)
   end subroutine

end module gtoint
