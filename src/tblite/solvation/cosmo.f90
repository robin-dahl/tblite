! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

#ifndef TBLITE_HAS_LIBCINT
#define TBLITE_HAS_LIBCINT 0
#endif

!> Standard COSMO/CPCM electrostatics using a moist cavity.
module tblite_solvation_cosmo
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   ! sVDW drop cavity:
   ! use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
   ! use moist_cavity_drop_lsf_svdw, only : moist_cavity_drop_lsf_svdw_type
   use moist_cavity_iswig, only : cavity_type_iswig, new_cavity_iswig
   use moist_model_component_pcm_solvers, only : solve_pcm_cholesky
   use moist_radii, only : radius_type, new_radii_custom_atoms
   use tblite_basis_type, only : basis_type
   use tblite_container_cache, only : container_cache
#if TBLITE_HAS_LIBCINT
   use tblite_integral_libcint, only : libcint_integral_type
#endif
   use tblite_scf_info, only : atom_resolved, orbital_resolved, scf_info
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: cosmo_input, cosmo_solvation, cosmo_cache, new_cosmo
   public :: cosmo_solvation_model

   type :: enum_cosmo_solvation_model
      integer :: cosmo = 1
      integer :: cpcm = 2
   end type enum_cosmo_solvation_model

   type(enum_cosmo_solvation_model), parameter :: cosmo_solvation_model = &
      & enum_cosmo_solvation_model()

   !> Input for the COSMO electrostatic term.
   type :: cosmo_input
      !> Conductor-like model variant.
      integer :: model = cosmo_solvation_model%cosmo
      !> Relative dielectric permittivity.
      real(wp) :: dielectric_const
      !> Optional element-resolved radii; indexed through structure%id.
      real(wp), allocatable :: rvdw(:)
      !> Scale factor applied to the cavity radii.
      real(wp) :: rscale = 1.0_wp
      !> Number of Lebedev points on each atomic sphere.
      integer :: nang = 110
      !> Include atom-resolved dipoles in the solute density.
      logical :: dipoles = .false.
      !> Include atom-resolved quadrupoles in the solute density.
      logical :: quadrupoles = .false.
      !> Use the complete AO density instead of atom-resolved multipoles.
      logical :: full_density = .false.
   end type cosmo_input

   interface cosmo_input
      module procedure :: create_cosmo_input
   end interface cosmo_input

   !> COSMO model with a multipolar representation of the solute density.
   type, extends(solvation_type) :: cosmo_solvation
      !> COSMO dielectric scaling, (epsilon-1)/(epsilon+1/2).
      real(wp) :: feps
      !> Relative dielectric permittivity.
      real(wp) :: dielectric_const
      !> Atom-resolved cavity radii.
      real(wp), allocatable :: rvdw(:)
      !> Number of Lebedev points per sphere.
      integer :: nang
      !> Include atom-resolved dipoles.
      logical :: dipoles = .false.
      !> Include atom-resolved quadrupoles.
      logical :: quadrupoles = .false.
      !> Use the complete AO density.
      logical :: full_density = .false.
      !> Gaussian AO basis needed for full-density coupling.
      type(basis_type) :: basis
   contains
      procedure :: update
      procedure :: variable_info
      procedure :: get_energy
      procedure :: get_potential
   end type cosmo_solvation

   !> Geometry-dependent COSMO data.
   type :: cosmo_cache
      ! sVDW drop cavity (previous implementation):
      ! type(cavity_type_drop) :: cavity
      type(cavity_type_iswig) :: cavity
      real(wp), allocatable :: amat(:, :)
      !> Potential kernels for monopoles, dipoles, and quadrupoles.
      real(wp), allocatable :: c0(:, :)
      real(wp), allocatable :: c1(:, :, :)
      real(wp), allocatable :: c2(:, :, :)
      real(wp), allocatable :: phi(:)
      real(wp), allocatable :: qsurf(:)
      !> Nuclear and AO-pair coupling to the Gaussian surface basis.
      real(wp), allocatable :: ucore(:, :)
      real(wp), allocatable :: bmat(:, :, :)
      logical :: ready = .false.
   end type cosmo_cache

contains

function create_cosmo_input(dielectric_const, model, rvdw, rscale, nang, &
      & dipoles, quadrupoles, full_density) result(self)
   real(wp), intent(in) :: dielectric_const
   integer, intent(in), optional :: model
   real(wp), intent(in), optional :: rvdw(:)
   real(wp), intent(in), optional :: rscale
   integer, intent(in), optional :: nang
   logical, intent(in), optional :: dipoles, quadrupoles
   logical, intent(in), optional :: full_density
   type(cosmo_input) :: self

   self%dielectric_const = dielectric_const
   if (present(model)) self%model = model
   if (present(rvdw)) self%rvdw = rvdw
   if (present(rscale)) self%rscale = rscale
   if (present(nang)) self%nang = nang
   if (present(dipoles)) self%dipoles = dipoles
   if (present(quadrupoles)) self%quadrupoles = quadrupoles
   if (present(full_density)) self%full_density = full_density
end function create_cosmo_input

!> Construct a COSMO model. The cavity itself is built during update.
subroutine new_cosmo(self, mol, input, error, basis)
   type(cosmo_solvation), intent(out) :: self
   type(structure_type), intent(in) :: mol
   type(cosmo_input), intent(in) :: input
   type(error_type), allocatable, intent(out) :: error
   type(basis_type), intent(in), optional :: basis

   integer :: iat

   if (input%dielectric_const <= 1.0_wp) then
      call fatal_error(error, "COSMO dielectric constant must be larger than one")
      return
   end if
   if (input%nang <= 0) then
      call fatal_error(error, "COSMO angular grid size must be positive")
      return
   end if

   self%dielectric_const = input%dielectric_const
   select case(input%model)
   case(cosmo_solvation_model%cosmo)
      self%label = "COSMO solvation model"
      self%feps = (input%dielectric_const - 1.0_wp) / &
         & (input%dielectric_const + 0.5_wp)
   case(cosmo_solvation_model%cpcm)
      self%label = "CPCM solvation model"
      self%feps = (input%dielectric_const - 1.0_wp) / input%dielectric_const
   case default
      call fatal_error(error, "Unknown COSMO/CPCM solvation model")
      return
   end select
   self%nang = input%nang
   self%dipoles = input%dipoles
   self%quadrupoles = input%quadrupoles
   self%full_density = input%full_density
   if (self%full_density) then
      if (.not.present(basis)) then
         call fatal_error(error, "Full-density COSMO requires a Gaussian AO basis")
         return
      end if
      self%basis = basis
   end if
   allocate(self%rvdw(mol%nat))
   if (allocated(input%rvdw)) then
      if (maxval(mol%id) > size(input%rvdw)) then
         call fatal_error(error, "COSMO radii do not cover all species")
         return
      end if
      self%rvdw = input%rscale * input%rvdw(mol%id)
   else
      do iat = 1, mol%nat
         self%rvdw(iat) = input%rscale * get_vdw_rad_cosmo(mol%num(mol%id(iat)))
      end do
   end if
end subroutine new_cosmo

!> Build the moist cavity, PCM matrix, and multipole-to-surface kernels.
subroutine update(self, mol, cache)
   class(cosmo_solvation), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache

   type(cosmo_cache), pointer :: ptr
   class(radius_type), allocatable :: radii
   ! sVDW drop cavity (previous implementation):
   ! type(moist_cavity_drop_lsf_svdw_type) :: svdw
   type(error_type), allocatable :: error
   integer :: igrid, iat
   real(wp) :: vec(3), r2, r1, r3, r5

   call taint(cache, ptr)
   ptr%ready = .false.

   call new_radii_custom_atoms(self%rvdw, radii, error)
   if (allocated(error)) return
   ! sVDW drop cavity (previous implementation):
   ! call svdw%new(blend_k=6.5_wp)
   ! call new_cavity_drop(ptr%cavity, verbose=0, nleb=self%nang, &
   !    & radius_model=radii, lsf_model=svdw, error=error)

   ! iSwiG cavity. The same custom atom radii and requested Lebedev grid are
   ! passed to moist; cut_a and cut_f retain moist's iSwiG defaults.
   call new_cavity_iswig(ptr%cavity, nleb=self%nang, radius_model=radii, &
      & error=error)
   if (allocated(error)) return
   call ptr%cavity%update(mol, error)
   if (allocated(error)) return

   if (allocated(ptr%amat)) deallocate(ptr%amat)
   if (allocated(ptr%c0)) deallocate(ptr%c0)
   if (allocated(ptr%c1)) deallocate(ptr%c1)
   if (allocated(ptr%c2)) deallocate(ptr%c2)
   if (allocated(ptr%phi)) deallocate(ptr%phi)
   if (allocated(ptr%qsurf)) deallocate(ptr%qsurf)
   allocate(ptr%amat(ptr%cavity%ngrid, ptr%cavity%ngrid))
   ! The standalone iSwiG cavity currently provides a collocation-style
   ! matrix which is not the Coulomb matrix of its Gaussian surface basis.
   ! Assemble the latter explicitly so that A, Ucore, and the AO-pair
   ! integrals all use the same normalized Gaussians.
   call build_iswig_amat(ptr%cavity%xyz, ptr%cavity%xi, ptr%cavity%f, &
      & ptr%amat)

   allocate(ptr%c0(ptr%cavity%ngrid, mol%nat))
   allocate(ptr%c1(3, ptr%cavity%ngrid, mol%nat))
   allocate(ptr%c2(6, ptr%cavity%ngrid, mol%nat))
   allocate(ptr%phi(ptr%cavity%ngrid), ptr%qsurf(ptr%cavity%ngrid))

   if (self%full_density) then
      if (allocated(ptr%ucore)) deallocate(ptr%ucore)
      if (allocated(ptr%bmat)) deallocate(ptr%bmat)
      allocate(ptr%ucore(ptr%cavity%ngrid, mol%nat))
      do iat = 1, mol%nat
         do igrid = 1, ptr%cavity%ngrid
            vec = ptr%cavity%xyz(:, igrid) - mol%xyz(:, iat)
            r3 = dot_product(vec, vec)
            if (r3 > epsilon(1.0_wp)) then
               r1 = sqrt(r3)
               ! Moist stores xi=sqrt(alpha) for the normalized Gaussian
               ! (alpha/pi)^(3/2)*exp(-alpha*r^2).
               ptr%ucore(igrid, iat) = erf(ptr%cavity%xi(igrid)*r1)/r1
            else
               ! Analytic r -> 0 limit of erf(xi*r)/r.
               ptr%ucore(igrid, iat) = 2.0_wp*ptr%cavity%xi(igrid)/sqrt(acos(-1.0_wp))
            end if
         end do
      end do
      call build_cosmo_density_integrals(mol, self%basis, ptr%cavity%xyz, &
         & ptr%cavity%xi, ptr%bmat, error)
      if (allocated(error)) return
   end if

   do iat = 1, mol%nat
      do igrid = 1, ptr%cavity%ngrid
         vec = ptr%cavity%xyz(:, igrid) - mol%xyz(:, iat)
         r2 = sum(vec**2)
         r1 = sqrt(r2)
         r3 = r1*r2
         r5 = r3*r2
         ptr%c0(igrid, iat) = 1.0_wp/r1
         ptr%c1(:, igrid, iat) = vec/r3
         ptr%c2(:, igrid, iat) = [vec(1)*vec(1), 2.0_wp*vec(1)*vec(2), &
            & vec(2)*vec(2), 2.0_wp*vec(1)*vec(3), &
            & 2.0_wp*vec(2)*vec(3), vec(3)*vec(3)]/r5
      end do
   end do
   ptr%ready = .true.
end subroutine update

!> Build A(i,j)=(g_i|g_j) for the normalized iSwiG surface Gaussians.
subroutine build_iswig_amat(xyz, xi, switch, amat)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: xi(:)
   real(wp), intent(in) :: switch(:)
   real(wp), intent(out) :: amat(:, :)

   integer :: igrid, jgrid
   real(wp) :: xiij, rij

   do igrid = 1, size(xi)
      amat(igrid, igrid) = sqrt(2.0_wp/acos(-1.0_wp))*xi(igrid) &
         & /switch(igrid)
      do jgrid = 1, igrid - 1
         xiij = xi(igrid)*xi(jgrid) &
            & /sqrt(xi(igrid)**2 + xi(jgrid)**2)
         rij = norm2(xyz(:, igrid) - xyz(:, jgrid))
         amat(igrid, jgrid) = erf(xiij*rij)/rij
         amat(jgrid, igrid) = amat(igrid, jgrid)
      end do
   end do
end subroutine build_iswig_amat

!> Form the solute potential and solve for apparent surface charges.
subroutine solve_surface(self, ptr, wfn)
   class(cosmo_solvation), intent(in) :: self
   type(cosmo_cache), intent(inout) :: ptr
   type(wavefunction_type), intent(in) :: wfn

   type(error_type), allocatable :: error
   integer :: iat, ic
   real(wp), allocatable :: density(:, :)

   if (self%full_density) then
      ptr%phi = matmul(ptr%ucore, wfn%n0at)
      density = sum(wfn%density, dim=3)
      do ic = 1, size(ptr%phi)
         ptr%phi(ic) = ptr%phi(ic) - sum(ptr%bmat(ic, :, :)*density)
      end do
   else
      ptr%phi = matmul(ptr%c0, wfn%qat(:, 1))
   end if
   do iat = 1, size(wfn%qat, 1)
      if (.not.self%full_density .and. self%dipoles) then
         do ic = 1, 3
            ptr%phi = ptr%phi + ptr%c1(ic, :, iat)*wfn%dpat(ic, iat, 1)
         end do
      end if
      if (.not.self%full_density .and. self%quadrupoles) then
         do ic = 1, 6
            ptr%phi = ptr%phi + ptr%c2(ic, :, iat)*wfn%qpat(ic, iat, 1)
         end do
      end if
   end do
   call solve_pcm_cholesky(ptr%amat, -self%feps*ptr%phi, ptr%qsurf, error)
   if (allocated(error)) ptr%ready = .false.
end subroutine solve_surface

subroutine get_energy(self, mol, cache, wfn, energies)
   class(cosmo_solvation), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(inout) :: energies(:)

   type(cosmo_cache), pointer :: ptr
   integer :: iat, ic
   real(wp) :: vat, vdp(3), vqp(6)

   call view(cache, ptr)
   if (.not.associated(ptr)) return
   if (.not.ptr%ready) return
   call solve_surface(self, ptr, wfn)
   if (.not.ptr%ready) return
   if (self%full_density) then
      energies(:) = energies(:) + 0.5_wp*dot_product(ptr%qsurf, ptr%phi)/real(mol%nat, wp)
   else
      do iat = 1, mol%nat
         vat = dot_product(ptr%c0(:, iat), ptr%qsurf)
         do ic = 1, 3
            vdp(ic) = dot_product(ptr%c1(ic, :, iat), ptr%qsurf)
         end do
         do ic = 1, 6
            vqp(ic) = dot_product(ptr%c2(ic, :, iat), ptr%qsurf)
         end do
         energies(iat) = energies(iat) + 0.5_wp*wfn%qat(iat, 1)*vat
         if (self%dipoles) energies(iat) = energies(iat) + &
            & 0.5_wp*dot_product(wfn%dpat(:, iat, 1), vdp)
         if (self%quadrupoles) energies(iat) = energies(iat) + &
            & 0.5_wp*dot_product(wfn%qpat(:, iat, 1), vqp)
      end do
   end if

end subroutine get_energy

subroutine get_potential(self, mol, cache, wfn, pot)
   class(cosmo_solvation), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   type(potential_type), intent(inout) :: pot

   type(cosmo_cache), pointer :: ptr
   integer :: iat, ic

   call view(cache, ptr)
   if (.not.associated(ptr)) return
   if (.not.ptr%ready) return
   call solve_surface(self, ptr, wfn)
   if (.not.ptr%ready) return
   if (self%full_density) then
      do ic = 1, size(ptr%qsurf)
         do iat = 1, size(pot%vmat, 3)
            pot%vmat(:, :, iat) = pot%vmat(:, :, iat) - &
               & ptr%qsurf(ic)*ptr%bmat(ic, :, :)
         end do
      end do
   else
      do iat = 1, mol%nat
         pot%vat(iat, 1) = pot%vat(iat, 1) + dot_product(ptr%c0(:, iat), ptr%qsurf)
         if (self%dipoles) then
            do ic = 1, 3
               pot%vdp(ic, iat, 1) = pot%vdp(ic, iat, 1) + &
                  & dot_product(ptr%c1(ic, :, iat), ptr%qsurf)
            end do
         end if
         if (self%quadrupoles) then
            do ic = 1, 6
               pot%vqp(ic, iat, 1) = pot%vqp(ic, iat, 1) + &
                  & dot_product(ptr%c2(ic, :, iat), ptr%qsurf)
            end do
         end if
      end do
   end if
end subroutine get_potential

pure function variable_info(self) result(info)
   class(cosmo_solvation), intent(in) :: self
   type(scf_info) :: info
   if (self%full_density) then
      info = scf_info(density=orbital_resolved)
   else
      info = scf_info(charge=atom_resolved, &
      & dipole=merge(atom_resolved, 0, self%dipoles), &
      & quadrupole=merge(atom_resolved, 0, self%quadrupoles))
   end if
end function variable_info

subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(cosmo_cache), pointer, intent(out) :: ptr
   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if
   if (.not.allocated(cache%raw)) allocate(cosmo_cache :: cache%raw)
   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(cosmo_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(cosmo_cache)
      ptr => target
   end select
end subroutine view

subroutine build_cosmo_density_integrals(mol, basis, xyz, xi, bmat, error)
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: basis
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in) :: xi(:)
   real(wp), allocatable, intent(out) :: bmat(:, :, :)
   type(error_type), allocatable, intent(out) :: error

#if TBLITE_HAS_LIBCINT
   type(libcint_integral_type) :: libcint

   if (size(xyz, 1) /= 3 .or. size(xyz, 2) /= size(xi)) then
      call fatal_error(error, "Invalid moist surface Gaussian data")
      return
   end if

   call libcint%initialize_integral(mol, basis)
   call libcint%surface_3c2e(mol, basis, xyz, xi, bmat)
#else
   call fatal_error(error, "Full-density COSMO requires libcint support")
#endif
end subroutine build_cosmo_density_integrals

end module tblite_solvation_cosmo
