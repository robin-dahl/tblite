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
   use moist_model_component_pcm_solvers, only : solve_pcm_iterative
   use moist_radii, only : radius_type, new_radii_custom_atoms
   use moist_radii_static, only : static_radius_type, new_gauss_radii
   use tblite_basis_type, only : basis_type
   use tblite_container_cache, only : container_cache
#if TBLITE_HAS_LIBCINT
   use tblite_integral_libcint, only : libcint_integral_type
#endif
   use tblite_scf_info, only : atom_resolved, not_used, orbital_resolved, scf_info
   use tblite_scf_potential, only : potential_type
   use tblite_solvation_type, only : solvation_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: cosmo_input, cosmo_solvation, cosmo_cache, new_cosmo
   public :: write_cpcm_file
   public :: get_gaussian_multipole_kernels
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
      !> Include atom-resolved monopoles in the solute density.
      logical :: monopoles = .true.
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
      !> Include atom-resolved monopoles.
      logical :: monopoles = .true.
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
      ! sVDW drop cavity:
      ! type(cavity_type_drop) :: cavity
      type(cavity_type_iswig) :: cavity
      real(wp), allocatable :: amat(:, :)
      !> Potential kernels for monopoles, dipoles, and quadrupoles.
      real(wp), allocatable :: c0(:, :)
      real(wp), allocatable :: c1(:, :)
      real(wp), allocatable :: c2(:, :)
      real(wp), allocatable :: phi(:)
      real(wp), allocatable :: qsurf(:)
      !> Nuclear and AO-pair coupling to the Gaussian surface basis.
      real(wp), allocatable :: ucore(:, :)
      real(wp), allocatable :: bmat(:, :, :)
      logical :: ready = .false.
   end type cosmo_cache

contains

function create_cosmo_input(dielectric_const, model, rvdw, rscale, nang, &
      & monopoles, dipoles, quadrupoles, full_density) result(self)
   real(wp), intent(in) :: dielectric_const
   integer, intent(in), optional :: model
   real(wp), intent(in), optional :: rvdw(:)
   real(wp), intent(in), optional :: rscale
   integer, intent(in), optional :: nang
   logical, intent(in), optional :: monopoles, dipoles, quadrupoles
   logical, intent(in), optional :: full_density
   type(cosmo_input) :: self

   self%dielectric_const = dielectric_const
   if (present(model)) self%model = model
   if (present(rvdw)) self%rvdw = rvdw
   if (present(rscale)) self%rscale = rscale
   if (present(nang)) self%nang = nang
   if (present(monopoles)) self%monopoles = monopoles
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

   type(static_radius_type) :: gauss_radii

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
   self%monopoles = input%monopoles
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
   if (allocated(input%rvdw)) then
      if (maxval(mol%id) > size(input%rvdw)) then
         call fatal_error(error, "COSMO radii do not cover all species")
         return
      end if
      allocate(self%rvdw(mol%nat))
      self%rvdw = input%rscale * input%rvdw(mol%id)
   else
      ! Resolve moist's ORCA-like default radii once, just as in ddX.
      call new_gauss_radii(gauss_radii)
      call gauss_radii%update(mol, error)
      if (allocated(error)) return
      allocate(self%rvdw(mol%nat), source=input%rscale*gauss_radii%f0)
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
   real(wp) :: vec(3), r1, r3, k0, k1(3), k2(6)

   call taint(cache, ptr)
   ptr%ready = .false.

   ! sVDW drop cavity (previous implementation):
   ! call svdw%new(blend_k=6.5_wp)
   ! call new_cavity_drop(ptr%cavity, verbose=0, nleb=self%nang, &
   !    & radius_model=radii, lsf_model=svdw, error=error)

   ! iSwiG cavity. self%rvdw was resolved in new_cosmo from either explicitly
   ! supplied radii or moist's ORCA-like Gaussian radii, as done for ddX.
   call new_radii_custom_atoms(self%rvdw, radii, error)
   if (allocated(error)) return
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
   call ptr%cavity%get_amat(ptr%amat, error)
   if (allocated(error)) return

   if (self%monopoles) allocate(ptr%c0(ptr%cavity%ngrid, mol%nat))
   if (self%dipoles) allocate(ptr%c1(ptr%cavity%ngrid, 3*mol%nat))
   if (self%quadrupoles) allocate(ptr%c2(ptr%cavity%ngrid, 6*mol%nat))
   allocate(ptr%phi(ptr%cavity%ngrid), source=0.0_wp)
   allocate(ptr%qsurf(ptr%cavity%ngrid), source=0.0_wp)

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

   if (self%monopoles .or. self%dipoles .or. self%quadrupoles) then
      !$omp parallel do collapse(2) default(none) &
      !$omp shared(self, ptr, mol) private(iat, igrid, vec, k0, k1, k2)
      do iat = 1, mol%nat
         do igrid = 1, ptr%cavity%ngrid
            vec = ptr%cavity%xyz(:, igrid) - mol%xyz(:, iat)
            call get_gaussian_multipole_kernels(vec, ptr%cavity%xi(igrid), &
               & k0, k1, k2)
            if (self%monopoles) ptr%c0(igrid, iat) = k0
            if (self%dipoles) ptr%c1(igrid, 3*iat-2:3*iat) = k1
            if (self%quadrupoles) ptr%c2(igrid, 6*iat-5:6*iat) = k2
         end do
      end do
      !$omp end parallel do
   end if
   ptr%ready = .true.
end subroutine update

subroutine get_potential(self, mol, cache, wfn, pot)
   class(cosmo_solvation), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   type(potential_type), intent(inout) :: pot

   type(cosmo_cache), pointer :: ptr
   integer :: iat, ic
   real(wp), allocatable :: vat(:), vdp(:, :), vqp(:, :)

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
      if (self%monopoles) then
         vat = matmul(transpose(ptr%c0), ptr%qsurf)
         pot%vat(:, 1) = pot%vat(:, 1) + vat
      end if
      if (self%dipoles) then
         vdp = reshape(matmul(transpose(ptr%c1), ptr%qsurf), [3, mol%nat])
         pot%vdp(:, :, 1) = pot%vdp(:, :, 1) + vdp
      end if
      if (self%quadrupoles) then
         vqp = reshape(matmul(transpose(ptr%c2), ptr%qsurf), [6, mol%nat])
         pot%vqp(:, :, 1) = pot%vqp(:, :, 1) + vqp
      end if
   end if
end subroutine get_potential

subroutine get_energy(self, mol, cache, wfn, energies)
   class(cosmo_solvation), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(inout) :: energies(:)

   type(cosmo_cache), pointer :: ptr
   integer :: iat
   real(wp), allocatable :: vat(:), vdp(:, :), vqp(:, :)

   call view(cache, ptr)
   if (.not.associated(ptr)) return
   if (.not.ptr%ready) return
   call solve_surface(self, ptr, wfn)
   if (.not.ptr%ready) return
   if (self%full_density) then
      energies(:) = energies(:) + 0.5_wp*dot_product(ptr%qsurf, ptr%phi)/real(mol%nat, wp)
   else
      if (self%monopoles) vat = matmul(transpose(ptr%c0), ptr%qsurf)
      if (self%dipoles) vdp = reshape(matmul(transpose(ptr%c1), ptr%qsurf), [3, mol%nat])
      if (self%quadrupoles) vqp = reshape(matmul(transpose(ptr%c2), ptr%qsurf), [6, mol%nat])
      do iat = 1, mol%nat
         if (self%monopoles) energies(iat) = energies(iat) + 0.5_wp*wfn%qat(iat, 1)*vat(iat)
         if (self%dipoles) energies(iat) = energies(iat) + &
            & 0.5_wp*dot_product(wfn%dpat(:, iat, 1), vdp(:, iat))
         if (self%quadrupoles) energies(iat) = energies(iat) + &
            & 0.5_wp*dot_product(wfn%qpat(:, iat, 1), vqp(:, iat))
      end do
   end if
   
end subroutine get_energy

!> Coulomb coupling of a monopole, dipole, and traceless Cartesian
!> quadrupole to one normalized spherical Gaussian surface function.
pure subroutine get_gaussian_multipole_kernels(vec, xi, c0, c1, c2)
   real(wp), intent(in) :: vec(3)
   real(wp), intent(in) :: xi
   real(wp), intent(out) :: c0
   real(wp), intent(out) :: c1(3)
   real(wp), intent(out) :: c2(6)

   real(wp) :: r1, r2, gaussian, erf_term
   real(wp) :: dipole_kernel, quadrupole_kernel
   real(wp), parameter :: sqrtpi = sqrt(acos(-1.0_wp))

   r2 = sum(vec**2)
   r1 = sqrt(r2)
   gaussian = exp(-(xi*r1)**2)
   erf_term = erf(xi*r1)

   c0 = erf_term/r1
   dipole_kernel = erf_term/(r1*r2) &
      & - 2.0_wp*xi*gaussian/(sqrtpi*r2)
   quadrupole_kernel = erf_term/(r1*r2*r2) &
      & - 2.0_wp*xi*gaussian/(sqrtpi*r2*r2) &
      & - 4.0_wp*xi**3*gaussian/(3.0_wp*sqrtpi*r2)
   c1 = dipole_kernel*vec
   c2 = [vec(1)*vec(1), 2.0_wp*vec(1)*vec(2), &
      & vec(2)*vec(2), 2.0_wp*vec(1)*vec(3), &
      & 2.0_wp*vec(2)*vec(3), vec(3)*vec(3)]*quadrupole_kernel
end subroutine get_gaussian_multipole_kernels

!> Form the solute potential and solve for apparent surface charges.
subroutine solve_surface(self, ptr, wfn)
   class(cosmo_solvation), intent(in) :: self
   type(cosmo_cache), intent(inout) :: ptr
   type(wavefunction_type), intent(in) :: wfn

   type(error_type), allocatable :: error, write_error
   integer :: ic
   real(wp), allocatable :: density(:, :)
   real(wp), parameter :: solver_tol = 1.0e-10_wp
   integer, parameter :: solver_maxiter = 1000

   if (self%full_density) then
      ptr%phi = matmul(ptr%ucore, wfn%n0at)
      density = sum(wfn%density, dim=3)
      do ic = 1, size(ptr%phi)
         ptr%phi(ic) = ptr%phi(ic) - sum(ptr%bmat(ic, :, :)*density)
      end do
   else
      ptr%phi = 0.0_wp
      if (self%monopoles) ptr%phi = matmul(ptr%c0, wfn%qat(:, 1))
      if (self%dipoles) ptr%phi = ptr%phi + &
         & matmul(ptr%c1, reshape(wfn%dpat(:, :, 1), [size(ptr%c1, 2)]))
      if (self%quadrupoles) ptr%phi = ptr%phi + &
         & matmul(ptr%c2, reshape(wfn%qpat(:, :, 1), [size(ptr%c2, 2)]))
   end if
   call solve_pcm_iterative(ptr%amat, -self%feps*ptr%phi, ptr%qsurf, &
      & solver_tol, solver_maxiter, error)
   if (allocated(error)) ptr%ready = .false.
   call write_cpcm_file(ptr%cavity, ptr%phi, ptr%qsurf, self%feps, write_error)
end subroutine solve_surface

!> Write the current iSwiG surface in ORCA-compatible CPCM table format.
subroutine write_cpcm_file(cavity, phi, qsurf, feps, error)
   type(cavity_type_iswig), intent(in) :: cavity
   real(wp), intent(in) :: phi(:)
   real(wp), intent(in) :: qsurf(:)
   real(wp), intent(in) :: feps
   type(error_type), allocatable, intent(out) :: error

   integer :: igrid, io, stat
   character(len=512) :: iomsg
   real(wp) :: charge

   if (size(phi) /= cavity%ngrid .or. size(qsurf) /= cavity%ngrid) then
      call fatal_error(error, "CPCM surface data dimensions do not match the cavity")
      return
   end if
   if (abs(feps) <= epsilon(1.0_wp)) then
      call fatal_error(error, "CPCM dielectric scaling factor is zero")
      return
   end if

   open(newunit=io, file="tblite.cpcm", status="replace", action="write", &
      & iostat=stat, iomsg=iomsg)
   if (stat /= 0) then
      call fatal_error(error, "Could not open tblite.cpcm: "//trim(iomsg))
      return
   end if

   write(io, '(a)', iostat=stat, iomsg=iomsg) &
      & "          X                 Y                 Z               area"// &
      & "            potential          charge            w_leb"// &
      & "             Switch_F          G_width       atom"
   if (stat == 0) then
      do igrid = 1, cavity%ngrid
         ! Moist stores the area, Lebedev weight, switching function, and
         ! Gaussian width in the convention printed by ORCA. ORCA prints
         ! conductor charges before application of the dielectric scaling.
         charge = qsurf(igrid)/feps
         write(io, '(3f18.9,3f18.9,f18.9,2f20.14,i8)', &
            & iostat=stat, iomsg=iomsg) cavity%xyz(:, igrid), &
            & cavity%a(igrid), phi(igrid), charge, &
            & cavity%wleb(igrid), cavity%f(igrid), &
            & cavity%xi(igrid), cavity%owner(igrid) - 1
         if (stat /= 0) exit
      end do
   end if

   close(io)
   if (stat /= 0) then
      call fatal_error(error, "Could not write tblite.cpcm: "//trim(iomsg))
      return
   end if
end subroutine write_cpcm_file

pure function variable_info(self) result(info)
   class(cosmo_solvation), intent(in) :: self
   type(scf_info) :: info
   if (self%full_density) then
      info = scf_info(density=orbital_resolved)
   else
      info = scf_info(charge=merge(atom_resolved, not_used, self%monopoles), &
      & dipole=merge(atom_resolved, not_used, self%dipoles), &
      & quadrupole=merge(atom_resolved, not_used, self%quadrupoles))
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
