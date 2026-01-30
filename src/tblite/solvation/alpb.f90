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

!> @file tblite/solvation/alpb.f90
!> Provides the analytical linearized Poission-Boltzmann model.

!> Analytical linearized Poisson-Boltzmann implicit solvation model.
!>
!> Implements a reaction field model of the Generalized Born type.
module tblite_solvation_alpb
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3
   use tblite_blas, only : dot, gemv, symv
   use tblite_container_cache, only : container_cache
   use tblite_mesh_lebedev, only : grid_size, get_angular_grid, list_bisection
   use tblite_scf_info, only : scf_info, atom_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_solvation_born, only : born_integrator, new_born_integrator
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type
   use tblite_solvation_cm5, only : get_cm5_charges
   use tblite_solvation_kernel, only : kernel_type, new_kernel, kernel_enum, kernel_enum_type
   implicit none
   private

   public :: new_alpb
   public :: born_kernel

   ! Alias for backward compatibility
   type(kernel_enum_type), parameter :: born_kernel = kernel_enum


   !> Input for ALPB solvation
   type, public :: alpb_input
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> Scaling factor for Born radii
      real(wp) :: born_scale = 1.0_wp
      !> Offset parameter for Born radii integration
      real(wp) :: born_offset = 0.0_wp
      !> Van-der-Waals radii for all species
      real(wp), allocatable :: rvdw(:)
      !> Dielectric descreening parameter
      real(wp), allocatable :: descreening(:)
      !> Interaction kernel
      integer :: kernel = born_kernel%p16
      !> Use analytical linearized Poisson-Boltzmann model
      logical :: alpb = .false.
      !> Solvent for parameter selection
      character(len=:), allocatable :: solvent
      !> Whether or not to use multipoles
      logical :: do_multipoles = .false.
   end type alpb_input

   !> Provide constructor for ALPB input
   interface alpb_input
      module procedure :: create_alpb_input
   end interface alpb_input

   !> Definition of ALPB/GBSA model
   type, public, extends(solvation_type) :: alpb_solvation
      !> Dielectric function
      real(wp) :: keps
      !> Analytical linearized Poisson-Boltzmann constant
      real(wp) :: alpbet
      !> Integrator for Born radii
      type(born_integrator) :: gbobc
      !> Kernel instance
      class(kernel_type), allocatable :: kernel
      !> Use CM5 charges (GFN1-xTB compatibility)
      logical :: useCM5 = .false.
      !> Whether or not to use multipoles
      logical :: do_multipoles
   contains
      !> Update cache from container
      procedure :: update
      !> Return dependency on density
      procedure :: variable_info
      !> Get solvation energy
      procedure :: get_energy
      !> Get solvation potential
      procedure :: get_potential
      !> Get solvation gradient
      procedure :: get_gradient
   end type alpb_solvation

   !> Provide constructor for ALPB solvation
   interface alpb_solvation
      module procedure :: create_alpb
   end interface alpb_solvation


   !> Restart data for ALPB calculation
   type, public :: alpb_cache
      !> Screening matrix
      real(wp), allocatable :: jmat(:, :)
      !> Scratch array for screening potential intermediates
      real(wp), allocatable :: vat(:)
      !> Born radii
      real(wp), allocatable :: rad(:)
      !> Derivatives of Born radii w.r.t. cartesian displacements
      real(wp), allocatable :: draddr(:, :, :)
      !> Scratch workspace for gradient construction
      real(wp), allocatable :: scratch(:)
      !> Workspace for atomic charges
      real(wp), allocatable :: qscratch(:)
      !> CM5 charges (only required for GFN1 compatibility)
      real(wp), allocatable :: cm5(:)
      !> CM5 charge derivatives
      real(wp), allocatable :: dcm5dr(:,:,:)
      !> Multipole interaction matrix for charges and dipoles
      real(wp), allocatable :: amat_sd(:, :, :)
      !> Multipole interaction matrix for dipoles and dipoles
      real(wp), allocatable :: amat_dd(:, :, :, :)
      !> Multipole interaction matrix for charges and quadrupoles
      real(wp), allocatable :: amat_sq(:, :, :)
      !> Multipole interaction matrix for dipoles and quadrupoles
      real(wp), allocatable :: amat_dq(:, :, :, :)
      !> Multipole interaction matrix for quadrupoles and quadrupoles
      real(wp), allocatable :: amat_qq(:, :, :, :)
      
   end type alpb_cache


   !> Identifier for container
   character(len=*), parameter :: label = "alpb/gbsa reaction field model"
   character(len=*), parameter :: multipole_label = "alpb/gbsa reaction field model with atomic multipole interactions"
   real(wp), parameter :: alpha_alpb = 0.571412_wp


contains


!> Consturctor for ALPB input to properly assign allocatable strings
function create_alpb_input(dielectric_const, solvent, alpb, kernel, do_multipoles) result(self)
   !> Dielectric constant
   real(wp), intent(in) :: dielectric_const
   !> Solvent for parameter selection
   character(len=*), intent(in), optional :: solvent
   !> Use analytical linearized Poisson-Boltzmann model
   logical, intent(in), optional :: alpb
   !> Interaction kernel
   integer, intent(in), optional :: kernel
   !> Whether or not to use multipoles
   logical, intent(in), optional :: do_multipoles

   type(alpb_input) :: self

   self%dielectric_const = dielectric_const

   if (present(solvent)) then 
      self%solvent = solvent
   end if

   if (present(alpb)) then 
      self%alpb = alpb
   end if

   if (present(kernel)) then 
      self%kernel = kernel
   end if

   if (do_multipoles) then 
      self%do_multipoles = .true.
   else 
      self%do_multipoles = .false.
   end if

end function create_alpb_input


!> Create new ALPB solvation model
subroutine new_alpb(self, mol, input, method)
   !> Instance of the solvation model
   type(alpb_solvation), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ALPB solvation
   type(alpb_input), intent(in) :: input
   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   real(wp), allocatable :: rvdw(:)

   self%do_multipoles = input%do_multipoles
   if (self%do_multipoles) then
      self%label = multipole_label
   else
      self%label = label
   end if
   self%alpbet = merge(alpha_alpb / input%dielectric_const, 0.0_wp, input%alpb)
   self%keps = (1.0_wp/input%dielectric_const - 1.0_wp) / (1.0_wp + self%alpbet)
   self%kernel = new_kernel(input%kernel, self%keps)
   if (allocated(input%solvent) .and. present(method)) then
      self%useCM5 = trim(method) == 'gfn1'
   endif

   if (allocated(input%rvdw)) then
      rvdw = input%rvdw
   else
      rvdw = get_vdw_rad_cosmo(mol%num)
   end if

   call new_born_integrator(self%gbobc, mol, rvdw, descreening=input%descreening, &
      & born_scale=input%born_scale, born_offset=input%born_offset)
end subroutine new_alpb


!> Type constructor for ALPB solvation
function create_alpb(mol, input, method) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ALPB solvation
   type(alpb_input), intent(in) :: input
   !> Method for parameter selection
   character(len=*), intent(in), optional :: method
   !> Instance of the solvation model
   type(alpb_solvation) :: self

   call new_alpb(self, mol, input, method)
end function create_alpb


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(alpb_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(alpb_cache), pointer :: ptr
   real(wp) :: adet

   call taint(cache, ptr)

   if (.not.allocated(ptr%jmat)) then
      allocate(ptr%jmat(mol%nat, mol%nat), source=0.0_wp)
   end if
   if (.not.allocated(ptr%vat)) then
      allocate(ptr%vat(mol%nat), source=0.0_wp)
   end if
   if (.not.allocated(ptr%rad)) then
      allocate(ptr%rad(mol%nat), source=0.0_wp)
   end if
   if (.not.allocated(ptr%draddr)) then
      allocate(ptr%draddr(3, mol%nat, mol%nat), source=0.0_wp)
   end if
   if (.not.allocated(ptr%qscratch))then
      allocate(ptr%qscratch(mol%nat), source=0.0_wp)
   endif 
   if (self%useCM5)then
      if (.not.allocated(ptr%cm5))then
         allocate(ptr%cm5(mol%nat), source=0.0_wp) 
      endif
      if (.not.allocated(ptr%dcm5dr))then
         allocate(ptr%dcm5dr(3, mol%nat, mol%nat), source=0.0_wp)
      endif
      call get_cm5_charges(mol, ptr%cm5, ptr%dcm5dr)
   endif
   if (self%useCM5.and..not.allocated(ptr%scratch))then
      allocate(ptr%scratch(mol%nat))
   endif
   call self%gbobc%get_rad(mol, ptr%rad, ptr%draddr)
   call self%kernel%add_kernel_mat(mol%nat, mol%xyz, ptr%rad, ptr%jmat)

   if (self%alpbet > 0.0_wp) then
      call get_adet(mol%nat, mol%xyz, self%gbobc%vdwr, adet)

      ptr%jmat(:mol%nat, :mol%nat) = ptr%jmat(:mol%nat, :mol%nat) &
         & + self%keps * self%alpbet / adet
   end if

   if (self%do_multipoles) then
      if (.not.allocated(ptr%amat_sd)) then
         allocate(ptr%amat_sd(3, mol%nat, mol%nat), source=0.0_wp)
      end if
      if (.not.allocated(ptr%amat_dd)) then
         allocate(ptr%amat_dd(3, mol%nat, 3, mol%nat), source=0.0_wp)
      end if
      if (.not.allocated(ptr%amat_sq)) then
         allocate(ptr%amat_sq(6, mol%nat, mol%nat), source=0.0_wp)
      end if
      if (.not.allocated(ptr%amat_dq)) then
         allocate(ptr%amat_dq(3, mol%nat, 6, mol%nat), source=0.0_wp)
      end if
      if (.not.allocated(ptr%amat_qq)) then
         allocate(ptr%amat_qq(6, mol%nat, 6, mol%nat), source=0.0_wp)
      end if

      call get_multipole_matrices(self, mol, mol%xyz, self%keps, ptr%rad, &
         & ptr%amat_sd, ptr%amat_dd, ptr%amat_sq, ptr%amat_dq, ptr%amat_qq)
   end if

end subroutine update


!> Get solvation energy
subroutine get_energy(self, mol, cache, wfn, energies)
   class(alpb_solvation), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(inout) :: energies(:)

   real(wp), allocatable :: vd(:, :), vq(:, :)
   type(alpb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%useCM5) then
      ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
   else
      ptr%qscratch(:) = wfn%qat(:, 1)
   end if

   ! charge-charge
   call symv(ptr%jmat, ptr%qscratch(:), ptr%vat, alpha=0.5_wp)
   energies(:) = energies + ptr%vat * ptr%qscratch(:) 

   if (self%do_multipoles) then
      allocate(vd(3, mol%nat), vq(6, mol%nat), source=0.0_wp)

      ! charge-dipole 
      call gemv(ptr%amat_sd, wfn%qat(:, 1), vd)                                  ! SD * q
      ! dipole-dipole
      call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), vd, beta=1.0_wp, alpha=0.5_wp)   ! + 1/2 DD * mu
      ! dipole-quadrupole
      call gemv(ptr%amat_dq, wfn%qpat(:, :, 1), vd, beta=1.0_wp, alpha=1.0_wp)   ! + DQ * Q   (NO 1/2)
      ! charge-quadrupole
      call gemv(ptr%amat_sq, wfn%qat(:, 1), vq)                                  ! SQ * q
      ! quadrupole-quadrupole
      call gemv(ptr%amat_qq, wfn%qpat(:, :, 1), vq, beta=1.0_wp, alpha=0.5_wp)   ! + 1/2 QQ * Q
   
      energies(:) = energies + sum(wfn%dpat(:, :, 1) * vd, 1) + sum(wfn%qpat(:, :, 1) * vq, 1)
   end if

end subroutine get_energy



!> Get solvation potential
subroutine get_potential(self, mol, cache, wfn, pot)
   class(alpb_solvation), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   type(potential_type), intent(inout) :: pot

   type(alpb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%useCM5) then
      ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
   else
      ptr%qscratch(:) = wfn%qat(:, 1)
   end if

   call symv(ptr%jmat, ptr%qscratch(:), pot%vat(:, 1), beta=1.0_wp)

   if (self%do_multipoles) then
      call gemv(ptr%amat_sd, wfn%qat(:, 1), pot%vdp(:, :, 1), beta=1.0_wp)
      call gemv(ptr%amat_sd, wfn%dpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

      call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), pot%vdp(:, :, 1), beta=1.0_wp)

      call gemv(ptr%amat_sq, wfn%qat(:, 1), pot%vqp(:, :, 1), beta=1.0_wp)
      call gemv(ptr%amat_sq, wfn%qpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

      call gemv(ptr%amat_dq, wfn%qpat(:, :, 1), pot%vdp(:, :, 1), beta=1.0_wp)           
      call gemv(ptr%amat_dq, wfn%dpat(:, :, 1), pot%vqp(:, :, 1), beta=1.0_wp, trans="T")

      call gemv(ptr%amat_qq, wfn%qpat(:, :, 1), pot%vqp(:, :, 1), beta=1.0_wp)
   end if

end subroutine get_potential


!> Get solvation gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the solvation model
   class(alpb_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the solvation free energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the solvation free energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   type(alpb_cache), pointer :: ptr
   real(wp) :: energy

   energy = 0.0_wp
   call view(cache, ptr)

   if(self%useCM5)then
      ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
   else
      ptr%qscratch(:) = wfn%qat(:, 1)
   endif

   call self%kernel%add_kernel_deriv(mol%nat, mol%xyz, ptr%qscratch(:), &
      & ptr%rad, ptr%draddr, energy, gradient)


   if (self%do_multipoles) then
      call self%kernel%add_kernel_deriv_multipole_contributions(mol%nat, mol%xyz, ptr%qscratch(:), &
         wfn%dpat(:,:,1), wfn%qpat(:,:,1), ptr%rad, ptr%draddr, gradient)
   end if

   if (self%alpbet > 0.0_wp) then
      call get_adet_deriv(mol%nat, mol%xyz, self%gbobc%vdwr, self%keps*self%alpbet, &
         & ptr%qscratch(:), gradient)
   end if

   if(self%useCM5)then
      call gemv(ptr%jmat, ptr%qscratch, ptr%scratch)
      call gemv(ptr%dcm5dr, ptr%scratch, gradient, beta=1.0_wp)
   endif

end subroutine get_gradient


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the solvation model
   class(alpb_solvation), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(alpb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(alpb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(alpb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(alpb_cache)
      ptr => target
   end select
end subroutine view


subroutine get_adet(nat, xyz, rad, aDet)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic radii
   real(wp), intent(in) :: rad(:)
   !> Shape descriptor of the structure
   real(wp), intent(out) :: aDet

   integer :: iat
   real(wp) :: r2, rad2, rad3, vol, vec(3), center(3), inertia(3, 3)
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   vol = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nat
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vol = vol + rad3
      center(:) = center + xyz(:, iat) * rad3
   end do
   center = center / vol

   inertia(:, :) = 0.0_wp
   do iat = 1, nat
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do

   aDet = sqrt(matdet_3x3(inertia)**(1.0_wp/3.0_wp)/(tof*vol))

end subroutine get_adet


subroutine get_adet_deriv(nAtom, xyz, rad, kEps, qvec, gradient)
   !> Number of atoms
   integer, intent(in) :: nAtom
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic radii
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kEps
   real(wp), intent(in) :: qvec(:)
   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iat
   real(wp) :: r2, rad2, rad3, vol, vec(3), center(3), inertia(3, 3), aDet
   real(wp) :: aDeriv(3, 3), qtotal
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   qtotal = 0.0_wp
   vol = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vol = vol + rad3
      center(:) = center + xyz(:, iat) * rad3
      qtotal = qtotal + qvec(iat)
   end do
   center = center / vol

   inertia(:, :) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do
   aDet = sqrt(matdet_3x3(inertia)**(1.0_wp/3.0_wp)/(tof*vol))

   aDeriv(:, :) = reshape([&
      & inertia(1,1)*(inertia(2,2)+inertia(3,3))-inertia(1,2)**2-inertia(1,3)**2, &
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
      & inertia(2,2)*(inertia(1,1)+inertia(3,3))-inertia(1,2)**2-inertia(2,3)**2, &
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
      & inertia(3,3)*(inertia(1,1)+inertia(2,2))-inertia(1,3)**2-inertia(2,3)**2],&
      & shape=[3, 3]) * (250.0_wp / (48.0_wp * vol**3 * aDet**5)) &
      & * (-0.5_wp * kEps * qtotal**2 / aDet**2)

   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      gradient(:, iat) = gradient(:, iat) + rad3 * matmul(aderiv, vec)
   end do

end subroutine get_adet_deriv

!> Build multipole interaction matrices
!>
!> Computes the interaction matrices for atomic multipoles (dipoles and quadrupoles)
!> based on derivatives of the solvation kernel K_ij between atom pairs.
!> These matrices encode how multipoles on atom i induce potentials/fields at atom j.
subroutine get_multipole_matrices(self, mol, xyz, keps, brad, &
      amat_sd, amat_dd, amat_sq, amat_dq, amat_qq)

   !> Instance of ALPB solvation model
   class(alpb_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in)  :: mol
   !> Cartesian coordinates
   real(wp), intent(in)              :: xyz(:, :)
   !> Dielectric screening factor
   real(wp), intent(in)              :: keps
   !> Born radii 
   real(wp), intent(in)              :: brad(:)
   !> Charge-dipole interaction matrix (3, nat, nat)
   real(wp), contiguous, intent(inout) :: amat_sd(:, :, :)
   !> Dipole-dipole interaction matrix (3, nat, 3, nat)
   real(wp), contiguous, intent(inout) :: amat_dd(:, :, :, :)
   !> Charge-quadrupole interaction matrix (6, nat, nat)
   real(wp), contiguous, intent(inout) :: amat_sq(:, :, :)
   !> Dipole-quadrupole interaction matrix (3, nat, 6, nat)
   real(wp), contiguous, intent(inout) :: amat_dq(:, :, :, :)
   !> Quadrupole-quadrupole interaction matrix (6, nat, 6, nat)
   real(wp), contiguous, intent(inout) :: amat_qq(:, :, :, :)

   !> Number of atoms
   integer :: nat
   !> Loop indices for atom pairs
   integer :: i, j
   !> Cartesian indices for tensor components
   integer :: alpha, beta, gamma, epsilon
   !> Packed indices for quadrupole components
   integer :: p, q
   !> Packing factors for symmetric tensor storage
   integer :: fab, fge

   !> Distance vector between atoms i and j
   real(wp) :: r(3)
   !> Interatomic distance
   real(wp) :: rij
   !> Inverse distance (1/r)
   real(wp) :: invr
   !> Inverse squared distance (1/r²)
   real(wp) :: invr2
   !> Unit vector along r_ij
   real(wp) :: uvec(3)

   !> First derivative of kernel: ∂K_ij/∂r_j
   real(wp) :: dKij_drj(3)
   !> Second derivative of kernel: ∂²K_ij/∂r_j²
   real(wp) :: d2Kij_drj2(3,3)
   !> Third derivative of kernel: ∂³K_ij/∂r_j³
   real(wp) :: d3Kij_drj3(3,3,3)
   !> Fourth derivative of kernel: ∂⁴K_ij/∂r_j⁴
   real(wp) :: d4Kij_drj4(3,3,3,3)

   !> Hessian times unit vector (for radial derivatives)
   real(wp) :: hu(3)
   !> Radial gradient component: ∇K · u
   real(wp) :: gpar
   !> Radial second derivative: u^T ∇∇K u
   real(wp) :: kpp
   !> Coefficient for quadrupole construction: (kpp + gpar/r) / 3
   real(wp) :: coef
   !> Isotropic correction term for quadrupole-quadrupole: coef/r²
   real(wp) :: s5

   !> Quadrupole tensor diagonal components
   real(wp) :: q11, q22, q33
   !> Quadrupole tensor off-diagonal components
   real(wp) :: q12, q13, q23
   !> Packed quadrupole tensor (6 components: xx, xy, yy, xz, yz, zz)
   real(wp) :: tc(6)

   !> 3×3 identity matrix
   real(wp) :: i3(3,3)
   !> Symmetry factor for fourth derivative terms
   real(wp) :: sym2
   !> Quadrupole-quadrupole interaction element
   real(wp) :: wabge

   !> Coordinates of atom i (source)
   real(wp) :: ri(3)
   !> Coordinates of atom j (response)
   real(wp) :: rj(3)
   !> Born radius of atom i
   real(wp) :: borni
   !> Born radius of atom j
   real(wp) :: bornj

   !> Packing arrays: map symmetric 3×3 matrix indices to 6-vector
   !> pa(p), pb(p) give (row, col) for packed index p
   !> Order: (1,1), (1,2), (2,2), (1,3), (2,3), (3,3) = (xx, xy, yy, xz, yz, zz)
   integer, parameter :: pa(6) = [1, 1, 2, 1, 2, 3]
   integer, parameter :: pb(6) = [1, 2, 2, 3, 3, 3]
   !> Packing factors: 1 for diagonal elements, 2 for off-diagonal (accounts for symmetry)
   integer, parameter :: pf(6) = [1, 2, 1, 2, 2, 1]

   nat = mol%nat

   ! Zero out the arrays before accumulating
   amat_sd = 0.0_wp
   amat_dd = 0.0_wp
   amat_sq = 0.0_wp
   amat_dq = 0.0_wp
   amat_qq = 0.0_wp

   ! Identity matrix
   i3 = 0.0_wp
   i3(1,1) = 1.0_wp
   i3(2,2) = 1.0_wp
   i3(3,3) = 1.0_wp

   do i = 1, nat
      do j = 1, nat
         if (i == j) cycle

         r(:) = xyz(:, i) - xyz(:, j)
         rij  = sqrt(dot_product(r, r))

         invr  = 1.0_wp / rij
         invr2 = invr * invr
         uvec(:) = r(:) * invr ! Unit vector along r_ij

         ! ------------------------------------------------------------
         ! all derivatives in this routine are defined w.r.t. atom j
         ! ------------------------------------------------------------
         rj(:) = xyz(:, j)
         ri(:) = xyz(:, i)
         bornj   = brad(j)
         borni   = brad(i)

         ! =================================================================================
         ! 1) Get charge-dipole interaction from first derivative of the interaction kernel
         ! =================================================================================
         call self%kernel%kernel_d1_pair(rj, ri, bornj, borni, dKij_drj)

         ! Index pattern (..., j, ..., i) here and in the following:
         ! Response on j due to source on i
         amat_sd(:, j, i) = amat_sd(:, j, i) + dKij_drj(:)

         ! ==================================================================================
         ! 2a) Get dipole-dipole interaction from second derivative of the interaction kernel
         !    \sum_{alpha,beta} - d2Kij / (d_rj,alpha d_rj,beta)
         !   = \sum_{alpha,beta} - dKij/(d_rj,alpha) dKij/(d_rj,beta)
         !                                |                       |
         !                         (x,y,z) dipole          (x,y,z) dipole
         ! ==================================================================================
         call self%kernel%kernel_d2_pair(rj, ri, bornj, borni, d2Kij_drj2)

         do alpha = 1, 3
            do beta = 1, 3
               ! Mapping dipole at i to a dipole-type interaction at j
               amat_dd(alpha, j, beta, i) = amat_dd(alpha, j, beta, i) - d2Kij_drj2(alpha,beta)
            end do
         end do

         ! ==================================================================================
         ! 2b) Get charge-quadrupole interaction from second derivative of the interaction kernel
         !    \sum_{alpha,beta} - d2Kij / (d_rj,alpha d_rj,beta)
         !                                       |
         !                        (xx,xy,xz,yy,yz,zz) quadrupole
         ! ==================================================================================

         ! gpar = grad · u,  kpp = u^t hess u
         gpar  = dot_product(dKij_drj, uvec)   ! Radial component of the gradient
         hu(:) = matmul(d2Kij_drj2, uvec)        
         kpp   = dot_product(uvec, hu)  ! Second derivative along the radial direction
         coef = (kpp + gpar * invr) / 3.0_wp
         ! -> For Coulomb this produces coef = 1/r^3 

         ! q = coef * u u^t, pack with doubled off-diagonals
         q11 = coef * uvec(1) * uvec(1)
         q22 = coef * uvec(2) * uvec(2)
         q33 = coef * uvec(3) * uvec(3)
         q12 = coef * uvec(1) * uvec(2)
         q13 = coef * uvec(1) * uvec(3)
         q23 = coef * uvec(2) * uvec(3)
         ! -> For Coulomb this produces vec_beta vec_gamma / r^5, as in the multipole.f90 implementation

         ! Convention: double the off-diagonal elements
         tc(1) = q11
         tc(2) = 2.0_wp * q12
         tc(3) = q22
         tc(4) = 2.0_wp * q13
         tc(5) = 2.0_wp * q23
         tc(6) = q33

         ! Mapping charge at i to a quadrupole-type interaction at j
         amat_sq(:, j, i) = amat_sq(:, j, i) + tc(:)

         ! ==================================================================================
         ! 3) Get dipole-quadrupole interaction from third derivative of the interaction kernel
         !    \sum_{alpha,beta,gamma} - d3Kij / (d_rj,alpha d_rj,beta drj,gamma)
         !   = \sum_{alpha,beta,gamma} - dKij/(d_rj,alpha) d2Kij/(d_rj,beta drj,gamma)
         !                                    |                           |
         !                             (x,y,z) dipole        (xx,xy,xz,yy,yz,zz) quadrupole
         ! ==================================================================================
         call self%kernel%kernel_d3_pair(rj, ri, bornj, borni, d3Kij_drj3)

         ! quadrupole cotribution in lower triangular order: 11, 12, 22, 13, 23, 33
         ! Off-diagonal terms are doubled, factor 1/3 from isotropic correction (trace removal)
         do alpha = 1, 3
            amat_dq(alpha, j, 1, i) = amat_dq(alpha, j, 1, i) + (-(1.0_wp/3.0_wp)) * d3Kij_drj3(alpha,1,1)
            amat_dq(alpha, j, 2, i) = amat_dq(alpha, j, 2, i) + (-(1.0_wp/3.0_wp)) * 2.0_wp * d3Kij_drj3(alpha,1,2)
            amat_dq(alpha, j, 3, i) = amat_dq(alpha, j, 3, i) + (-(1.0_wp/3.0_wp)) * d3Kij_drj3(alpha,2,2)
            amat_dq(alpha, j, 4, i) = amat_dq(alpha, j, 4, i) + (-(1.0_wp/3.0_wp)) * 2.0_wp * d3Kij_drj3(alpha,1,3)
            amat_dq(alpha, j, 5, i) = amat_dq(alpha, j, 5, i) + (-(1.0_wp/3.0_wp)) * 2.0_wp * d3Kij_drj3(alpha,2,3)
            amat_dq(alpha, j, 6, i) = amat_dq(alpha, j, 6, i) + (-(1.0_wp/3.0_wp)) * d3Kij_drj3(alpha,3,3)
         end do

         ! ==================================================================================
         ! 4) Get quadrupole-quadrupole interaction from fourth derivative of the interaction kernel
         !    \sum_{alpha,beta,gamma, epsilon} - d4Kij / (d_rj,alpha d_rj,beta drj,gamma, drj,epsilon)
         !   = \sum_{alpha,beta,gamma,epsilon} - d2Kij/(d_rj,alpha d_rj,beta) d2Kij/(drj,gamma, drj,epsilon)
         !                                                    |                           |
         !                                     (xx,xy,xz,yy,yz,zz) quadrupole   (xx,xy,xz,yy,yz,zz) quadrupole
         ! ==================================================================================
         call self%kernel%kernel_d4_pair(rj, ri, bornj, borni, d4Kij_drj4)

         ! Isotropic correction coefficient
         s5 = coef * invr2

         ! pa(p) and pb(p) give the Cartesian for (alpha, beta) to match lower-triangular packing  
         ! pf(p) gives the factor to account for doubled off-diagonals in the packed representation
         do p = 1, 6
            alpha   = pa(p)
            beta   = pb(p)
            fab = pf(p)

            do q = 1, 6
               gamma   = pa(q)
               epsilon = pb(q)
               fge = pf(q)

               ! δab​δcd ​+ δac​δbd ​+ δad​δbc (via the identity)
               sym2 = i3(alpha,beta)*i3(gamma,epsilon) + i3(alpha,gamma)*i3(beta,epsilon) + i3(alpha,epsilon)*i3(beta,gamma)

               ! wabge ​= 1/3 ∂_abge​ Kij + 1/2​ (δab​δcd ​+ δac​δbd ​+ δad​δbc)*s5
               wabge = (1.0_wp/3.0_wp) * d4Kij_drj4(alpha,beta,gamma,epsilon) + 0.5_wp * sym2 * s5

               ! Mapping quadrupole at i to a quadrupole-type interaction at j
               ! Multiply by (fab * fcd) to account for the double off-diagonal packing of the quadrupole components
               amat_qq(p, j, q, i) = amat_qq(p, j, q, i) + real(fab*fge, wp) * wabge
            end do
         end do

      end do
   end do

end subroutine get_multipole_matrices






end module tblite_solvation_alpb
