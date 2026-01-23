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
   public :: get_multipole_matrix

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

   real(wp), parameter :: alpha_alpb = 0.571412_wp


contains


!> Consturctor for ALPB input to properly assign allocatable strings
function create_alpb_input(dielectric_const, solvent, alpb, kernel) result(self)
   !> Dielectric constant
   real(wp), intent(in) :: dielectric_const
   !> Solvent for parameter selection
   character(len=*), intent(in), optional :: solvent
   !> Use analytical linearized Poisson-Boltzmann model
   logical, intent(in), optional :: alpb
   !> Interaction kernel
   integer, intent(in), optional :: kernel

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

   self%label = label
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
      allocate(ptr%jmat(mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%vat)) then
      allocate(ptr%vat(mol%nat))
   end if
   if (.not.allocated(ptr%rad)) then
      allocate(ptr%rad(mol%nat))
   end if
   if (.not.allocated(ptr%draddr)) then
      allocate(ptr%draddr(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%qscratch))then
      allocate(ptr%qscratch(mol%nat))
   endif 
   if (self%useCM5)then
      if (.not.allocated(ptr%cm5))then
         allocate(ptr%cm5(mol%nat)) 
      endif
      if (.not.allocated(ptr%dcm5dr))then
         allocate(ptr%dcm5dr(3, mol%nat, mol%nat))
      endif
      call get_cm5_charges(mol, ptr%cm5, ptr%dcm5dr)
   endif
   if (self%useCM5.and..not.allocated(ptr%scratch))then
      allocate(ptr%scratch(mol%nat))
   endif
   if (.not.allocated(ptr%amat_sd)) then
      allocate(ptr%amat_sd(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%amat_dd)) then
      allocate(ptr%amat_dd(3, mol%nat, 3, mol%nat))
   end if
   if (.not.allocated(ptr%amat_sq)) then
      allocate(ptr%amat_sq(6, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%amat_dq)) then
      allocate(ptr%amat_dq(3, mol%nat, 6, mol%nat))
   end if
   if (.not.allocated(ptr%amat_qq)) then
      allocate(ptr%amat_qq(6, mol%nat, 6, mol%nat))
   end if
   
   call self%gbobc%get_rad(mol, ptr%rad, ptr%draddr)
   ptr%jmat(:, :) = 0.0_wp
   call self%kernel%add_kernel_mat(mol%nat, mol%xyz, ptr%rad, ptr%jmat)

   if (self%alpbet > 0.0_wp) then
      call get_adet(mol%nat, mol%xyz, self%gbobc%vdwr, adet)

      ptr%jmat(:mol%nat, :mol%nat) = ptr%jmat(:mol%nat, :mol%nat) &
         & + self%keps * self%alpbet / adet
   end if

   ! Compute multipole interaction matrix
   call get_multipole_matrix(self, mol, mol%xyz, self%keps, ptr%rad, ptr%draddr, &
      & ptr%amat_sd, ptr%amat_dd, ptr%amat_sq, ptr%amat_dq, ptr%amat_qq)
end subroutine update


!> Get solvation energy
subroutine get_energy(self, mol, cache, wfn, energies)
   class(alpb_solvation), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(container_cache), intent(inout) :: cache
   type(wavefunction_type), intent(in) :: wfn
   real(wp), intent(inout) :: energies(:)

   real(wp), allocatable :: vs(:), vd(:, :), vq(:, :)
   type(alpb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%useCM5) then
      ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
   else
      ptr%qscratch(:) = wfn%qat(:, 1)
   end if

   allocate(vs(mol%nat), vd(3, mol%nat), vq(6, mol%nat))
   vd(:,:) = 0.0_wp
   vq(:,:) = 0.0_wp

   ! charge-charge
   call symv(ptr%jmat, ptr%qscratch(:), ptr%vat, alpha=0.5_wp)

   call gemv(ptr%amat_sd, wfn%qat(:, 1), vd)                                  ! SD * q
   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), vd, beta=1.0_wp, alpha=0.5_wp)   ! + 1/2 DD * mu
   call gemv(ptr%amat_dq, wfn%qpat(:, :, 1), vd, beta=1.0_wp, alpha=1.0_wp)   ! + DQ * Q   (NO 1/2)

   call gemv(ptr%amat_sq, wfn%qat(:, 1), vq)                                  ! SQ * q
   call gemv(ptr%amat_qq, wfn%qpat(:, :, 1), vq, beta=1.0_wp, alpha=0.5_wp)   ! + 1/2 QQ * Q

   energies(:) = energies &
      + ptr%vat * ptr%qscratch(:) &
      + sum(wfn%dpat(:, :, 1) * vd, 1) &
      + sum(wfn%qpat(:, :, 1) * vq, 1)

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

   call gemv(ptr%amat_sd, wfn%qat(:, 1), pot%vdp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_sd, wfn%dpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

   call gemv(ptr%amat_dd, wfn%dpat(:, :, 1), pot%vdp(:, :, 1), beta=1.0_wp)

   call gemv(ptr%amat_sq, wfn%qat(:, 1), pot%vqp(:, :, 1), beta=1.0_wp)
   call gemv(ptr%amat_sq, wfn%qpat(:, :, 1), pot%vat(:, 1), beta=1.0_wp, trans="T")

   call gemv(ptr%amat_dq, wfn%qpat(:, :, 1), pot%vdp(:, :, 1), beta=1.0_wp)           
   call gemv(ptr%amat_dq, wfn%dpat(:, :, 1), pot%vqp(:, :, 1), beta=1.0_wp, trans="T")

   call gemv(ptr%amat_qq, wfn%qpat(:, :, 1), pot%vqp(:, :, 1), beta=1.0_wp)

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

   if (self%alpbet > 0.0_wp) then
      call get_adet_deriv(mol%nat, mol%xyz, self%gbobc%vdwr, self%kEps*self%alpbet, &
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


! subroutine get_multipole_matrix(self, mol, xyz, keps, brad, brdr, amat_sd, amat_dd, amat_sq)
!    class(alpb_solvation), intent(in) :: self
!    type(structure_type), intent(in)  :: mol
!    real(wp), intent(in)              :: xyz(:, :)
!    real(wp), intent(in)              :: keps
!    real(wp), intent(in)              :: brad(:)
!    real(wp), contiguous, intent(in)  :: brdr(:, :, :)
!    real(wp), contiguous, intent(inout) :: amat_sd(:, :, :)      ! (3,nat,nat)
!    real(wp), intent(inout)             :: amat_dd(:, :, :, :)    ! (3,nat,3,nat)
!    real(wp), intent(inout)             :: amat_sq(:, :, :)       ! (6,nat,nat)

!    integer :: i, j, nat
!    real(wp), allocatable :: dKdr_ij(:, :)          ! (3,nat)
!    real(wp), allocatable :: d2Kdr2_ij(:, :, :, :)  ! (3,nat,3,nat)

!    ! kept for interface compatibility
!    real(wp), allocatable :: brdr2(:, :, :, :, :)   ! (3,nat,3,nat,nat)
!    real(wp), allocatable :: temp(:), temp2(:,:,:)

!    real(wp), parameter :: tiny_r = 1.0e-14_wp
!    real(wp) :: R(3), rij, invr
!    real(wp) :: u(3), g(3), Hu(3), gpar, kpp, coef
!    real(wp) :: Q11, Q22, Q33, Q12, Q13, Q23
!    real(wp) :: tc(6)

!    nat = mol%nat

!    allocate(dKdr_ij(3, nat), source=0.0_wp)
!    allocate(d2Kdr2_ij(3, nat, 3, nat), source=0.0_wp)

!    allocate(brdr2(3, nat, 3, nat, nat), source=0.0_wp)
!    allocate(temp(nat), source=0.0_wp)
!    allocate(temp2(3, nat, nat), source=0.0_wp)

!    call self%gbobc%get_rad(mol, temp, temp2, brdr2)

!    do i = 1, nat
!       do j = 1, nat
!          if (i == j) cycle

!          R(:)  = xyz(:, i) - xyz(:, j)
!          rij   = sqrt(dot_product(R, R))
!          if (rij <= tiny_r) cycle
!          invr  = 1.0_wp / rij
!          u(:)  = R(:) * invr

!          ! ---- First derivatives (for SD and for SQ coefficient) ----
!          call self%kernel%compute_kernel_dkdr_ij(nat, xyz, brad, brdr, i, j, dKdr_ij)

!          g(:) = dKdr_ij(:, j)

!          ! SD: A_sd(:,j,i) += dK/dr_j
!          amat_sd(:, j, i) = amat_sd(:, j, i) + g(:)

!          ! ---- Second derivatives (for DD and for SQ coefficient) ----
!          call self%kernel%compute_kernel_d2kdr2_ij(nat, xyz, brad, brdr, brdr2, i, j, d2Kdr2_ij)

!          ! DD: A_dd(:,j,:,i) -= d2K/(dr_j dr_j)
!          amat_dd(:, j, :, i) = amat_dd(:, j, :, i) - d2Kdr2_ij(:, j, :, j)

!          ! ---- SQ (universal): build Q = c * u u^T, pack into tc ----
!          ! radial projections:
!          gpar = dot_product(g, u)

!          Hu(:) = matmul(d2Kdr2_ij(:, j, :, j), u)
!          kpp   = dot_product(u, Hu)

!          ! universal coefficient:
!          coef  = (kpp + gpar * invr) / 3.0_wp

!          ! symmetric tensor Q = coef * (u u^T)
!          Q11 = coef * u(1) * u(1)
!          Q22 = coef * u(2) * u(2)
!          Q33 = coef * u(3) * u(3)
!          Q12 = coef * u(1) * u(2)
!          Q13 = coef * u(1) * u(3)
!          Q23 = coef * u(2) * u(3)

!          ! pack in legacy order (xx, 2xy, yy, 2xz, 2yz, zz)
!          tc(1) = Q11
!          tc(2) = 2.0_wp * Q12
!          tc(3) = Q22
!          tc(4) = 2.0_wp * Q13
!          tc(5) = 2.0_wp * Q23
!          tc(6) = Q33

!          amat_sq(:, j, i) = amat_sq(:, j, i) + tc
!       end do
!    end do

!    deallocate(dKdr_ij, d2Kdr2_ij, brdr2, temp, temp2)
! end subroutine get_multipole_matrix

! subroutine get_multipole_matrix(self, mol, xyz, keps, brad, brdr, amat_sd, amat_dd, amat_sq, amat_dq)
!    class(alpb_solvation), intent(in) :: self
!    type(structure_type), intent(in)  :: mol
!    real(wp), intent(in)              :: xyz(:, :)
!    real(wp), intent(in)              :: keps
!    real(wp), intent(in)              :: brad(:)
!    real(wp), contiguous, intent(in)  :: brdr(:, :, :)

!    real(wp), contiguous, intent(inout) :: amat_sd(:, :, :)      ! (3,nat,nat)
!    real(wp), intent(inout)             :: amat_dd(:, :, :, :)    ! (3,nat,3,nat)
!    real(wp), intent(inout)             :: amat_sq(:, :, :)       ! (6,nat,nat)
!    real(wp), intent(inout)             :: amat_dq(:, :, :, :)    ! (3,nat,6,nat)

!    integer :: i, j, nat, a

!    real(wp), allocatable :: dKdr_ij(:, :)                 ! (3,nat)
!    real(wp), allocatable :: d2Kdr2_ij(:, :, :, :)         ! (3,nat,3,nat)
!    real(wp), allocatable :: d3Kdr3_ij(:, :, :, :, :, :)   ! (3,nat,3,nat,3,nat)

!    ! kept for interface compatibility
!    real(wp), allocatable :: brdr2(:, :, :, :, :)          ! (3,nat,3,nat,nat)
!    real(wp), allocatable :: brdr3(:, :, :, :, :, :, :)       ! (3,nat,3,nat,nat,nat)  <-- one more dimension than brdr2
!    real(wp), allocatable :: temp(:), temp2(:,:,:)

!    real(wp), parameter :: tiny_r = 1.0e-14_wp
!    real(wp) :: R(3), rij, invr
!    real(wp) :: uvec(3), g(3), Hu(3), gpar, kpp, coef
!    real(wp) :: Q11, Q22, Q33, Q12, Q13, Q23
!    real(wp) :: tc(6)

!    real(wp) :: U_dq(3,3,3)   ! dipole-quadrupole tensor U_{a,bc}

!    nat = mol%nat

!    allocate(dKdr_ij(3, nat), source=0.0_wp)
!    allocate(d2Kdr2_ij(3, nat, 3, nat), source=0.0_wp)
!    allocate(d3Kdr3_ij(3, nat, 3, nat, 3, nat), source=0.0_wp)

!    allocate(brdr2(3, nat, 3, nat, nat), source=0.0_wp)
!    allocate(brdr3(3, nat, 3, nat, 3, nat, nat), source=0.0_wp)   ! frozen radii => stays zero
!    allocate(temp(nat), source=0.0_wp)
!    allocate(temp2(3, nat, nat), source=0.0_wp)

!    call self%gbobc%get_rad(mol, temp, temp2, brdr2)

!    do i = 1, nat
!       do j = 1, nat
!          if (i == j) cycle

!          R(:)  = xyz(:, i) - xyz(:, j)
!          rij   = sqrt(dot_product(R, R))
!          if (rij <= tiny_r) cycle

!          invr    = 1.0_wp / rij
!          uvec(:) = R(:) * invr

!          ! ---- First derivatives (SD and SQ coefficient) ----
!          call self%kernel%compute_kernel_dkdr_ij(nat, xyz, brad, brdr, i, j, dKdr_ij)
!          g(:) = dKdr_ij(:, j)

!          ! SD: A_sd(:,j,i) += dK/dr_j
!          amat_sd(:, j, i) = amat_sd(:, j, i) + g(:)

!          ! ---- Second derivatives (DD and SQ coefficient) ----
!          call self%kernel%compute_kernel_d2kdr2_ij(nat, xyz, brad, brdr, brdr2, i, j, d2Kdr2_ij)

!          ! DD: A_dd(:,j,:,i) -= d2K/(dr_j dr_j)
!          amat_dd(:, j, :, i) = amat_dd(:, j, :, i) - d2Kdr2_ij(:, j, :, j)

!          ! ---- SQ (legacy-compatible radial form): Q = coef * u u^T ----
!          gpar  = dot_product(g, uvec)
!          Hu(:) = matmul(d2Kdr2_ij(:, j, :, j), uvec)
!          kpp   = dot_product(uvec, Hu)

!          coef  = (kpp + gpar * invr) / 3.0_wp

!          Q11 = coef * uvec(1) * uvec(1)
!          Q22 = coef * uvec(2) * uvec(2)
!          Q33 = coef * uvec(3) * uvec(3)
!          Q12 = coef * uvec(1) * uvec(2)
!          Q13 = coef * uvec(1) * uvec(3)
!          Q23 = coef * uvec(2) * uvec(3)

!          tc(1) = Q11
!          tc(2) = 2.0_wp * Q12
!          tc(3) = Q22
!          tc(4) = 2.0_wp * Q13
!          tc(5) = 2.0_wp * Q23
!          tc(6) = Q33

!          amat_sq(:, j, i) = amat_sq(:, j, i) + tc

!          ! ---- Third derivatives (DQ) ----
!          call self%kernel%compute_kernel_d3Kdr3_ij(nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)

!          ! U_{a,bc} = -(1/3) * d^3K / (dx_{j,a} dx_{j,b} dx_{j,c})
!          U_dq(:,:,:) = - (1.0_wp/3.0_wp) * d3Kdr3_ij(:, j, :, j, :, j)

!          ! Pack (bc) -> p in (11,12,22,13,23,33) with off-diagonals doubled
!          do a = 1, 3
!             amat_dq(a, j, 1, i) = amat_dq(a, j, 1, i) + U_dq(a,1,1)
!             amat_dq(a, j, 2, i) = amat_dq(a, j, 2, i) + 2.0_wp * U_dq(a,1,2)
!             amat_dq(a, j, 3, i) = amat_dq(a, j, 3, i) + U_dq(a,2,2)
!             amat_dq(a, j, 4, i) = amat_dq(a, j, 4, i) + 2.0_wp * U_dq(a,1,3)
!             amat_dq(a, j, 5, i) = amat_dq(a, j, 5, i) + 2.0_wp * U_dq(a,2,3)
!             amat_dq(a, j, 6, i) = amat_dq(a, j, 6, i) + U_dq(a,3,3)
!          end do

!       end do
!    end do

!    deallocate(dKdr_ij, d2Kdr2_ij, d3Kdr3_ij, brdr2, brdr3, temp, temp2)

! end subroutine get_multipole_matrix


subroutine get_multipole_matrix(self, mol, xyz, keps, brad, brdr, amat_sd, amat_dd, amat_sq, amat_dq, amat_qq)
   ! Convention-consistent multipole interaction matrices built purely from kernel derivatives.
   !
   ! This routine is "universal" in the sense that it works for any kernel K as long as the
   ! derivative providers return the correct derivatives w.r.t. coordinates of atom j.
   !
   ! SD, DD, DQ are direct derivative definitions (same as before).
   ! SQ is kept in the same legacy radial form used in your existing implicit routine:
   !    Q_bc = coef * u_b u_c,  with  coef = (kpp + gpar/r) / 3
   !
   ! QQ is constructed to match the explicit Coulomb W_abcd convention:
   !    W_abcd = (1/3) * d4K_abcd  +  (1/2) * sym(δδ)_abcd * s5
   ! with s5 = coef / r^2 (=> 1/r^5 for Coulomb), and sym(δδ)_abcd = δ_ab δ_cd + δ_ac δ_bd + δ_ad δ_bc.
   !
   ! Packing convention for symmetric pairs:
   !   p = 1..6 corresponds to (11,12,22,13,23,33) and off-diagonals are doubled in storage.

   class(alpb_solvation), intent(in) :: self
   type(structure_type), intent(in)  :: mol
   real(wp), intent(in)              :: xyz(:, :)
   real(wp), intent(in)              :: keps
   real(wp), intent(in)              :: brad(:)
   real(wp), contiguous, intent(in)  :: brdr(:, :, :)

   real(wp), contiguous, intent(inout) :: amat_sd(:, :, :)        ! (3,nat,nat)
   real(wp), intent(inout)             :: amat_dd(:, :, :, :)      ! (3,nat,3,nat)
   real(wp), intent(inout)             :: amat_sq(:, :, :)         ! (6,nat,nat)
   real(wp), intent(inout)             :: amat_dq(:, :, :, :)      ! (3,nat,6,nat)
   real(wp), intent(inout)             :: amat_qq(:, :, :, :)      ! (6,nat,6,nat)

   integer :: i, j, nat, a, b, c, d, p, q
   integer :: fab, fcd

   real(wp), parameter :: tiny_r = 1.0e-14_wp
   real(wp) :: R(3), rij, invr
   real(wp) :: uvec(3), g(3), Hu(3), gpar, kpp, coef
   real(wp) :: Q11, Q22, Q33, Q12, Q13, Q23
   real(wp) :: tc(6)

   real(wp) :: U_dq(3,3,3)

   real(wp) :: I3(3,3)
   real(wp) :: sym2, Wabcd, s5

   ! Derivative work arrays
   real(wp), allocatable :: dKdr_ij(:, :)                         ! (3,nat)
   real(wp), allocatable :: d2Kdr2_ij(:, :, :, :)                 ! (3,nat,3,nat)
   real(wp), allocatable :: d3Kdr3_ij(:, :, :, :, :, :)           ! (3,nat,3,nat,3,nat)
   real(wp), allocatable :: d4Kdr4_ij(:, :, :, :, :, :, :, :)     ! (3,nat,3,nat,3,nat,3,nat)

   ! Born radii derivative tensors (kept for interface compatibility)
   real(wp), allocatable :: brdr2(:, :, :, :, :)                  ! (3,nat,3,nat,nat)
   real(wp), allocatable :: brdr3(:, :, :, :, :, :, :)            ! (3,nat,3,nat,3,nat,nat)
   real(wp), allocatable :: brdr4(:, :, :, :, :, :, :, :, :)      ! (3,nat,3,nat,3,nat,3,nat,nat)
   real(wp), allocatable :: temp(:), temp2(:,:,:)

   integer, parameter :: pa(6) = [1, 1, 2, 1, 2, 3]
   integer, parameter :: pb(6) = [1, 2, 2, 3, 3, 3]
   integer, parameter :: pf(6) = [1, 2, 1, 2, 2, 1]

   nat = mol%nat

   allocate(dKdr_ij(3, nat), source=0.0_wp)
   allocate(d2Kdr2_ij(3, nat, 3, nat), source=0.0_wp)
   allocate(d3Kdr3_ij(3, nat, 3, nat, 3, nat), source=0.0_wp)
   allocate(d4Kdr4_ij(3, nat, 3, nat, 3, nat, 3, nat), source=0.0_wp)

   allocate(brdr2(3, nat, 3, nat, nat), source=0.0_wp)
   allocate(brdr3(3, nat, 3, nat, 3, nat, nat), source=0.0_wp)          ! frozen radii => stays zero
   allocate(brdr4(3, nat, 3, nat, 3, nat, 3, nat, nat), source=0.0_wp)  ! frozen radii => stays zero

   allocate(temp(nat), source=0.0_wp)
   allocate(temp2(3, nat, nat), source=0.0_wp)

   call self%gbobc%get_rad(mol, temp, temp2, brdr2)

   I3 = 0.0_wp
   I3(1,1)=1.0_wp; I3(2,2)=1.0_wp; I3(3,3)=1.0_wp

   do i = 1, nat
      do j = 1, nat
         if (i == j) cycle

         R(:)  = xyz(:, i) - xyz(:, j)
         rij   = sqrt(dot_product(R, R))
         if (rij <= tiny_r) cycle

         invr    = 1.0_wp / rij
         uvec(:) = R(:) * invr

         ! ---- First derivatives (SD and also used to form coef) ----
         call self%kernel%compute_kernel_dkdr_ij(nat, xyz, brad, brdr, i, j, dKdr_ij)
         g(:) = dKdr_ij(:, j)

         ! SD: A_sd(:,j,i) += dK/dr_j
         amat_sd(:, j, i) = amat_sd(:, j, i) + g(:)

         ! ---- Second derivatives (DD and used to form coef and SQ) ----
         call self%kernel%compute_kernel_d2kdr2_ij(nat, xyz, brad, brdr, brdr2, i, j, d2Kdr2_ij)

         ! DD: A_dd(:,j,:,i) -= d2K/(dr_j dr_j)
         amat_dd(:, j, :, i) = amat_dd(:, j, :, i) - d2Kdr2_ij(:, j, :, j)

         ! ---- SQ (legacy-compatible radial form): Q = coef * u u^T ----
         gpar  = dot_product(g, uvec)
         Hu(:) = matmul(d2Kdr2_ij(:, j, :, j), uvec)
         kpp   = dot_product(uvec, Hu)

         coef  = (kpp + gpar * invr) / 3.0_wp

         Q11 = coef * uvec(1) * uvec(1)
         Q22 = coef * uvec(2) * uvec(2)
         Q33 = coef * uvec(3) * uvec(3)
         Q12 = coef * uvec(1) * uvec(2)
         Q13 = coef * uvec(1) * uvec(3)
         Q23 = coef * uvec(2) * uvec(3)

         tc(1) = Q11
         tc(2) = 2.0_wp * Q12
         tc(3) = Q22
         tc(4) = 2.0_wp * Q13
         tc(5) = 2.0_wp * Q23
         tc(6) = Q33

         amat_sq(:, j, i) = amat_sq(:, j, i) + tc

         ! ---- Third derivatives (DQ) ----
         call self%kernel%compute_kernel_d3Kdr3_ij(nat, xyz, brad, brdr, brdr2, brdr3, i, j, d3Kdr3_ij)

         ! U_{a,bc} = -(1/3) * d^3K / (dx_{j,a} dx_{j,b} dx_{j,c})
         U_dq(:,:,:) = - (1.0_wp/3.0_wp) * d3Kdr3_ij(:, j, :, j, :, j)

         ! Pack (bc) -> p in (11,12,22,13,23,33) with off-diagonals doubled
         do a = 1, 3
            amat_dq(a, j, 1, i) = amat_dq(a, j, 1, i) + U_dq(a,1,1)
            amat_dq(a, j, 2, i) = amat_dq(a, j, 2, i) + 2.0_wp * U_dq(a,1,2)
            amat_dq(a, j, 3, i) = amat_dq(a, j, 3, i) + U_dq(a,2,2)
            amat_dq(a, j, 4, i) = amat_dq(a, j, 4, i) + 2.0_wp * U_dq(a,1,3)
            amat_dq(a, j, 5, i) = amat_dq(a, j, 5, i) + 2.0_wp * U_dq(a,2,3)
            amat_dq(a, j, 6, i) = amat_dq(a, j, 6, i) + U_dq(a,3,3)
         end do

         ! ---- Fourth derivatives (QQ, convention-consistent) ----
         call self%kernel%compute_kernel_d4Kdr4_ij(nat, xyz, brad, brdr, brdr2, brdr3, brdr4, i, j, d4Kdr4_ij)

         ! s5 = coef / r^2  (for Coulomb: coef=1/r^3 => s5=1/r^5)
         s5 = coef * invr * invr

         do p = 1, 6
            a   = pa(p)
            b   = pb(p)
            fab = pf(p)

            do q = 1, 6
               c   = pa(q)
               d   = pb(q)
               fcd = pf(q)

               ! sym2 = δ_ab δ_cd + δ_ac δ_bd + δ_ad δ_bc
               sym2 = I3(a,b)*I3(c,d) + I3(a,c)*I3(b,d) + I3(a,d)*I3(b,c)

               ! Match explicit Coulomb convention for any kernel:
               !   W = (1/3) * d4K_abcd + (1/2) * sym(δδ)_abcd * s5
               Wabcd = (1.0_wp/3.0_wp) * d4Kdr4_ij(a, j, b, j, c, j, d, j) &
                     + 0.5_wp * sym2 * s5

               amat_qq(p, j, q, i) = amat_qq(p, j, q, i) + real(fab*fcd, wp) * Wabcd
            end do
         end do

      end do
   end do

   deallocate(dKdr_ij, d2Kdr2_ij, d3Kdr3_ij, d4Kdr4_ij)
   deallocate(brdr2, brdr3, brdr4, temp, temp2)

end subroutine get_multipole_matrix














end module tblite_solvation_alpb
