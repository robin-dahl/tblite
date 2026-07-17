! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later

!> Native tblite integral implementation.
module tblite_integral_native
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type, cgto_type
   use tblite_integral_dipole, only : dipole_3d, dipole_cgto_diat
   use tblite_integral_multipole, only : maxl2, msao, mlao, lmap, smap, sqrtpi3, lx, &
      & overlap_1d, multipole_3d, multipole_grad_3d, shift_operator, &
      & multipole_cgto_diat
   use tblite_integral_trafo, only : transform0, transform1, transform2
   use tblite_integral_handler, only : integral_handler
   implicit none
   private

   public :: native_integral_type
   public :: dipole_cgto, multipole_cgto, multipole_grad_cgto
   public :: get_dipole_integrals, get_multipole_integrals

   interface get_dipole_integrals
      module procedure :: get_dipole_integrals_lat
      module procedure :: get_dipole_integrals_diat_lat
   end interface get_dipole_integrals

   interface get_multipole_integrals
      module procedure :: get_multipole_integrals_lat
      module procedure :: get_multipole_integrals_diat_lat
   end interface get_multipole_integrals

   type, extends(integral_handler) :: native_integral_type
   contains
      procedure :: initialize_integral => initialize_native
      procedure :: multipole_cgto => multipole_native
      procedure :: multipole_grad_cgto => multipole_gradient_native
      procedure :: dipole_cgto => dipole_native
   end type native_integral_type

contains

pure subroutine dipole_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), pre
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   d3d(:, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 1
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call dipole_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3d(:, mlj, mli) = d3d(:, mlj, mli) + cc*dip
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, d3d, dpint, .true., .true.)

end subroutine dipole_cgto

pure subroutine multipole_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, qpint)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integrals for the given pair i  and j
   real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), quad(6), pre, tr
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: q3d(6, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   d3d(:, :, :) = 0.0_wp
   q3d(:, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 2
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call multipole_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip, quad)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3d(:, mlj, mli) = d3d(:, mlj, mli) + cc*dip
               q3d(:, mlj, mli) = q3d(:, mlj, mli) + cc*quad
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, d3d, dpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, q3d, qpint, .true., .true.)

   ! remove trace from quadrupole integrals (transfrom to spherical harmonics and back)
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         tr = 0.5_wp * (qpint(1, mlj, mli) + qpint(3, mlj, mli) + qpint(6, mlj, mli))
         qpint(1, mlj, mli) = 1.5_wp * qpint(1, mlj, mli) - tr
         qpint(2, mlj, mli) = 1.5_wp * qpint(2, mlj, mli)
         qpint(3, mlj, mli) = 1.5_wp * qpint(3, mlj, mli) - tr
         qpint(4, mlj, mli) = 1.5_wp * qpint(4, mlj, mli)
         qpint(5, mlj, mli) = 1.5_wp * qpint(5, mlj, mli)
         qpint(6, mlj, mli) = 1.5_wp * qpint(6, mlj, mli) - tr
      end do
   end do

end subroutine multipole_cgto

pure subroutine multipole_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, qpint, &
      & doverlap, ddpintj, dqpintj, ddpinti, dqpinti)
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integrals for the given pair i  and j
   real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpinti(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: dqpinti(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpintj(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: dqpintj(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), quad(6)
   real(wp) :: pre, grad(3), ddip(3, 3), dquad(3, 6), tr, dtr(3)
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3dj(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: q3dj(6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3di(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dq3di(3, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3dj(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dq3dj(3, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   d3dj(:, :, :) = 0.0_wp
   q3dj(:, :, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp
   dd3di(:, :, :, :) = 0.0_wp
   dq3di(:, :, :, :) = 0.0_wp
   dd3dj(:, :, :, :) = 0.0_wp
   dq3dj(:, :, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 3
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call multipole_grad_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip, quad, grad, ddip, dquad)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3dj(:, mlj, mli) = d3dj(:, mlj, mli) + cc*dip
               q3dj(:, mlj, mli) = q3dj(:, mlj, mli) + cc*quad
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad
               dd3dj(:, :, mlj, mli) = dd3dj(:, :, mlj, mli) + cc*ddip
               dq3dj(:, :, mlj, mli) = dq3dj(:, :, mlj, mli) + cc*dquad
            end do
         end do
      end do
   end do

   do mli = 1, mlao(cgtoi%ang)
      do mlj = 1, mlao(cgtoj%ang)
         call shift_operator(vec, s3d(mlj, mli), d3dj(:, mlj, mli), q3dj(:, mlj, mli), &
            & ds3d(:, mlj, mli), dd3dj(:, :, mlj, mli), dq3dj(:, :, mlj, mli), &
            & dd3di(:, :, mlj, mli), dq3di(:, :, mlj, mli))
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, d3dj, dpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, q3dj, qpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dd3dj, ddpintj, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dq3dj, dqpintj, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dd3di, ddpinti, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dq3di, dqpinti, .true., .true.)

   ! remove trace from quadrupole integrals (transfrom to spherical harmonics and back)
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         tr = 0.5_wp * (qpint(1, mlj, mli) + qpint(3, mlj, mli) + qpint(6, mlj, mli))
         qpint(1, mlj, mli) = 1.5_wp * qpint(1, mlj, mli) - tr
         qpint(2, mlj, mli) = 1.5_wp * qpint(2, mlj, mli)
         qpint(3, mlj, mli) = 1.5_wp * qpint(3, mlj, mli) - tr
         qpint(4, mlj, mli) = 1.5_wp * qpint(4, mlj, mli)
         qpint(5, mlj, mli) = 1.5_wp * qpint(5, mlj, mli)
         qpint(6, mlj, mli) = 1.5_wp * qpint(6, mlj, mli) - tr
      end do
   end do

   ! remove trace from quadrupole integral derivatives
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         dtr = dqpinti(:, 1, mlj, mli) + dqpinti(:, 3, mlj, mli) + dqpinti(:, 6, mlj, mli)
         dqpinti(:, 1, mlj, mli) = 1.5_wp * dqpinti(:, 1, mlj, mli) - 0.5_wp * dtr
         dqpinti(:, 2, mlj, mli) = 1.5_wp * dqpinti(:, 2, mlj, mli)
         dqpinti(:, 3, mlj, mli) = 1.5_wp * dqpinti(:, 3, mlj, mli) - 0.5_wp * dtr
         dqpinti(:, 4, mlj, mli) = 1.5_wp * dqpinti(:, 4, mlj, mli)
         dqpinti(:, 5, mlj, mli) = 1.5_wp * dqpinti(:, 5, mlj, mli)
         dqpinti(:, 6, mlj, mli) = 1.5_wp * dqpinti(:, 6, mlj, mli) - 0.5_wp * dtr
      end do
   end do

   ! remove trace from quadrupole integral derivatives
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         dtr = dqpintj(:, 1, mlj, mli) + dqpintj(:, 3, mlj, mli) + dqpintj(:, 6, mlj, mli)
         dqpintj(:, 1, mlj, mli) = 1.5_wp * dqpintj(:, 1, mlj, mli) - 0.5_wp * dtr
         dqpintj(:, 2, mlj, mli) = 1.5_wp * dqpintj(:, 2, mlj, mli)
         dqpintj(:, 3, mlj, mli) = 1.5_wp * dqpintj(:, 3, mlj, mli) - 0.5_wp * dtr
         dqpintj(:, 4, mlj, mli) = 1.5_wp * dqpintj(:, 4, mlj, mli)
         dqpintj(:, 5, mlj, mli) = 1.5_wp * dqpintj(:, 5, mlj, mli)
         dqpintj(:, 6, mlj, mli) = 1.5_wp * dqpintj(:, 6, mlj, mli) - 0.5_wp * dtr
      end do
   end do

end subroutine multipole_grad_cgto



!> Evaluate overlap for a molecular structure
subroutine get_dipole_integrals_lat(mol, trans, cutoff, bas, overlap, dpint)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), dtmp(:, :)

   overlap(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmp(3, msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap, dpint) private(r2, vec, stmp, dtmp) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call dipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))

                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_dipole_integrals_lat

!> Evaluate dipole integrals for a molecular structure
!> with diatomic-frame-scaled overlap elements
subroutine get_dipole_integrals_diat_lat(mol, &
   & trans, cutoff, bas, scal_fac, overlap, overlap_diat, dpint)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Scaling factors for the diatomic frame
   real(wp), intent(in) :: scal_fac(:,:)
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Overlap matrix with diatomic frame-scaled elements
   real(wp), intent(out) :: overlap_diat(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)
   !> Scaling factors for the diatomic frame for the three differnt bonding motifs
   !> (sigma, pi, delta)
   real(wp) :: ksig, kpi, kdel

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), dtmp(:, :), stmp_diat(:)

   if (size(scal_fac,1) /= 3) then
      error stop 'Error: scal_fac must have the dimension of 3, &
      & since it covers the three different types of bonding'
   end if

   overlap(:, :) = 0.0_wp
   overlap_diat(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), stmp_diat(msao(bas%maxl)**2), &
   & dtmp(3, msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap, overlap_diat, dpint, scal_fac) &
   !$omp private(r2, vec, stmp, dtmp, stmp_diat) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao) &
   !$omp private(ksig, kpi, kdel)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            if (iat /= jat) then
               !> Determine scaling factors from atom parameters
               ksig = 2.0_wp / (1.0_wp / scal_fac(1,mol%num(mol%id(iat))) &
               & + 1.0_wp / scal_fac(1,mol%num(mol%id(jat))) )
               kpi = 2.0_wp / (1.0_wp / scal_fac(2,mol%num(mol%id(iat))) &
               & + 1.0_wp / scal_fac(2,mol%num(mol%id(jat))) )
               kdel = 2.0_wp / (1.0_wp / scal_fac(3,mol%num(mol%id(iat))) &
               & + 1.0_wp / scal_fac(3,mol%num(mol%id(jat))) )
            end if
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  stmp = 0.0_wp
                  stmp_diat = 0.0_wp
                  dtmp = 0.0_wp
                  if (iat /= jat) then
                     call dipole_cgto_diat(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                        & r2, vec, bas%intcut, ksig, kpi, kdel, &
                        & stmp, stmp_diat, dtmp)
                  else
                     call dipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                        & r2, vec, bas%intcut, stmp, dtmp)
                     stmp_diat = stmp
                  endif
                  
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, msao(bas%cgto(jsh, jzp)%ang)
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))
                        
                        overlap_diat(jj+jao, ii+iao) = overlap_diat(jj+jao, ii+iao) &
                           & + stmp_diat(jao + nao*(iao-1))
                        
                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_dipole_integrals_diat_lat


!> Evaluate overlap for a molecular structure
subroutine get_multipole_integrals_lat(mol, trans, cutoff, bas, overlap, dpint, qpint)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)
   !> Quadrupole moment integral matrix
   real(wp), intent(out) :: qpint(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)

   overlap(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp
   qpint(:, :, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmp(3, msao(bas%maxl)**2), qtmp(6, msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap, dpint, qpint) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao) &
   !$omp private(r2, vec, stmp, dtmp, qtmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp, qtmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))

                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, jao + nao*(iao-1))

                        qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                           & + qtmp(:, jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_multipole_integrals_lat

!> Evaluate multipole integrals for a molecular structure
!> with diatomic-frame-scaled overlap elements
subroutine get_multipole_integrals_diat_lat(mol, & 
   & trans, cutoff, bas, scal_fac, overlap, overlap_diat, &
   & dpint, qpint)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Scaling factors for the diatomic frame
   real(wp), intent(in) :: scal_fac(:,:)
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Overlap matrix with diatomic frame-scaled elements
   real(wp), intent(out) :: overlap_diat(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)
   !> Quadrupole moment integral matrix
   real(wp), intent(out) :: qpint(:, :, :)

   !> Scaling factors for the diatomic frame for sigma-, pi-, delta-bonding
   real(wp) :: ksig, kpi, kdel

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :), stmp_diat(:)

   if (size(scal_fac,1) /= 3) then
      error stop 'Error: scal_fac must have the dimension of 3, &
      & since it covers the three different types of bonding'
   end if

   overlap(:, :) = 0.0_wp
   overlap_diat(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp
   qpint(:, :, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmp(3, msao(bas%maxl)**2), qtmp(6, msao(bas%maxl)**2), &
      & stmp_diat(msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap, overlap_diat, dpint, qpint, scal_fac) &
   !$omp private(r2, vec, stmp, dtmp, qtmp, stmp_diat) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao) &
   !$omp private(ksig, kpi, kdel)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            if (iat /= jat) then
               !> Determine scaling factors from atom parameters
               ksig = 2.0_wp / (1.0_wp / scal_fac(1,mol%num(mol%id(iat))) &
               & + 1.0_wp / scal_fac(1,mol%num(mol%id(jat))) )
               kpi = 2.0_wp / (1.0_wp / scal_fac(2,mol%num(mol%id(iat))) &
               & + 1.0_wp / scal_fac(2,mol%num(mol%id(jat))) )
               kdel = 2.0_wp / (1.0_wp / scal_fac(3,mol%num(mol%id(iat))) &
               & + 1.0_wp / scal_fac(3,mol%num(mol%id(jat))) )
            end if
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  stmp = 0.0_wp
                  stmp_diat = 0.0_wp
                  dtmp = 0.0_wp
                  qtmp = 0.0_wp
                  if (iat /= jat) then
                     call multipole_cgto_diat(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                        & r2, vec, bas%intcut, ksig, kpi, kdel, stmp, &
                        & stmp_diat, dtmp, qtmp)
                  else
                     call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                        & r2, vec, bas%intcut, stmp, dtmp, qtmp)
                     stmp_diat = stmp
                  endif

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))

                        overlap_diat(jj+jao, ii+iao) = overlap_diat(jj+jao, ii+iao) &
                           & + stmp_diat(jao + nao*(iao-1))

                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, jao + nao*(iao-1))

                        qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                           & + qtmp(:, jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_multipole_integrals_diat_lat

subroutine initialize_native(self, mol, basis)
   !> Native integral evaluator
   class(native_integral_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: basis
end subroutine initialize_native

subroutine multipole_native(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, dipole, quadrupole)
   !> Native integral evaluator
   class(native_integral_type), intent(in) :: self
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Global shell index of the contracted Gaussian function on center j
   integer, intent(in) :: jsh
   !> Global shell index of the contracted Gaussian function on center i
   integer, intent(in) :: ish
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i and j
   real(wp), intent(out) :: overlap(:)
   !> Dipole moment integrals for the given pair i and j
   real(wp), intent(out) :: dipole(:, :)
   !> Quadrupole moment integrals for the given pair i and j
   real(wp), intent(out) :: quadrupole(:, :)
   
   call multipole_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dipole, quadrupole)
end subroutine multipole_native

subroutine multipole_gradient_native(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, &
      & dipole, quadrupole, doverlap, ddipole_j, dquadrupole_j, ddipole_i, &
      & dquadrupole_i)
   !> Native integral evaluator
   class(native_integral_type), intent(in) :: self
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Global shell index of the contracted Gaussian function on center j
   integer, intent(in) :: jsh
   !> Global shell index of the contracted Gaussian function on center i
   integer, intent(in) :: ish
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i and j
   real(wp), intent(out) :: overlap(:)
   !> Dipole moment integrals for the given pair i and j
   real(wp), intent(out) :: dipole(:, :)
   !> Quadrupole moment integrals for the given pair i and j
   real(wp), intent(out) :: quadrupole(:, :)
   !> Overlap integral gradient for the given pair i and j
   real(wp), intent(out) :: doverlap(:, :)
   !> Dipole moment integral gradient with respect to center j
   real(wp), intent(out) :: ddipole_j(:, :, :)
   !> Quadrupole moment integral gradient with respect to center j
   real(wp), intent(out) :: dquadrupole_j(:, :, :)
   !> Dipole moment integral gradient with respect to center i
   real(wp), intent(out) :: ddipole_i(:, :, :)
   !> Quadrupole moment integral gradient with respect to center i
   real(wp), intent(out) :: dquadrupole_i(:, :, :)
   
   call multipole_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dipole, quadrupole, doverlap, &
      & ddipole_j, dquadrupole_j, ddipole_i, dquadrupole_i)
end subroutine multipole_gradient_native

subroutine dipole_native(self, cgtoj, cgtoi, jsh, ish, r2, vec, intcut, overlap, dipole)
   !> Native integral evaluator
   class(native_integral_type), intent(in) :: self
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Global shell index of the contracted Gaussian function on center j
   integer, intent(in) :: jsh
   !> Global shell index of the contracted Gaussian function on center i
   integer, intent(in) :: ish
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i and j
   real(wp), intent(out) :: overlap(:)
   !> Dipole moment integrals for the given pair i and j
   real(wp), intent(out) :: dipole(:, :)
   
   call dipole_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dipole)
end subroutine dipole_native

end module tblite_integral_native
