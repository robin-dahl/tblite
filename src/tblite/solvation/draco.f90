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

!> @file tblite/solvation/draco.f90
!> Reexport of all DRACO related modules

!> Definition of an abstract base class for defining solvation models
module tblite_solvation_draco
    use tblite_solvation_draco_type, only : draco_type, new_draco
    implicit none
    private
 
    public :: draco_type, new_draco
 
 end module tblite_solvation_draco
 