!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module eos_incomp
!
! Equation of state from Bui 2021 for incompressible solids
!
! :References:
!   Bui 2021
!
! Implementation from 
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - K       : *elastic bulk modulus, ???*
!
! :Dependencies: infile_utils
!
 implicit none

 public :: init_eos_incomp, equationofstate_incomp
 public :: eos_info_incomp, read_options_eos_incomp, write_options_eos_incomp

 private

contains
!-----------------------------------------------------------------------
!+
!  initialise the equation of state
!+
!-----------------------------------------------------------------------
subroutine init_eos_incomp(ierr)
 use units,         only:in_code_units
 integer, intent(out) :: ierr

 ierr = 0

end subroutine init_eos_incomp

!-----------------------------------------------------------------------
!+
!  EoS from Tillotson (1962) ; Implementation from Benz et al. (1986),
!  Melosh & Asphaug (1993) and Kegerreis et al. (2019)
!  notes:
!   - u_iv (incipient vaporisation) is called u_s  in Benz et al. 1986
!   - u_cv (complete vaporisation) is called u_s' in Benz et al. 1986
!+
!-----------------------------------------------------------------------
subroutine equationofstate_incomp(rho,pressure,spsound,gamma)
 use granular_variables,    only: rhos,K
 real, intent(inout) :: rho
 real, intent(out)   :: pressure, spsound, gamma

 real :: rhos_inv

 ! Everything has been input in cgs

 rhos_inv = 1. / rhos ! rho s is already in code units, so just need to get the inverse
 pressure = K * (rho * rhos_inv-1.)
 gamma = 1.
 spsound = sqrt(K * rhos_inv)

end subroutine equationofstate_incomp

!----------------------------------------------------------------
!+
!  print eos information
!+
!----------------------------------------------------------------
subroutine eos_info_incomp(iprint)
 use granular_variables,    only: rhos,K
 integer, intent(in) :: iprint

 write(iprint,"(/,a)") ' Incompressible EoS'
 write(iprint,"(a,1pg10.3,a)") '  rhos = ',rhos,' g/cm^3 [solid density of material]'
 write(iprint,"(a,1pg10.3,a)") '    K = ',K, ' unitless  [bulk elastic modulus of material]'

end subroutine eos_info_incomp

!-----------------------------------------------------------------------
!+
!  reads equation of state options from the input file
!+
!-----------------------------------------------------------------------
subroutine read_options_eos_incomp(db,nerr)
 use infile_utils, only:inopts,read_inopt
 use units,        only:in_code_units
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr


end subroutine read_options_eos_incomp

!-----------------------------------------------------------------------
!+
!  writes equation of state options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_eos_incomp(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

end subroutine write_options_eos_incomp

end module eos_incomp
