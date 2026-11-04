!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module granular
!
! Routines and variables related to the μ(I) model for granular viscosity.
! Holds the material parameters (also used by eos_incomp) and controls read/write
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   -
!
! :Dependencies: infile_utils, io, units
!
 implicit none

 real, public :: ds,mus,mu2,I0,rhos,K
 character(len=12) :: ds_in, K_in, rhos_in

 public :: granularshearfunc
 public :: set_defaults_granular, write_options_granular, read_options_granular, set_granular_density, set_granular_K
 public :: init_granular

 private

contains

!----------------------------------------------------------------
!+
!  set default values for the viscosity coefficients
!+
!----------------------------------------------------------------
subroutine set_defaults_granular

   ! Units that don't call 'in_code_units' are unitless.
   ! The importance of these values is taken from Bui 2021 (eqn 64 and 65)
   ! The default values are for glass balls, taken from Jop et al., 2006

   ! Original values from Forterre 2003
   mus = 0.38 ! static friction coefficient, unitless
   mu2 = 0.64 ! critical friction angle at high I, unitless
   I0 = 0.279 ! a material constant (unnamed?), unitless

   ! Dimensional values are stored as strings and only converted to code units
   ! in init_granular, once the units are known (i.e. after the dump is read)
   ds_in = '0.053*cm' ! grain diameter
   ds = -1.

   ! Set the default of rhos to be -1, as it then needs to be set by the setup
   ! A recommended value would be setting this to 2.5
   rhos_in = '-1.0*g/cm^3'
   rhos = -1.  ! negative flags "not set" to check_setup

   ! Following values need to be calibrated and are estimates roughly by Ha Bui
   K_in = '-1.0*kpa'
   K = -1.  ! negative flags "not set" to check_setup

end subroutine set_defaults_granular

!----------------------------------------------------------------
!+
!  routine to write physical granularity options to input file
!+
!----------------------------------------------------------------
subroutine write_options_granular(iwritein)
  use infile_utils, only:write_inopt
  integer, intent(in) :: iwritein

  write(iwritein,"(/,a)") '# options controlling granular flow model'
  call write_inopt(ds_in,'ds','grain diameter of the material (code units or e.g. 0.053*cm)',iwritein)
  call write_inopt(mus,'mus','static friction coefficient (e.g 0.38), unitless',iwritein)
  call write_inopt(mu2,'mu2','material constant for yielding (e.g 0.64), unitless',iwritein)
  call write_inopt(I0,'I0','material constant for shear (e.g 0.279), unitless',iwritein)

  call write_inopt(rhos_in,'rhos','solid density of the material (code units or e.g 2.5*g/cm^3)',iwritein)
  call write_inopt(K_in,'K','bulk modulus (code units or e.g 101.*kpa)',iwritein)

 end subroutine write_options_granular

 !----------------------------------------------------------------
 !+
 !  routine to read physical granularity options from input file
 !+
 !----------------------------------------------------------------
 subroutine read_options_granular(db,nerr)
  use infile_utils, only:inopts,read_inopt
  type(inopts), intent(inout) :: db(:)
  integer,      intent(inout) :: nerr

  ! Only the strings are read here: the units are not known yet when the
  ! infile is read, so the conversion to code units is done in init_granular
  call read_inopt(ds_in,'ds',db,errcount=nerr)
  call read_inopt(rhos_in,'rhos',db,errcount=nerr)
  call read_inopt(K_in,'K',db,errcount=nerr)

  call read_inopt(mu2,'mu2',db,errcount=nerr,min=0.,max=1.,default=mu2)
  call read_inopt(I0,'I0',db,errcount=nerr,min=0.,max=1.,default=I0)
  call read_inopt(mus,'mus',db,errcount=nerr,min=0.,max=1.,default=mus)

 end subroutine read_options_granular

 !----------------------------------------------------------------
 !+
 !  convert the granular parameters to code units
 !  (must be called AFTER the units are known, i.e. after the dump is read)
 !+
 !----------------------------------------------------------------
 subroutine init_granular(ierr)
  use io,           only:error
  use units,        only:in_code_units
  integer, intent(out) :: ierr
  character(len=*), parameter :: label = 'init_granular'
  integer :: ierr1

  ierr = 0

  ds = in_code_units(ds_in, ierr1, 'length')
  if (ierr1 /= 0) then
     call error(label,'could not convert units of length for ds = '//trim(ds_in))
     ierr = ierr + 1
  endif

  rhos = in_code_units(rhos_in, ierr1, 'density')
  if (ierr1 /= 0) then
     call error(label,'could not convert units of density for rhos = '//trim(rhos_in))
     ierr = ierr + 1
  endif

  K = in_code_units(K_in, ierr1, 'pressure')
  if (ierr1 /= 0) then
     call error(label,'could not convert units of pressure for K = '//trim(K_in))
     ierr = ierr + 1
  endif

 end subroutine init_granular

 !----------------------------------------------------------------
 !+
 !  routine to set the reference density
 !+
 !----------------------------------------------------------------
 subroutine set_granular_density(rhos_string)
   use units,  only:in_code_units
   character(len=*), intent(in) :: rhos_string
   integer :: ierr

   ! This subroutine is for setting rhos from outside of the file (e.g, setting it in the setup)
   ! rhos_string should be e.g 2.5*g/cm^3 so 'in_code_units' can convert it
   rhos_in = rhos_string
   rhos = in_code_units(rhos_in, ierr, 'density')
 end subroutine set_granular_density

 !----------------------------------------------------------------
 !+
 !  routine to set the K
 !+
 !----------------------------------------------------------------
 subroutine set_granular_K(K_string)
   use units,  only:in_code_units
   character(len=*), intent(in) :: K_string
   integer :: ierr

   ! This subroutine is for setting K from outside of the file (e.g, setting it in the setup)
   ! K_string should be e.g "101.*kpa" so in_code_units can convert it
   K_in = K_string
   K = in_code_units(K_string, ierr, 'pressure')
 end subroutine set_granular_K

!----------------------------------------------------------------
!+
!  function which returns shear parameter \nu as per
!  the mu I model as explained in (Bui 2021)
!+
!----------------------------------------------------------------
real function granularshearfunc(strain, granularpressure, rho1i)
 real, intent(in) :: strain(6), granularpressure, rho1i
 real :: I, muI, gammadot

 ! -- mu I model
 ! Go through how shear viscosity should be calculated in the mu I model
 ! NOTE this gives eta/rho (a value called nu in Price 2018)

 ! strain is gammadot_ij = dv_i/dx_j + dv_j/dx_i (= 2*epsilondot_ij, see strain_from_dvdx),
 ! so the shear rate is |gammadot| = sqrt(0.5*gammadot:gammadot) = sqrt(2*epsilondot:epsilondot)
 gammadot = sqrt(0.5*(strain(1)**2 + 2.*strain(2)**2 + 2.*strain(3)**2 + strain(4)**2 + 2.*strain(5)**2 + strain(6)**2))

 if(granularpressure > 0) then
   I = ds * gammadot / (granularpressure/rhos)**(0.5)

   ! If I is very low or zero, then muI --> mus
   if(I == 0) then
     muI = mus
   else
     muI = mus + (mu2 - mus) / ((I0/I) + 1)
   endif

   ! Ensure no division by zero if the strain rate is zero
   if(gammadot /= 0) then
    granularshearfunc = muI * granularpressure / gammadot * rho1i
   else
    granularshearfunc = 0
   endif

 else
   ! If pressure is negative, then the shear stress should be zero as per Bui 2021
   granularshearfunc = 0
 endif

end function granularshearfunc

end module granular
