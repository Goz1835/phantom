!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module granular_variables
!
! Holds the variables to be read by granular.f90 and eos_incomp, and controls read/write
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - 
!
! :Dependencies: dim, eos, infile_utils, io, part, timestep, units
!
 implicit none

 real, public :: ds,mus,mu2,I0,rhos,K
 character(len=12) :: ds_in,K_in 
 
 public :: set_defaults_granular, write_options_granular, read_options_granular

 private

contains

!----------------------------------------------------------------
!+
!  set default values for the viscosity coefficients
!+
!----------------------------------------------------------------
subroutine set_defaults_granular
   use units, only:in_code_units, unit_density
   use physcon, only:g_per_cc

   integer :: ierr

   ! Units that don't call 'in_code_units' are unitless.
   ! The importance of these values is taken from Bui 2021 (eqn 64 and 65)
   ! The default values are for glass balls, taken from Jop et al., 2006
   
   ! Original values from Forterre 2003
   mus = 0.38 ! static friction coefficient, unitless
   mu2 = 0.64 ! critical friction angle at high I, unitless
   I0 = 0.279 ! a material constant (unnamed?), unitless

   ds_in = '0.053*cm' ! grain diameter   
   ds = in_code_units(ds_in, ierr, 'length')
   ierr = 0
   
   ! Set the default of rhos to be -1, as it then needs to be set by the setup
   ! A recommended value would be setting this to 2.5
   rhos = -1.0
   print*,rhos," rhos granular_variables",rhos*unit_density

   ! Following values need to be calibrated and are estimates roughly by Ha Bui
   K_in = "100.*kPa" ! elastic bulk modulus
   K = in_code_units(K_in, ierr, 'pressure')
   ierr = 0


end subroutine set_defaults_granular

!----------------------------------------------------------------
!+
!  routine to write physical granularity options to input file
!+
!----------------------------------------------------------------
subroutine write_options_granular(iwritein)
  use infile_utils, only:write_inopt
  use units,        only:unit_density
  integer, intent(in) :: iwritein
 
  write(iwritein,"(/,a)") '# options controlling granular flow model'
  call write_inopt(ds_in,'ds','grain diameter of the material (code units or e.g. 0.053*cm)',iwritein)
  call write_inopt(mus,'mus','static friction coefficient (e.g 0.38), unitless',iwritein)
  call write_inopt(mu2,'mu2','material constant for yielding (e.g 0.64), unitless',iwritein)
  call write_inopt(I0,'I0','material constant for shear (e.g 0.279), unitless',iwritein)
  call write_inopt(rhos * unit_density,'rhos','solid density of the material (in g/cm^3)',iwritein)

  call write_inopt(K_in,'K','bulk modulus (code units or e.g. 100.*kPa)',iwritein)
 
 end subroutine write_options_granular
 
 !----------------------------------------------------------------
 !+
 !  routine to read physical granularity options from input file
 !+
 !----------------------------------------------------------------
 subroutine read_options_granular(db,nerr)
  use io,           only:error
  use infile_utils, only:inopts,read_inopt
  use units,        only:in_code_units, unit_density
  type(inopts), intent(inout) :: db(:)
  integer,      intent(inout) :: nerr
  character(len=*), parameter :: label = 'read_infile'
  integer :: ierr = 0
  real :: rhos_cgs
 
  call read_inopt(ds_in,'ds',db,errcount=nerr)
  call read_inopt(rhos_cgs,'rhos',db,errcount=nerr)
  call read_inopt(K_in,'K',db,errcount=nerr)

  call read_inopt(mu2,'mu2',db,errcount=nerr,min=0.,max=1.,default=mu2)
  call read_inopt(I0,'I0',db,errcount=nerr,min=0.,max=1.,default=I0)
  call read_inopt(mus,'mus',db,errcount=nerr,min=0.,max=1.,default=mus)
 
  ds = in_code_units(ds_in, ierr, 'length')
  if (ierr /= 0) then
     call error(label,'could not convert units of length for ds')
     nerr = nerr + 1
  endif  
  ierr = 0

  ! Convert rhos into code units from the cgs version
  rhos = rhos_cgs / unit_density

  K = in_code_units(K_in, ierr, 'pressure')
  if (ierr /= 0) then
     call error(label,'could not convert units of pressure for K')
     nerr = nerr + 1
  endif  

 end subroutine read_options_granular


end module granular_variables
