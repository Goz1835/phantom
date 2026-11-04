!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module viscosity
!
! Routines related to physical viscosity
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - bulkvisc   : *magnitude of bulk viscosity*
!   - irealvisc  : *physical viscosity type (0=none,1=const,2=Shakura/Sunyaev, 4=mu I model)*
!   - shearparam : *magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)*
!
! :Dependencies: dim, eos, infile_utils, io, part, timestep
!
 implicit none
 integer, public :: irealvisc
 real, public :: shearparam, bulkvisc, HoverR

 public :: shearfunc,set_defaults_viscosity,dt_viscosity,viscinfo
 public :: write_options_viscosity,read_options_viscosity

 private

contains

!----------------------------------------------------------------
!+
!  set default values for the viscosity coefficients
!+
!----------------------------------------------------------------
subroutine set_defaults_viscosity
 use granular,     only:set_defaults_granular

 irealvisc = 0   ! Physical viscosity
 shearparam = 0.1  ! alphadisc (if irealvisc=2) or nu if irealvisc=1
 bulkvisc = 0.0   ! bulk viscosity parameter in code units

 call set_defaults_granular

end subroutine set_defaults_viscosity

!----------------------------------------------------------------
!+
!  function which returns shear parameter \nu given
!  the alpha disc parameter, positions and sound speed
!  ie. \nu = alpha*cs*H
!
!  for irealvisc = 2 alpha is specified in the input file
!  for irealvisc = 1 \nu is specified directly in the input file
!+
!----------------------------------------------------------------
real function shearfunc(xi,yi,zi,spsoundi, ignore_sandcastle, pri, strain, rho1i)
 use part,   only:xyzmh_ptmass
 use eos,    only:polyk,qfacdisc
 use units,  only:in_code_units
 use granular, only:granularshearfunc
 real, intent(in) :: xi,yi,zi,spsoundi
 integer, intent(in) :: ignore_sandcastle
 real, intent(in) :: pri, strain(6), rho1i
 real :: rsph2,omega1,H,r1,r2

 select case(irealvisc)
 case(2)
!
!--Shakura & Sunyaev disc viscosity
!  nu = alpha*cs*H
!
!  but H = cs/Omega, so nu = alpha*cs**2/Omega
!
!  (Omega is assumed to be Keplerian ie. r**-3/2)
!
    rsph2 = xi*xi + yi*yi + zi*zi
    omega1 = rsph2**0.75  ! omega1 = 1/omega
    shearfunc = shearparam*omega1*spsoundi**2

 case(1)
!
!--constant shear viscosity
!
    shearfunc = shearparam

 case(0)
!
!--default is zero
!
   shearfunc = 0

 case(3)
!
!--Shakura & Sunyaev for ieos=14
!

    r1=sqrt((xi-xyzmh_ptmass(1,1))**2+(yi-xyzmh_ptmass(2,1))**2 + (zi-xyzmh_ptmass(3,1))**2)
    r2=sqrt((xi-xyzmh_ptmass(1,2))**2+(yi-xyzmh_ptmass(2,2))**2 + (zi-xyzmh_ptmass(3,2))**2)

    if (r2<r1) then
       H=sqrt(polyk)*r2**(-qfacdisc+1.5)
    else
       H=sqrt(polyk)*r1**(-qfacdisc+1.5)
    endif

    shearfunc=shearparam*spsoundi*H
 case(4)
!
! -- mu I model
! Go through how shear viscosity should be calculated in the mu I model
   if(ignore_sandcastle == 0) then
    shearfunc = granularshearfunc(strain, pri, rho1i)
    !shearfunc = 0
   else
    shearfunc = 0
   endif


 case default

    stop 'invalid choice for physical viscosity'

 end select

end function shearfunc

!----------------------------------------------------------------
!+
!  set explicit timestep for physical viscosity
!+
!----------------------------------------------------------------
real function dt_viscosity(xi,yi,zi,hi,spsoundi)
 use timestep,only: C_force,bignumber
 real, intent(in) :: xi,yi,zi,hi,spsoundi
 real             :: viscnu

 viscnu = shearfunc(xi,yi,zi,spsoundi, 1, 1.1, [1.1,1.1,1.1,1.1,1.1,1.1], 1.1)

 if (viscnu > tiny(viscnu)) then
    dt_viscosity = 0.4*C_force*hi*hi/viscnu
 else
    dt_viscosity = bignumber
 endif

end function dt_viscosity

!----------------------------------------------------------------
!+
!  prints info about physical viscosity settings into the header
!+
!----------------------------------------------------------------
subroutine viscinfo(ivisc,iprint)
 use dim, only:maxp,maxdvdx
 integer, intent(in) :: ivisc,iprint

 select case(ivisc)
 case(4)
   write(iprint,"(a,es10.3)") ' μ(I) model = ',shearparam
 case(3)
    write(iprint,"(a,es10.3)") ' Shakura-Sunyaev viscosity for binary systems, alpha_SS = ',shearparam
 case(2)
    write(iprint,"(a,es10.3)") ' Shakura-Sunyaev viscosity, alpha_SS = ',shearparam
 case(1)
    write(iprint,"(a,es10.3)") ' Constant physical viscosity, nu = ',shearparam
 case(0)
 case default
    write(iprint,"(a,es10.3)") ' Unknown setting for physical viscosity, nu = ',shearparam
 end select
 if (ivisc /= 0) then
    if (maxdvdx==maxp) then
       write(iprint,"(a,/)") ' (computed using two first derivatives)'
    else
       write(iprint,"(a,/)") ' (computed using direct second derivatives)'
    endif
 endif

end subroutine viscinfo

!----------------------------------------------------------------
!+
!  routine to write physical viscosity options to input file
!+
!----------------------------------------------------------------
subroutine write_options_viscosity(iwritein)
 use infile_utils, only:write_inopt
 use granular,     only:write_options_granular
 integer, intent(in) :: iwritein

 write(iwritein,"(/,a)") '# options controlling physical viscosity'
 call write_inopt(irealvisc,'irealvisc','physical viscosity type (0=none,1=const,2=Shakura/Sunyaev)',iwritein)
 call write_inopt(shearparam,'shearparam','magnitude of shear viscosity (irealvisc=1) or alpha_SS (irealvisc=2)',iwritein)
 call write_inopt(bulkvisc,'bulkvisc','magnitude of bulk viscosity',iwritein)

 if(irealvisc == 4) call write_options_granular(iwritein)

end subroutine write_options_viscosity

!----------------------------------------------------------------
!+
!  routine to read physical viscosity options from input file
!+
!----------------------------------------------------------------
subroutine read_options_viscosity(db,nerr)
 use io,           only:error
 use infile_utils, only:inopts,read_inopt
 use granular,     only:read_options_granular
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 character(len=*), parameter :: label = 'read_infile'

 call read_inopt(irealvisc,'irealvisc',db,errcount=nerr,min=0,max=12,default=0)
 call read_inopt(shearparam,'shearparam',db,errcount=nerr,min=0.,default=0.1)
 call read_inopt(bulkvisc,'bulkvisc',db,errcount=nerr,default=0.0,min=0.)
 if (irealvisc==2 .and. shearparam > 1) call error(label,'alpha > 1 for shakura-sunyaev viscosity')

 if(irealvisc == 4) call read_options_granular(db, nerr)

end subroutine read_options_viscosity

end module viscosity
