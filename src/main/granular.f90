!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module granular
!
! Routines related to the μ(I) model for granular viscosity
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

 public :: granularshearfunc

 private

contains

!----------------------------------------------------------------
!+
!  function which returns shear parameter \nu as per
!  the mu I model as explained in (Bui 2021)
!+
!----------------------------------------------------------------
real function granularshearfunc(strain, granularpressure, rho1i)
 use granular_variables,    only: ds,mus,mu2,I0,rhos

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
