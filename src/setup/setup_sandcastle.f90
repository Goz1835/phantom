!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup for a collapsing sandcastle
!
! :References: None
!
! :Owner: Jamie Lawson
!
! :Runtime parameters:
!   -
!
! :Dependencies: boundary, dim, infile_utils, io, io_control, kernel,
!   mpidomain, mpiutils, part, physcon, setup_params, timestep, unifdis,
!   units
!
 implicit none
 public :: setpart

 private
 !--private module variables
 integer :: npartb, nparts, nlayers
 real    :: boxsize, height, radius

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:maxvxyzu,gr
 use setup_params, only:rhozero
 use unifdis,      only:set_unifdis
 use io,           only:master,fatal
 use physcon,      only:pi,metre,kg,seconds,g_per_cc
 use timestep,     only:tmax,dtmax
 use kernel,       only:hfact_default
 use part,         only:igas,periodic,set_particle_type,iboundary
 use mpiutils,     only:reduceall_mpi
 use mpidomain,    only:i_belong
 use viscosity,    only:irealvisc
 use granular_variables, only:rhos
 use eos,          only:ieos
 use infile_utils, only:get_options,infile_exists
 use units,        only:set_units,udist,unit_density,utime
 use externalforces,only:iext_gravity,grav_accel
 use options,      only:iexternalforce
 use granular_variables,     only:rhos
 integer,           intent(in)    :: id
 integer,           intent(out)   :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 real,              intent(out)   :: vxyzu(:,:)
 real                             :: deltab,deltas,totmass
 integer                          :: i,ierr

 call set_units(dist=metre,time=seconds,mass=kg)
 !
 ! General parameters
 !
 time        = 0.
 hfact       = hfact_default
 rhozero     = 3.5*g_per_cc/unit_density
 rhos        = rhozero ! Set rhos in the equation of state to be the rho zero that is determined here (in code units)
 gamma       = 1
 polyk       = 0.
 iexternalforce = iext_gravity
 grav_accel  = 9.81*metre/udist/(seconds/utime)
 irealvisc = 4
 ieos = 26
 !
 ! Default setup parameters
 !
 boxsize     = 10.*metre/udist
 nparts      = 40
 npartb      = 80
 height      = 5.0*metre/udist
 radius      = 1.0*metre/udist
 nlayers     = 3
 !
 ! Infile
 !
 if (.not. infile_exists(fileprefix)) then
    tmax      = 1.2
    dtmax     = 0.005
 endif
 !
 ! Read setup parameters from setup file
 !
 print "(/,1x,63('-'),1(/,1x,a),/,1x,63('-'),/)", 'Sandcastles.'

 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
                  read_setupfile,write_setupfile)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
 deltab= 2*boxsize/npartb
 deltas = 2*boxsize/nparts
 !
 ! Put particles on grid
 !
 call set_unifdis('cubic',id,master,-0.5*boxsize,0.5*boxsize,&
                  -0.5*boxsize,0.5*boxsize,-nlayers*deltab,0.0,&
                  deltab,hfact,npart,xyzh,periodic,mask=i_belong)
 !
 ! Finalise particle properties
 !
 npartoftype(:)    = 0
 npartoftype(iboundary) = npart
 totmass           = rhozero*(boxsize*boxsize*nlayers*deltab)
 massoftype(iboundary)        = totmass/reduceall_mpi('+',npartoftype(iboundary))
 if (id==master) print*,' boundary particle mass = ',massoftype(iboundary)

 do i=1,npart
    vxyzu(:,i) = 0.
    call set_particle_type(i,iboundary)
 enddo

 !
 ! making sandcastle
 !
 call set_unifdis('closepacked',id,master,-radius,radius,&
                  -radius,radius,0.0,height,&
                  deltas,hfact,npart,xyzh,periodic,mask=i_belong, rcylmin=0.0, rcylmax=radius)

 !
 ! Finalise particle properties
 !
 npartoftype(igas) = npart-npartoftype(iboundary)
 totmass           = rhozero*(3.14159*radius**2*height)
 massoftype(igas)      = totmass/reduceall_mpi('+',npartoftype(igas))
 if (id==master) print*,' gas particle mass = ',massoftype(igas)

 do i=npartoftype(iboundary)+1,npart
   vxyzu(:,i) = 0.
   call set_particle_type(i,igas)
enddo

 rhos = rhozero

end subroutine setpart

!----------------------------------------------------------------
!+
!  write parameters to setup file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20

 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for Sandcastle setup routine'
 call write_inopt(npartb, 'npartb' ,'number of particles in x-direction for boundary',iunit)
 call write_inopt(nparts, 'nparts' ,'number of particles in x-direction for sandcastle',iunit)
 call write_inopt(boxsize,'boxsize','size of the box'   ,iunit)
 call write_inopt(height,'height','height of sandcastle',iunit)
 call write_inopt(radius,'radius','radius of sandcastle',iunit)
 call write_inopt(nlayers,'nlayers','number of particle layers in boundary',iunit)
 close(iunit)

end subroutine write_setupfile
!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)
 integer :: nerr

 nerr = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(npartb ,'npartb' ,db,min=8,errcount=nerr)
 call read_inopt(nparts ,'nparts' ,db,min=8,errcount=nerr)
 call read_inopt(boxsize,'boxsize',db,min=0.,errcount=nerr)
 call read_inopt(radius,'radius',db,min=0.,max=boxsize,errcount=nerr)
 call read_inopt(height,'height',db,min=0.,max=boxsize,errcount=nerr)
 call read_inopt(nlayers,'nlayers',db,min=1,max=100,errcount=nerr)
 call close_db(db)

 if (nerr > 0) then
    print "(1x,i2,a)",nerr,' error(s) during read of setup file: re-writing...'
    ierr = nerr
 endif

end subroutine read_setupfile

end module setup
