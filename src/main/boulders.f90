!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module boulders
!
! Rigid boulders made of SPH boundary particles, following section 4.4 of
! Bui & Nguyen (2021). Each boulder k moves as a rigid body through the
! motion of its centre of mass:
!
!   M_k dV_k/dt = sum_i f_i                       (eq. 133)
!   dL_k/dt     = sum_i (r_i - R_k) x f_i         (eq. 134, L_k = I_k Omega_k)
!   dr_i/dt     = V_k + Omega_k x (r_i - R_k)     (eq. 137)
!
! f_i is the force on boundary particle i: the penalty contact force from
! regolith and other boulders (eqs. 135-136) plus gravity (self-gravity
! from the tree and external/sink forces).
!
! Boulder membership is stored per particle in iboulder (part.F90). The
! rigid-body state is rebuilt from the particles on every startup.
!
! :References:
!   Bui & Nguyen (2021), Computers and Geotechnics 138, 104315
!   Bui et al. (2014), Int. J. Numer. Anal. Meth. Geomech. 38, 1321
!
! :Owner: Daniel Price
!
! :Runtime parameters:
!   - boulder_damp          : *damping ratio of the contact dashpot*
!   - boulder_kn_fac        : *contact stiffness in units of K*dx (K = granular bulk modulus)*
!   - boulder_muw           : *interface friction coefficient mu_w*
!   - boulder_soil_reaction : *apply the contact reaction force to the regolith*
!
! :Dependencies: dim, granular_variables, infile_utils, io, part, timestep,
!   vectorutils
!
 implicit none

 public :: init_boulders,get_boulder_forces
 public :: kick_boulders,kick_boulders_ext,drift_boulders
 public :: predict_boulders,correct_boulders,set_boulder_vpred
 public :: write_options_boulders,read_options_boulders

 ! runtime parameters
 real,    public :: boulder_kn_fac = 1.
 real,    public :: boulder_damp   = 0.3
 real,    public :: boulder_muw    = 0.5
 logical, public :: boulder_soil_reaction = .true.

 private

 logical :: initialised = .false.
 integer :: nwarn_overlap = 0
 real    :: dxb      = 0.   ! initial particle spacing = contact distance (h in eq. 135)
 real    :: kcontact = 0.   ! contact stiffness K_ai
 real    :: pmassb   = 0.   ! mass of one boulder particle
 real, allocatable :: mb(:),rmaxb(:)
 real, allocatable :: xb(:,:),vb(:,:),angb(:,:),omegab(:,:)
 real, allocatable :: fb(:,:),taub(:,:),vbhalf(:,:),angbhalf(:,:)
 real, allocatable :: qb(:,:,:),ibinv(:,:,:)

contains

!----------------------------------------------------------------
!+
!  build the rigid-body state of each boulder from its particles
!  (called on every startup, so also after a restart)
!+
!----------------------------------------------------------------
subroutine init_boulders(npart,xyzh,vxyzu)
 use dim,         only:mpi,ind_timesteps,gr,use_apr
 use io,          only:fatal,warning,iprint,id,master
 use part,        only:iboulder,xyzboulder0,nboulders,massoftype,iboundary,isdead_or_accreted
 use granular_variables, only:get_granular_K_code
 use vectorutils, only:matrixinvert3D,cross_product3D
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 integer :: i,j,k,ierr
 integer :: iref(nboulders),icen(nboulders),np(nboulders)
 real    :: dxsum(3,nboulders),vsum(3,nboulders),rcen(nboulders)
 real    :: r(3),dv(3),dl(3),ib(3,3),r2,dvmax,vmax,d2,d2min,kbulk

 if (nboulders <= 0) return
 if (mpi)           call fatal('init_boulders','rigid boulders are not implemented with MPI')
 if (ind_timesteps) call fatal('init_boulders','rigid boulders require global timestepping (IND_TIMESTEPS=no)')
 if (gr)            call fatal('init_boulders','rigid boulders are not implemented with GR')
 if (use_apr)       call fatal('init_boulders','rigid boulders are not implemented with APR')

 call allocate_boulders(nboulders)
 pmassb = massoftype(iboundary)
 !
 ! centre of mass and velocity, summing offsets from a reference
 ! particle to avoid roundoff with large (e.g. heliocentric) positions
 !
 iref  = 0
 np    = 0
 dxsum = 0.
 vsum  = 0.
 do i=1,npart
    k = iboulder(i)
    if (k <= 0 .or. isdead_or_accreted(xyzh(4,i))) cycle
    if (k > nboulders) call fatal('init_boulders','boulder id exceeds nboulders',ival=k)
    if (iref(k)==0) iref(k) = i
    np(k) = np(k) + 1
    dxsum(:,k) = dxsum(:,k) + (xyzh(1:3,i) - xyzh(1:3,iref(k)))
    vsum(:,k)  = vsum(:,k) + vxyzu(1:3,i)
 enddo
 if (any(np <= 0)) call fatal('init_boulders','a boulder has no particles')

 do k=1,nboulders
    mb(k)   = np(k)*pmassb
    xb(:,k) = xyzh(1:3,iref(k)) + dxsum(:,k)/np(k)
    vb(:,k) = vsum(:,k)/np(k)
 enddo
 !
 ! body-frame offsets (Q = identity at startup), inertia tensor and
 ! angular momentum about the centre of mass
 !
 angb  = 0.
 rmaxb = 0.
 rcen  = huge(rcen)
 icen  = 0
 do k=1,nboulders
    qb(:,:,k) = 0.
    qb(1,1,k) = 1.; qb(2,2,k) = 1.; qb(3,3,k) = 1.
 enddo
 ibinv = 0.
 do i=1,npart
    k = iboulder(i)
    if (k <= 0 .or. isdead_or_accreted(xyzh(4,i))) cycle
    r = (xyzh(1:3,i) - xyzh(1:3,iref(k))) - dxsum(:,k)/np(k)
    xyzboulder0(:,i) = r
    r2 = dot_product(r,r)
    do j=1,3
       ibinv(:,j,k) = ibinv(:,j,k) - pmassb*r(:)*r(j)
       ibinv(j,j,k) = ibinv(j,j,k) + pmassb*r2
    enddo
    dv = vxyzu(1:3,i) - vb(:,k)
    call cross_product3D(r,dv,dl)
    angb(:,k) = angb(:,k) + pmassb*dl
    rmaxb(k) = max(rmaxb(k),sqrt(r2))
    if (r2 < rcen(k)) then
       rcen(k) = r2
       icen(k) = i
    endif
 enddo
 do k=1,nboulders
    ib = ibinv(:,:,k)
    call matrixinvert3D(ib,ibinv(:,:,k),ierr)
    if (ierr /= 0) then
       call warning('init_boulders','singular inertia tensor: boulder will not rotate',ival=k)
       ibinv(:,:,k) = 0.
    endif
 enddo
 call update_omega(angb,omegab)
 !
 ! project particle velocities onto the rigid motion V + Omega x r
 !
 dvmax = 0.
 vmax  = 0.
 do i=1,npart
    k = iboulder(i)
    if (k <= 0 .or. isdead_or_accreted(xyzh(4,i))) cycle
    call cross_product3D(omegab(:,k),xyzboulder0(:,i),dv)
    dv = vb(:,k) + dv
    dvmax = max(dvmax,sqrt(dot_product(dv - vxyzu(1:3,i),dv - vxyzu(1:3,i))))
    vmax  = max(vmax,sqrt(dot_product(dv,dv)))
    vxyzu(1:3,i) = dv
 enddo
 !
 ! contact distance = initial particle spacing, measured from the
 ! particle nearest the centre of each boulder (rigid lattice keeps it)
 !
 d2min = huge(d2min)
 do i=1,npart
    k = iboulder(i)
    if (k <= 0 .or. isdead_or_accreted(xyzh(4,i))) cycle
    if (i == icen(k)) cycle
    r = xyzboulder0(:,i) - xyzboulder0(:,icen(k))
    d2 = dot_product(r,r)
    if (d2 > 0.) d2min = min(d2min,d2)
 enddo
 if (d2min >= huge(d2min)) call fatal('init_boulders','could not determine particle spacing in boulders')
 dxb = sqrt(d2min)

 kbulk = get_granular_K_code()
 if (kbulk <= 0.) call fatal('init_boulders','rigid boulders need the granular bulk modulus K > 0 (ieos=26)')
 kcontact = boulder_kn_fac*kbulk*dxb
 if (kcontact <= 0.) call fatal('init_boulders','contact stiffness <= 0: check boulder_kn_fac')

 fb = 0.
 taub = 0.
 vbhalf = vb
 angbhalf = angb
 initialised = .true.

 if (id==master) then
    write(iprint,"(/,a,i4,a)") ' Rigid boulders (Bui & Nguyen 2021, sec. 4.4): ',nboulders,' boulders'
    write(iprint,"(a,es12.4,a,es12.4,a,es12.4)") ' contact distance = ',dxb,' stiffness K_ai = ',kcontact,&
                                                 ' dt_contact = ',sqrt(massoftype(iboundary)/kcontact)
    do k=1,nboulders
       write(iprint,"(a,i4,a,i7,a,es11.4,a,es11.4,a,3(es11.3),a,3(es11.3))") &
             ' boulder ',k,': np = ',np(k),' M = ',mb(k),' R = ',rmaxb(k),' V = ',vb(:,k),' Omega = ',omegab(:,k)
    enddo
    if (vmax > 0.) write(iprint,"(a,es10.3)") ' max relative velocity correction to rigid motion = ',dvmax/vmax
    write(iprint,*)
 endif

end subroutine init_boulders

!----------------------------------------------------------------
!+
!  contact forces (eqs. 135-136) between boulders and regolith and
!  between boulders, plus the total force and torque on each
!  boulder (eqs. 133-134) and the boulder timestep constraint
!+
!----------------------------------------------------------------
subroutine get_boulder_forces(npart,xyzh,vel,fxyzu,dtb)
 use part,        only:iboulder,xyzboulder0,nboulders,massoftype,igas,iphase,iamtype,&
                       isdead_or_accreted,maxphase,maxp
 use timestep,    only:C_force
 use io,          only:warning
 use vectorutils, only:cross_product3D
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:),vel(:,:)
 real,    intent(inout) :: fxyzu(:,:)
 real,    intent(out)   :: dtb
 integer :: ncell(nboulders),icelloff(nboulders),ioff(nboulders+1)
 integer :: ncs(nboulders),nbb(nboulders,nboulders)
 real    :: xmin(3,nboulders),fc(3,nboulders),tc(3,nboulders)
 integer, allocatable :: headcell(:),nextp(:),ilist(:)
 integer :: i,a,k,l,ip,ncelltot,ic,ix,iy,iz,icell(3),ipos(nboulders)
 real    :: ma,mred,r(3),dx(3),dv(3),f(3),dl(3),d2,dist,deltamax,omegamax,rcut2,dxb2
 logical :: isgas

 dtb = huge(dtb)
 if (.not.initialised .or. nboulders <= 0) return

 ma   = massoftype(igas)
 dxb2 = dxb*dxb
 !
 ! list of particles in each boulder
 !
 ioff = 0
 do i=1,npart
    k = iboulder(i)
    if (k > 0) ioff(k+1) = ioff(k+1) + 1
 enddo
 ioff(1) = 0
 do k=1,nboulders
    ioff(k+1) = ioff(k+1) + ioff(k)
 enddo
 allocate(ilist(max(ioff(nboulders+1),1)))
 ipos = ioff(1:nboulders)
 do i=1,npart
    k = iboulder(i)
    if (k > 0) then
       ipos(k) = ipos(k) + 1
       ilist(ipos(k)) = i
    endif
 enddo
 !
 ! a regular grid of cell size dx over each boulder, with a linked list
 ! of the boulder particles in each cell
 !
 ncelltot = 0
 do k=1,nboulders
    ncell(k) = max(1,ceiling(2.*(rmaxb(k) + dxb)/dxb) + 1)
    xmin(:,k) = xb(:,k) - (rmaxb(k) + dxb)
    icelloff(k) = ncelltot
    ncelltot = ncelltot + ncell(k)**3
 enddo
 allocate(headcell(ncelltot),nextp(npart))
 headcell = 0
 nextp    = 0
 do k=1,nboulders
    do ip=ioff(k)+1,ioff(k+1)
       i = ilist(ip)
       icell = min(max(int(floor((xyzh(1:3,i) - xmin(:,k))/dxb)) + 1,1),ncell(k))
       ic = icelloff(k) + icell(1) + (icell(2)-1)*ncell(k) + (icell(3)-1)*ncell(k)**2
       nextp(i) = headcell(ic)
       headcell(ic) = i
    enddo
 enddo
 !
 ! pass 1: count the contact points N_C of each boulder (eq. 135)
 !
 ncs = 0
 !$omp parallel do default(none) schedule(guided) &
 !$omp shared(npart,xyzh,iphase,iboulder,nboulders,xb,rmaxb,dxb,dxb2,xmin,ncell,icelloff,headcell,nextp) &
 !$omp shared(maxp,maxphase) &
 !$omp private(a,k,isgas,rcut2,dx,icell,ix,iy,iz,ic,i,d2) &
 !$omp reduction(+:ncs)
 do a=1,npart
    if (isdead_or_accreted(xyzh(4,a)) .or. iboulder(a) > 0) cycle
    isgas = .true.
    if (maxphase==maxp) isgas = (iamtype(iphase(a))==igas)
    if (.not.isgas) cycle
    do k=1,nboulders
       rcut2 = (rmaxb(k) + dxb)**2
       dx = xyzh(1:3,a) - xb(:,k)
       if (dot_product(dx,dx) > rcut2) cycle
       icell = int(floor((xyzh(1:3,a) - xmin(:,k))/dxb)) + 1
       do iz=max(icell(3)-1,1),min(icell(3)+1,ncell(k))
          do iy=max(icell(2)-1,1),min(icell(2)+1,ncell(k))
             do ix=max(icell(1)-1,1),min(icell(1)+1,ncell(k))
                ic = icelloff(k) + ix + (iy-1)*ncell(k) + (iz-1)*ncell(k)**2
                i = headcell(ic)
                do while (i > 0)
                   dx = xyzh(1:3,i) - xyzh(1:3,a)
                   d2 = dot_product(dx,dx)
                   if (d2 < dxb2) ncs(k) = ncs(k) + 1
                   i = nextp(i)
                enddo
             enddo
          enddo
       enddo
    enddo
 enddo
 !$omp end parallel do
 !
 ! pass 2: soil-boulder contact forces, divided by N_C
 !
 fc = 0.
 tc = 0.
 deltamax = 0.
 !$omp parallel do default(none) schedule(guided) &
 !$omp shared(npart,xyzh,vel,fxyzu,iphase,iboulder,xyzboulder0,nboulders,xb,rmaxb,mb,qb) &
 !$omp shared(dxb,dxb2,xmin,ncell,icelloff,headcell,nextp,ncs,ma,boulder_soil_reaction) &
 !$omp shared(maxp,maxphase) &
 !$omp private(a,k,isgas,rcut2,dx,dv,icell,ix,iy,iz,ic,i,d2,dist,mred,f,r,dl) &
 !$omp reduction(+:fc,tc) &
 !$omp reduction(max:deltamax)
 do a=1,npart
    if (isdead_or_accreted(xyzh(4,a)) .or. iboulder(a) > 0) cycle
    isgas = .true.
    if (maxphase==maxp) isgas = (iamtype(iphase(a))==igas)
    if (.not.isgas) cycle
    do k=1,nboulders
       if (ncs(k) <= 0) cycle
       rcut2 = (rmaxb(k) + dxb)**2
       dx = xyzh(1:3,a) - xb(:,k)
       if (dot_product(dx,dx) > rcut2) cycle
       mred = ma*mb(k)/(ma + mb(k))
       icell = int(floor((xyzh(1:3,a) - xmin(:,k))/dxb)) + 1
       do iz=max(icell(3)-1,1),min(icell(3)+1,ncell(k))
          do iy=max(icell(2)-1,1),min(icell(2)+1,ncell(k))
             do ix=max(icell(1)-1,1),min(icell(1)+1,ncell(k))
                ic = icelloff(k) + ix + (iy-1)*ncell(k) + (iz-1)*ncell(k)**2
                i = headcell(ic)
                do while (i > 0)
                   dx = xyzh(1:3,i) - xyzh(1:3,a)
                   d2 = dot_product(dx,dx)
                   if (d2 < dxb2 .and. d2 > 0.) then
                      dist = sqrt(d2)
                      deltamax = max(deltamax,dxb - dist)
                      dv = vel(1:3,i) - vel(1:3,a)
                      call contact_force(dx,dv,dist,mred,f)
                      f = f/ncs(k)
                      fc(:,k) = fc(:,k) + f
                      r = matmul(qb(:,:,k),xyzboulder0(:,i))
                      call cross_product3D(r,f,dl)
                      tc(:,k) = tc(:,k) + dl
                      if (boulder_soil_reaction) fxyzu(1:3,a) = fxyzu(1:3,a) - f/ma
                   endif
                   i = nextp(i)
                enddo
             enddo
          enddo
       enddo
    enddo
 enddo
 !$omp end parallel do
 !
 ! boulder-boulder contact forces, divided by the number of contact
 ! points between each pair so that the forces are equal and opposite
 !
 nbb = 0
 do k=1,nboulders-1
    do l=k+1,nboulders
       call boulder_pair_contacts(k,l,.false.,xyzh,vel,ioff,ilist,headcell,nextp,xmin,ncell,icelloff,&
                                  nbb(k,l),fc,tc,deltamax)
    enddo
 enddo
 do k=1,nboulders-1
    do l=k+1,nboulders
       if (nbb(k,l) > 0) call boulder_pair_contacts(k,l,.true.,xyzh,vel,ioff,ilist,headcell,nextp,xmin,ncell,&
                                                    icelloff,nbb(k,l),fc,tc,deltamax)
    enddo
 enddo
 !
 ! total force and torque on each boulder: contact forces plus gravity
 ! (self-gravity from the tree, computed in force for boulder particles)
 !
 do k=1,nboulders
    fb(:,k)   = fc(:,k)
    taub(:,k) = tc(:,k)
    do ip=ioff(k)+1,ioff(k+1)
       i = ilist(ip)
       f = pmassb*fxyzu(1:3,i)
       r = matmul(qb(:,:,k),xyzboulder0(:,i))
       call cross_product3D(r,f,dl)
       fb(:,k)   = fb(:,k) + f
       taub(:,k) = taub(:,k) + dl
    enddo
 enddo
 deallocate(headcell,nextp,ilist)
 !
 ! timestep constraint from the contact springs and the rotation rate
 !
 dtb = C_force*sqrt(ma/kcontact)
 omegamax = 0.
 do k=1,nboulders
    omegamax = max(omegamax,sqrt(dot_product(omegab(:,k),omegab(:,k))))
 enddo
 if (omegamax > 0.) dtb = min(dtb,0.1/omegamax)

 if (deltamax > 0.5*dxb .and. nwarn_overlap < 10) then
    nwarn_overlap = nwarn_overlap + 1
    call warning('boulders','contact overlap exceeds half the particle spacing: increase boulder_kn_fac',&
                 val=deltamax/dxb)
 endif

end subroutine get_boulder_forces

!----------------------------------------------------------------
!+
!  count (apply=.false.) or apply (apply=.true.) the contact forces
!  between the particles of boulders k and l
!+
!----------------------------------------------------------------
subroutine boulder_pair_contacts(k,l,apply,xyzh,vel,ioff,ilist,headcell,nextp,xmin,ncell,icelloff,&
                                 nkl,fc,tc,deltamax)
 use part,        only:xyzboulder0
 use vectorutils, only:cross_product3D
 integer, intent(in)    :: k,l
 logical, intent(in)    :: apply
 real,    intent(in)    :: xyzh(:,:),vel(:,:),xmin(:,:)
 integer, intent(in)    :: ioff(:),ilist(:),headcell(:),nextp(:),ncell(:),icelloff(:)
 integer, intent(inout) :: nkl
 real,    intent(inout) :: fc(:,:),tc(:,:),deltamax
 integer :: ip,i,j,ix,iy,iz,ic,icell(3)
 real    :: dx(3),dv(3),f(3),r(3),dl(3),d2,dist,mred,rcut2

 dx = xb(:,k) - xb(:,l)
 if (dot_product(dx,dx) > (rmaxb(k) + rmaxb(l) + dxb)**2) return
 rcut2 = (rmaxb(l) + dxb)**2
 mred  = mb(k)*mb(l)/(mb(k) + mb(l))
 do ip=ioff(k)+1,ioff(k+1)
    i = ilist(ip)
    dx = xyzh(1:3,i) - xb(:,l)
    if (dot_product(dx,dx) > rcut2) cycle
    icell = int(floor((xyzh(1:3,i) - xmin(:,l))/dxb)) + 1
    do iz=max(icell(3)-1,1),min(icell(3)+1,ncell(l))
       do iy=max(icell(2)-1,1),min(icell(2)+1,ncell(l))
          do ix=max(icell(1)-1,1),min(icell(1)+1,ncell(l))
             ic = icelloff(l) + ix + (iy-1)*ncell(l) + (iz-1)*ncell(l)**2
             j = headcell(ic)
             do while (j > 0)
                dx = xyzh(1:3,i) - xyzh(1:3,j)
                d2 = dot_product(dx,dx)
                if (d2 < dxb*dxb .and. d2 > 0.) then
                   if (apply) then
                      ! force on particle i (boulder k) from particle j (boulder l)
                      dist = sqrt(d2)
                      deltamax = max(deltamax,dxb - dist)
                      dv = vel(1:3,i) - vel(1:3,j)
                      call contact_force(dx,dv,dist,mred,f)
                      f = f/nkl
                      fc(:,k) = fc(:,k) + f
                      fc(:,l) = fc(:,l) - f
                      r = matmul(qb(:,:,k),xyzboulder0(:,i))
                      call cross_product3D(r,f,dl)
                      tc(:,k) = tc(:,k) + dl
                      r = matmul(qb(:,:,l),xyzboulder0(:,j))
                      call cross_product3D(r,-f,dl)
                      tc(:,l) = tc(:,l) + dl
                   else
                      nkl = nkl + 1
                   endif
                endif
                j = nextp(j)
             enddo
          enddo
       enddo
    enddo
 enddo

end subroutine boulder_pair_contacts

!----------------------------------------------------------------
!+
!  penalty contact force on particle i from particle a (eqs. 135-136,
!  before division by the number of contact points). dx = x_i - x_a,
!  dv = v_i - v_a and dist = |dx| < dxb. The tangential force is a
!  velocity-regularised Coulomb friction (no contact history).
!+
!----------------------------------------------------------------
subroutine contact_force(dx,dv,dist,mred,f)
 real, intent(in)  :: dx(3),dv(3),dist,mred
 real, intent(out) :: f(3)
 real :: nhat(3),vt(3),cdamp,vn,fn,ft,vtmag

 nhat  = dx/dist
 cdamp = 2.*boulder_damp*sqrt(kcontact*mred)
 vn    = dot_product(dv,nhat)
 ! normal spring-dashpot, never attractive
 fn    = max(0.,kcontact*(dxb - dist) - cdamp*vn)
 f     = fn*nhat
 ! Coulomb friction |f_t| <= mu_w |f_n|, opposing the tangential velocity
 vt    = dv - vn*nhat
 vtmag = sqrt(dot_product(vt,vt))
 if (vtmag > tiny(vtmag)) then
    ft = min(boulder_muw*fn,cdamp*vtmag)
    f  = f - ft*vt/vtmag
 endif

end subroutine contact_force

!----------------------------------------------------------------
!+
!  kick boulders using the force and torque from the last call
!  to get_boulder_forces
!+
!----------------------------------------------------------------
subroutine kick_boulders(dt,npart,vxyzu)
 use part, only:nboulders
 real,    intent(in)    :: dt
 integer, intent(in)    :: npart
 real,    intent(inout) :: vxyzu(:,:)
 integer :: k

 if (.not.initialised .or. nboulders <= 0) return
 do k=1,nboulders
    vb(:,k)   = vb(:,k) + dt*fb(:,k)/mb(k)
    angb(:,k) = angb(:,k) + dt*taub(:,k)
 enddo
 call update_omega(angb,omegab)
 call set_rigid_velocities(npart,vxyzu,vb,omegab)

end subroutine kick_boulders

!----------------------------------------------------------------
!+
!  kick boulders with the external and sink forces (fext) during
!  substepping (eqs. 133-134 for the "fast" forces)
!+
!----------------------------------------------------------------
subroutine kick_boulders_ext(dt,npart,xyzh,vxyzu,fext)
 use part,        only:iboulder,xyzboulder0,nboulders,isdead_or_accreted
 use vectorutils, only:cross_product3D
 real,    intent(in)    :: dt
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:),fext(:,:)
 real,    intent(inout) :: vxyzu(:,:)
 real    :: fk(3,nboulders),tk(3,nboulders),f(3),r(3),dl(3)
 integer :: i,k

 if (.not.initialised .or. nboulders <= 0) return
 fk = 0.
 tk = 0.
 do i=1,npart
    k = iboulder(i)
    if (k <= 0 .or. isdead_or_accreted(xyzh(4,i))) cycle
    f = pmassb*fext(1:3,i)
    r = matmul(qb(:,:,k),xyzboulder0(:,i))
    call cross_product3D(r,f,dl)
    fk(:,k) = fk(:,k) + f
    tk(:,k) = tk(:,k) + dl
 enddo
 do k=1,nboulders
    vb(:,k)   = vb(:,k) + dt*fk(:,k)/mb(k)
    angb(:,k) = angb(:,k) + dt*tk(:,k)
 enddo
 call update_omega(angb,omegab)
 call set_rigid_velocities(npart,vxyzu,vb,omegab)

end subroutine kick_boulders_ext

!----------------------------------------------------------------
!+
!  move boulders rigidly over dt (exact integral of eq. 137 for
!  constant V and Omega): translate the centre and rotate the
!  orientation, then place the particles at R + Q r0
!+
!----------------------------------------------------------------
subroutine drift_boulders(dt,npart,xyzh,vxyzu)
 use part,        only:iboulder,xyzboulder0,nboulders,isdead_or_accreted
 use vectorutils, only:rotatevec,cross_product3D
 real,    intent(in)    :: dt
 integer, intent(in)    :: npart
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 integer :: i,j,k
 real    :: w,r(3),dv(3)

 if (.not.initialised .or. nboulders <= 0) return
 do k=1,nboulders
    xb(:,k) = xb(:,k) + dt*vb(:,k)
    w = sqrt(dot_product(omegab(:,k),omegab(:,k)))
    if (w*abs(dt) > 0.) then
       do j=1,3
          call rotatevec(qb(:,j,k),omegab(:,k),w*dt)
       enddo
       call orthonormalise(qb(:,:,k))
    endif
 enddo
 call update_omega(angb,omegab)

 !$omp parallel do default(none) &
 !$omp shared(npart,xyzh,vxyzu,iboulder,xyzboulder0,xb,vb,qb,omegab) &
 !$omp private(i,k,r,dv)
 do i=1,npart
    k = iboulder(i)
    if (k <= 0) cycle
    if (isdead_or_accreted(xyzh(4,i))) cycle
    r = matmul(qb(:,:,k),xyzboulder0(:,i))
    xyzh(1:3,i) = xb(:,k) + r
    call cross_product3D(omegab(:,k),r,dv)
    vxyzu(1:3,i) = vb(:,k) + dv
 enddo
 !$omp end parallel do

end subroutine drift_boulders

!----------------------------------------------------------------
!+
!  predict the boulder velocities to the end of the step for the
!  force evaluation (stores the half-step values for the corrector)
!+
!----------------------------------------------------------------
subroutine predict_boulders(hdt,npart,vpred)
 use part, only:nboulders
 real,    intent(in)    :: hdt
 integer, intent(in)    :: npart
 real,    intent(inout) :: vpred(:,:)
 real    :: vp(3,nboulders),angp(3,nboulders),omp(3,nboulders)
 integer :: k

 if (.not.initialised .or. nboulders <= 0) return
 vbhalf   = vb
 angbhalf = angb
 do k=1,nboulders
    vp(:,k)   = vb(:,k) + hdt*fb(:,k)/mb(k)
    angp(:,k) = angb(:,k) + hdt*taub(:,k)
 enddo
 call update_omega(angp,omp)
 call set_rigid_velocities(npart,vpred,vp,omp)

end subroutine predict_boulders

!----------------------------------------------------------------
!+
!  corrector: second half kick from the half-step values, so that
!  it can be repeated if the step needs to iterate. Adds the boulder
!  particles to the velocity convergence check.
!+
!----------------------------------------------------------------
subroutine correct_boulders(hdt,npart,xyzh,vxyzu,vpred,errmax,v2mean,np)
 use part, only:iboulder,nboulders,isdead_or_accreted
 real,    intent(in)    :: hdt
 integer, intent(in)    :: npart
 real,    intent(in)    :: xyzh(:,:),vpred(:,:)
 real,    intent(inout) :: vxyzu(:,:),errmax,v2mean
 integer, intent(inout) :: np
 integer :: i,k

 if (.not.initialised .or. nboulders <= 0) return
 do k=1,nboulders
    vb(:,k)   = vbhalf(:,k) + hdt*fb(:,k)/mb(k)
    angb(:,k) = angbhalf(:,k) + hdt*taub(:,k)
 enddo
 call update_omega(angb,omegab)
 call set_rigid_velocities(npart,vxyzu,vb,omegab)

 do i=1,npart
    if (iboulder(i) <= 0) cycle
    if (isdead_or_accreted(xyzh(4,i))) cycle
    errmax = max(errmax,dot_product(vxyzu(1:3,i) - vpred(1:3,i),vxyzu(1:3,i) - vpred(1:3,i)))
    v2mean = v2mean + dot_product(vxyzu(1:3,i),vxyzu(1:3,i))
    np     = np + 1
 enddo

end subroutine correct_boulders

!----------------------------------------------------------------
!+
!  copy the corrected boulder velocities into vpred when the
!  step iterates
!+
!----------------------------------------------------------------
subroutine set_boulder_vpred(npart,vxyzu,vpred)
 use part, only:iboulder,nboulders
 integer, intent(in)    :: npart
 real,    intent(in)    :: vxyzu(:,:)
 real,    intent(inout) :: vpred(:,:)
 integer :: i

 if (.not.initialised .or. nboulders <= 0) return
 do i=1,npart
    if (iboulder(i) > 0) vpred(1:3,i) = vxyzu(1:3,i)
 enddo

end subroutine set_boulder_vpred

!----------------------------------------------------------------
!+
!  set the velocity of each boulder particle to V + Omega x r
!+
!----------------------------------------------------------------
subroutine set_rigid_velocities(npart,vel,vk,omegak)
 use part,        only:iboulder,xyzboulder0
 use vectorutils, only:cross_product3D
 integer, intent(in)    :: npart
 real,    intent(inout) :: vel(:,:)
 real,    intent(in)    :: vk(:,:),omegak(:,:)
 integer :: i,k
 real    :: r(3),dv(3)

 !$omp parallel do default(none) &
 !$omp shared(npart,vel,vk,omegak,iboulder,xyzboulder0,qb) &
 !$omp private(i,k,r,dv)
 do i=1,npart
    k = iboulder(i)
    if (k <= 0) cycle
    r = matmul(qb(:,:,k),xyzboulder0(:,i))
    call cross_product3D(omegak(:,k),r,dv)
    vel(1:3,i) = vk(:,k) + dv
 enddo
 !$omp end parallel do

end subroutine set_rigid_velocities

!----------------------------------------------------------------
!+
!  angular velocity from angular momentum, Omega = Q I_body^-1 Q^T L
!+
!----------------------------------------------------------------
subroutine update_omega(angk,omegak)
 use part, only:nboulders
 real, intent(in)  :: angk(:,:)
 real, intent(out) :: omegak(:,:)
 integer :: k

 do k=1,nboulders
    omegak(:,k) = matmul(qb(:,:,k),matmul(ibinv(:,:,k),matmul(transpose(qb(:,:,k)),angk(:,k))))
 enddo

end subroutine update_omega

!----------------------------------------------------------------
!+
!  re-orthonormalise a rotation matrix (columns) by Gram-Schmidt
!+
!----------------------------------------------------------------
subroutine orthonormalise(q)
 use vectorutils, only:cross_product3D
 real, intent(inout) :: q(3,3)

 q(:,1) = q(:,1)/sqrt(dot_product(q(:,1),q(:,1)))
 q(:,2) = q(:,2) - dot_product(q(:,1),q(:,2))*q(:,1)
 q(:,2) = q(:,2)/sqrt(dot_product(q(:,2),q(:,2)))
 call cross_product3D(q(:,1),q(:,2),q(:,3))

end subroutine orthonormalise

!----------------------------------------------------------------
!+
!  allocate the per-boulder arrays
!+
!----------------------------------------------------------------
subroutine allocate_boulders(nb)
 integer, intent(in) :: nb

 if (allocated(mb)) deallocate(mb,rmaxb,xb,vb,angb,omegab,fb,taub,vbhalf,angbhalf,qb,ibinv)
 allocate(mb(nb),rmaxb(nb),xb(3,nb),vb(3,nb),angb(3,nb),omegab(3,nb))
 allocate(fb(3,nb),taub(3,nb),vbhalf(3,nb),angbhalf(3,nb),qb(3,3,nb),ibinv(3,3,nb))

end subroutine allocate_boulders

!----------------------------------------------------------------
!+
!  write runtime options to the input file
!+
!----------------------------------------------------------------
subroutine write_options_boulders(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling rigid boulders'
 call write_inopt(boulder_kn_fac,'boulder_kn_fac','contact stiffness in units of K*dx (K = granular bulk modulus)',iunit)
 call write_inopt(boulder_damp,'boulder_damp','damping ratio of the contact dashpot',iunit)
 call write_inopt(boulder_muw,'boulder_muw','interface friction coefficient mu_w',iunit)
 call write_inopt(boulder_soil_reaction,'boulder_soil_reaction','apply the contact reaction force to the regolith',iunit)

end subroutine write_options_boulders

!----------------------------------------------------------------
!+
!  read runtime options from the input file (all optional)
!+
!----------------------------------------------------------------
subroutine read_options_boulders(db,nerr)
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(boulder_kn_fac,'boulder_kn_fac',db,errcount=nerr,min=0.,default=boulder_kn_fac)
 call read_inopt(boulder_damp,'boulder_damp',db,errcount=nerr,min=0.,default=boulder_damp)
 call read_inopt(boulder_muw,'boulder_muw',db,errcount=nerr,min=0.,default=boulder_muw)
 call read_inopt(boulder_soil_reaction,'boulder_soil_reaction',db,errcount=nerr,default=boulder_soil_reaction)

end subroutine read_options_boulders

end module boulders
