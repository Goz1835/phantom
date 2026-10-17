!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setboulders
!
! Place spherical rigid boulders inside a body made of SPH particles on a
! closepacked lattice, following section 4.4 of Bui & Nguyen (2021). The
! existing particles inside each sphere are converted into boundary
! particles of that boulder; no particles are created or removed.
!
!  - one boulder is placed at the centroid of the body
!  - several boulders are spread evenly through the body (centroidal
!    Voronoi tessellation by Lloyd iterations, seeded by farthest-point
!    sampling)
!  - a boulder centre is only accepted if the whole sphere lies inside
!    the body with at least one layer of regolith around it, and boulders
!    are kept at least two particle spacings apart, so boulders are never
!    partially generated and never overlap
!
! :References: Bui & Nguyen (2021), Computers and Geotechnics 138, 104315
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, part
!
 implicit none

 public :: set_boulders

 private

 ! a regular grid with a linked list of the particles in each cell
 type cellgrid
    real :: xmin(3),size
    integer :: n(3)
    integer, allocatable :: head(:),next(:)
 end type cellgrid

 integer, parameter :: maxlloyd = 200

contains

!----------------------------------------------------------------
!+
!  convert particles i1..i2 of a body into nb spherical boulders of
!  radius rb. The particles must come from a closepacked lattice with
!  h = hfact*dx (e.g. from set_shape)
!+
!----------------------------------------------------------------
subroutine set_boulders(nb,rb,i1,i2,xyzh,npartoftype,massoftype,hfact,ierr)
 use part, only:iboulder,nboulders,igas,iboundary,set_particle_type
 use io,   only:fatal,warning,id,master
 integer, intent(in)    :: nb,i1,i2
 real,    intent(in)    :: rb,hfact
 real,    intent(in)    :: xyzh(:,:)
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 integer, intent(out)   :: ierr
 type(cellgrid) :: fine,coarse
 integer, allocatable :: isurf(:),icand(:),nconv(:)
 real,    allocatable :: x(:,:)
 real    :: xc(3),dx,rcrit,centre(3,nb),seeds(3,nb),seeds2(3,nb),dsurf,sepmin,sepreq,r(3)
 integer :: n,nsurf,ncand,i,j,k,nconvtot

 ierr = 0
 if (nb <= 0) return
 if (nb > 127) call fatal('set_boulders','at most 127 boulders are supported')
 n = i2 - i1 + 1
 if (n < 1) call fatal('set_boulders','no particles to make boulders from')
 if (rb <= 0.) call fatal('set_boulders','boulder radius must be positive')
 !
 ! particle spacing: set_unifdis sets h = hfact*dx, and for the
 ! closepacked lattice dx is the nearest-neighbour distance
 !
 dx = xyzh(4,i1)/hfact
 if (rb < dx) call fatal('set_boulders','boulder radius is smaller than the particle spacing',val=dx)
 if (rb < 2.*xyzh(4,i1)) call warning('set_boulders','boulder radius < 2h: boulders are poorly resolved')
 !
 ! positions relative to the centroid of the body (avoids roundoff
 ! with large heliocentric positions)
 !
 allocate(x(3,n))
 xc = 0.
 do i=1,n
    xc = xc + (xyzh(1:3,i1+i-1) - xyzh(1:3,i1))
 enddo
 xc = xyzh(1:3,i1) + xc/n
 do i=1,n
    x(:,i) = xyzh(1:3,i1+i-1) - xc
 enddo
 !
 ! surface particles: fewer than the 12 nearest neighbours of an
 ! interior closepacked particle
 !
 call build_grid(fine,x,n,1.05*dx)
 allocate(isurf(n))
 nsurf = 0
 do i=1,n
    if (count_neighbours(fine,x,x(:,i),1.05*dx) - 1 < 12) then
       nsurf = nsurf + 1
       isurf(nsurf) = i
    endif
 enddo
 if (nsurf == 0) call fatal('set_boulders','could not find the surface of the body')
 !
 ! a centre is admissible if it is inside the body and no surface
 ! particle lies within rcrit, so the whole sphere exists with at least
 ! one full layer of regolith outside the surface layer
 !
 rcrit = rb + 2.*dx
 call build_grid(coarse,x,n,rcrit,isurf(1:nsurf))

 if (nb == 1) then
    !
    ! one boulder at the centroid of the body
    !
    centre(:,1) = 0.
    if (count_neighbours(fine,x,centre(:,1),dx) == 0) &
       call fatal('set_boulders','the centroid of the body lies outside it: cannot place a central boulder')
    if (.not.admissible(centre(:,1))) then
       dsurf = surface_distance(centre(:,1))
       if (id==master) print "(a,es12.4,a,es12.4)",' ERROR: boulder radius ',rb,&
                             ' does not fit at the centre; largest radius that fits = ',max(dsurf - 2.*dx,0.)
       call fatal('set_boulders','boulder does not fit at the centre of the body')
    endif
 else
    !
    ! candidate centres: particle positions where a boulder fits
    !
    allocate(icand(n))
    ncand = 0
    do i=1,n
       if (admissible(x(:,i))) then
          ncand = ncand + 1
          icand(ncand) = i
       endif
    enddo
    if (ncand == 0) then
       dsurf = 0.
       do i=1,n
          dsurf = max(dsurf,surface_distance(x(:,i)))
       enddo
       if (id==master) print "(a,es12.4,a,es12.4)",' ERROR: boulder radius ',rb,&
                             ' does not fit anywhere in the body; largest radius that fits = ',max(dsurf - 2.*dx,0.)
       call fatal('set_boulders','boulder radius too large for the body')
    endif
    !
    ! seeds spread by farthest-point sampling, then Lloyd iterations
    ! towards a centroidal Voronoi tessellation (equal-volume regions)
    !
    call farthest_point_seeds(seeds,.true.)
    centre = seeds
    call lloyd_iterations(centre)
    sepreq = 2.*rb + 2.*dx
    if (min_separation(centre) < sepreq) then
       ! push overlapping boulders apart, keeping each centre admissible
       if (id==master) print "(a)",' boulders overlap after Lloyd iterations: separating them'
       call separate_centres(centre)
    endif
    if (min_separation(centre) < sepreq) then
       ! fall back to farthest-point placement, which greedily maximises the
       ! separation, starting from the centre or from the edge of the body
       if (id==master) print "(a)",' trying farthest-point placement'
       call farthest_point_seeds(seeds2,.false.)
       if (min_separation(seeds2) > min_separation(seeds)) seeds = seeds2
       centre = seeds
       if (min_separation(centre) < sepreq) call separate_centres(centre)
    endif
    sepmin = min_separation(centre)
    if (sepmin < sepreq) then
       if (id==master) print "(2(a,es12.4))",' ERROR: closest boulder centres are ',sepmin,&
                             ' apart but need to be at least ',sepreq
       call fatal('set_boulders','boulders do not fit without overlapping: use fewer or smaller boulders')
    endif
    deallocate(icand)
 endif
 !
 ! convert the particles inside each sphere into boundary particles of
 ! that boulder, and make sure every other particle has a type
 !
 allocate(nconv(nb))
 nconv = 0
 do i=1,n
    j = i1 + i - 1
    iboulder(j) = 0
    do k=1,nb
       r = x(:,i) - centre(:,k)
       if (dot_product(r,r) <= rb*rb) then
          iboulder(j) = int(k,kind=kind(iboulder))
          nconv(k) = nconv(k) + 1
          exit
       endif
    enddo
    if (iboulder(j) > 0) then
       call set_particle_type(j,iboundary)
    else
       call set_particle_type(j,igas)
    endif
 enddo
 if (any(nconv <= 0)) call fatal('set_boulders','a boulder contains no particles')
 nconvtot = sum(nconv)
 npartoftype(igas)      = npartoftype(igas) - nconvtot
 npartoftype(iboundary) = npartoftype(iboundary) + nconvtot
 massoftype(iboundary)  = massoftype(igas)
 nboulders = nb

 if (id==master) then
    write(*,"(/,a,i4,a,es11.4,a,es11.4)") ' Rigid boulders: ',nb,' boulders of radius ',rb,', particle spacing ',dx
    write(*,"(a,3(es14.6))") ' centroid of the body = ',xc
    do k=1,nb
       write(*,"(a,i4,a,3(es11.3),a,i7,a,es11.4,a,es11.4)") ' boulder ',k,': offset from centroid = ',centre(:,k),&
             ' np = ',nconv(k),' mass = ',nconv(k)*massoftype(iboundary),&
             ' distance to surface = ',surface_distance(centre(:,k))
    enddo
    if (nb > 1) write(*,"(a,es11.4,a,es11.4,a)") ' minimum separation of boulder centres = ',min_separation(centre),&
                                                 ' (required >= ',2.*rb + 2.*dx,')'
    write(*,"(a,i8,a,/)") ' converted ',nconvtot,' particles to boundary particles'
 endif

contains

!----------------------------------------------------------------
! a centre is admissible if it is inside the body (a particle lies
! within dx of it) and no surface particle lies within rcrit
!----------------------------------------------------------------
logical function admissible(c)
 real, intent(in) :: c(3)

 admissible = (count_neighbours(fine,x,c,dx) > 0)
 if (admissible) admissible = (count_neighbours(coarse,x,c,rcrit) == 0)

end function admissible

!----------------------------------------------------------------
! distance from c to the nearest surface particle
!----------------------------------------------------------------
real function surface_distance(c)
 real, intent(in) :: c(3)
 integer :: is
 real    :: d2,dr(3)

 d2 = huge(d2)
 do is=1,nsurf
    dr = x(:,isurf(is)) - c
    d2 = min(d2,dot_product(dr,dr))
 enddo
 surface_distance = sqrt(d2)

end function surface_distance

!----------------------------------------------------------------
! choose nb candidate positions by farthest-point sampling, starting
! from the candidate nearest the centroid (from_centre) or farthest
! from it
!----------------------------------------------------------------
subroutine farthest_point_seeds(cseed,from_centre)
 real,    intent(out) :: cseed(3,nb)
 logical, intent(in)  :: from_centre
 real    :: d2min(ncand),d2,dr(3)
 integer :: ic,kk,inext

 inext = 1
 if (from_centre) then
    d2 = huge(d2)
 else
    d2 = -1.
 endif
 do ic=1,ncand
    dr = x(:,icand(ic))
    if ((from_centre .and. dot_product(dr,dr) < d2) .or. &
        (.not.from_centre .and. dot_product(dr,dr) > d2)) then
       d2 = dot_product(dr,dr)
       inext = ic
    endif
 enddo
 d2min = huge(d2)
 do kk=1,nb
    cseed(:,kk) = x(:,icand(inext))
    do ic=1,ncand
       dr = x(:,icand(ic)) - cseed(:,kk)
       d2min(ic) = min(d2min(ic),dot_product(dr,dr))
    enddo
    inext = maxloc(d2min,dim=1)
 enddo

end subroutine farthest_point_seeds

!----------------------------------------------------------------
! Lloyd iterations: move each centre to the centroid of the particles
! closest to it, snapping it to the nearest candidate position if a
! boulder would not fit there
!----------------------------------------------------------------
subroutine lloyd_iterations(c)
 real, intent(inout) :: c(3,nb)
 real    :: csum(3,nb),cnew(3),dr(3),d2,d2best,shift
 integer :: nk(nb),it,ip,kk,kbest,ic,icbest

 do it=1,maxlloyd
    csum = 0.
    nk   = 0
    do ip=1,n
       d2best = huge(d2best)
       kbest  = 1
       do kk=1,nb
          dr = x(:,ip) - c(:,kk)
          d2 = dot_product(dr,dr)
          if (d2 < d2best) then
             d2best = d2
             kbest  = kk
          endif
       enddo
       csum(:,kbest) = csum(:,kbest) + x(:,ip)
       nk(kbest) = nk(kbest) + 1
    enddo
    shift = 0.
    do kk=1,nb
       if (nk(kk) == 0) cycle
       cnew = csum(:,kk)/nk(kk)
       if (.not.admissible(cnew)) then
          d2best = huge(d2best)
          icbest = 1
          do ic=1,ncand
             dr = x(:,icand(ic)) - cnew
             d2 = dot_product(dr,dr)
             if (d2 < d2best) then
                d2best = d2
                icbest = ic
             endif
          enddo
          cnew = x(:,icand(icbest))
       endif
       shift = max(shift,sqrt(dot_product(cnew - c(:,kk),cnew - c(:,kk))))
       c(:,kk) = cnew
    enddo
    if (shift < 1.e-3*dx) exit
 enddo
 if (id==master) print "(a,i4,a)",' boulder placement: ',min(it,maxlloyd),' Lloyd iterations'

end subroutine lloyd_iterations

!----------------------------------------------------------------
! push centres that are closer than sepreq apart away from each other,
! snapping any centre that leaves the admissible region back to the
! nearest candidate position
!----------------------------------------------------------------
subroutine separate_centres(c)
 real, intent(inout) :: c(3,nb)
 real    :: dc(3,nb),dr(3),d,d2,d2best,cnew(3)
 integer :: it,kk,ll,ic,icbest
 integer, parameter :: maxsep = 500

 do it=1,maxsep
    if (min_separation(c) >= sepreq) exit
    dc = 0.
    do kk=1,nb-1
       do ll=kk+1,nb
          dr = c(:,kk) - c(:,ll)
          d  = sqrt(dot_product(dr,dr))
          if (d < sepreq) then
             if (d > 0.) then
                dr = dr/d
             else
                dr = (/1.,0.,0./)
             endif
             dc(:,kk) = dc(:,kk) + 0.5*(sepreq - d)*dr
             dc(:,ll) = dc(:,ll) - 0.5*(sepreq - d)*dr
          endif
       enddo
    enddo
    do kk=1,nb
       cnew = c(:,kk) + dc(:,kk)
       if (.not.admissible(cnew)) then
          d2best = huge(d2best)
          icbest = 1
          do ic=1,ncand
             dr = x(:,icand(ic)) - cnew
             d2 = dot_product(dr,dr)
             if (d2 < d2best) then
                d2best = d2
                icbest = ic
             endif
          enddo
          cnew = x(:,icand(icbest))
       endif
       c(:,kk) = cnew
    enddo
 enddo

end subroutine separate_centres

end subroutine set_boulders

!----------------------------------------------------------------
!+
!  minimum distance between any two of the centres
!+
!----------------------------------------------------------------
real function min_separation(c)
 real, intent(in) :: c(:,:)
 integer :: k,l
 real    :: dr(3)

 min_separation = huge(min_separation)
 do k=1,size(c,2)-1
    do l=k+1,size(c,2)
       dr = c(:,k) - c(:,l)
       min_separation = min(min_separation,sqrt(dot_product(dr,dr)))
    enddo
 enddo

end function min_separation

!----------------------------------------------------------------
!+
!  build a grid of the given cell size over the particles in x
!  (or only the subset ilist if present)
!+
!----------------------------------------------------------------
subroutine build_grid(g,x,n,cellsize,ilist)
 type(cellgrid), intent(out) :: g
 real,           intent(in)  :: x(:,:),cellsize
 integer,        intent(in)  :: n
 integer,        intent(in), optional :: ilist(:)
 integer :: m,ii,i,ic
 real    :: xmax(3)

 if (present(ilist)) then
    m = size(ilist)
 else
    m = n
 endif
 g%xmin = huge(0.)
 xmax   = -huge(0.)
 do ii=1,m
    i = ii
    if (present(ilist)) i = ilist(ii)
    g%xmin = min(g%xmin,x(:,i))
    xmax   = max(xmax,x(:,i))
 enddo
 g%size = cellsize
 g%n    = int((xmax - g%xmin)/cellsize) + 1
 allocate(g%head(product(g%n)),g%next(n))
 g%head = 0
 g%next = 0
 do ii=1,m
    i = ii
    if (present(ilist)) i = ilist(ii)
    ic = cell_index(g,cell_coords(g,x(:,i)))
    g%next(i) = g%head(ic)
    g%head(ic) = i
 enddo

end subroutine build_grid

!----------------------------------------------------------------
!+
!  integer cell coordinates of a position (not clamped to the grid)
!+
!----------------------------------------------------------------
function cell_coords(g,p) result(ic3)
 type(cellgrid), intent(in) :: g
 real,           intent(in) :: p(3)
 integer :: ic3(3)

 ic3 = int(floor((p - g%xmin)/g%size)) + 1

end function cell_coords

integer function cell_index(g,ic3)
 type(cellgrid), intent(in) :: g
 integer,        intent(in) :: ic3(3)
 integer :: jc(3)

 jc = min(max(ic3,1),g%n)
 cell_index = jc(1) + (jc(2)-1)*g%n(1) + (jc(3)-1)*g%n(1)*g%n(2)

end function cell_index

!----------------------------------------------------------------
!+
!  number of particles in the grid within radius (<= cell size) of p
!+
!----------------------------------------------------------------
integer function count_neighbours(g,x,p,radius)
 type(cellgrid), intent(in) :: g
 real,           intent(in) :: x(:,:),p(3),radius
 integer :: ic3(3),ix,iy,iz,i
 real    :: dr(3),r2

 count_neighbours = 0
 r2  = radius*radius
 ic3 = cell_coords(g,p)
 do iz=max(ic3(3)-1,1),min(ic3(3)+1,g%n(3))
    do iy=max(ic3(2)-1,1),min(ic3(2)+1,g%n(2))
       do ix=max(ic3(1)-1,1),min(ic3(1)+1,g%n(1))
          i = g%head(cell_index(g,(/ix,iy,iz/)))
          do while (i > 0)
             dr = x(:,i) - p
             if (dot_product(dr,dr) < r2) count_neighbours = count_neighbours + 1
             i = g%next(i)
          enddo
       enddo
    enddo
 enddo

end function count_neighbours

end module setboulders
