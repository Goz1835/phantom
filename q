[1mdiff --git a/src/main/config.F90 b/src/main/config.F90[m
[1mindex c63951495..f609a1846 100644[m
[1m--- a/src/main/config.F90[m
[1m+++ b/src/main/config.F90[m
[36m@@ -161,6 +161,10 @@[m [mmodule dim[m
  integer, parameter :: nalpha = 3[m
 #endif[m
 [m
[32m+[m[32m ! incompressible equation of state variables[m
[32m+[m[32m logical, parameter :: do_evolveddensity = .true.[m
[32m+[m[32m integer :: max_evolved_density = 0[m
[32m+[m
  ! default is to only store divv[m
  integer :: ndivcurlv = 1[m
  logical :: curlv = .false.[m
[36m@@ -434,6 +438,9 @@[m [msubroutine update_max_sizes(n,ntot)[m
     maxprad = maxp[m
     maxlum = maxp[m
  endif[m
[32m+[m[32m if(do_evolveddensity) then[m
[32m+[m[32m    max_evolved_density = maxp[m
[32m+[m[32m endif[m
 ! Very convoluted, but follows original logic...[m
  maxphase = maxan[m
  maxgradh = maxan[m
[1mdiff --git a/src/main/cons2prim.f90 b/src/main/cons2prim.f90[m
[1mindex 9bf1518fa..5244540df 100644[m
[1m--- a/src/main/cons2prim.f90[m
[1m+++ b/src/main/cons2prim.f90[m
[36m@@ -275,7 +275,8 @@[m [msubroutine cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,&[m
                                 Bevol,Bxyz,dustevol,dustfrac,alphaind)[m
  use part,              only:isdead_or_accreted,massoftype,igas,rhoh,igasP,iradP,iradxi,ics,imu,iX,iZ,&[m
                              iohm,ihall,nden_nimhd,eta_nimhd,iambi,get_partinfo,iphase,this_is_a_test,&[m
[31m-                             ndustsmall,itemp,ikappa,idmu,idgamma,icv,aprmassoftype,apr_level[m
[32m+[m[32m                             ndustsmall,itemp,ikappa,idmu,idgamma,icv,aprmassoftype,apr_level,&[m
[32m+[m[32m                             rhomat[m
  use part,              only:nucleation,igamma[m
  use eos,               only:equationofstate,ieos,eos_outputs_mu,done_init_eos,init_eos,gmw,X_in,Z_in,gamma[m
  use radiation_utils,   only:radiation_equation_of_state,get_opacity[m
[36m@@ -296,6 +297,7 @@[m [msubroutine cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,&[m
  real         :: rhoi,spsound,p_on_rhogas,rhogas,gasfrac,pmassi,uui[m
  real         :: Bxi,Byi,Bzi,psii,xi_limiteri,Bi,temperaturei,mui,X_i,Z_i,gammai[m
  real         :: xi,yi,zi,hi[m
[32m+[m[32m real         :: rhomati[m
  logical      :: iactivei,iamgasi,iamdusti[m
 [m
  iactivei = .true.[m
[36m@@ -311,6 +313,7 @@[m [msubroutine cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,&[m
  mui    = gmw[m
  X_i    = X_in[m
  Z_i    = Z_in[m
[32m+[m[32m rhomati = 0.[m
 [m
 !$omp parallel do default (none) &[m
 !$omp shared(xyzh,vxyzu,npart,rad,eos_vars,radprop,Bevol,Bxyz,apr_level) &[m
[36m@@ -318,10 +321,11 @@[m [msubroutine cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,&[m
 !$omp shared(alpha,alphamax,iphase,maxphase,maxp,massoftype,aprmassoftype) &[m
 !$omp shared(use_dustfrac,dustfrac,dustevol,this_is_a_test,ndustsmall,alphaind,dvdx) &[m
 !$omp shared(iopacity_type,use_var_comp,do_nucleation,update_muGamma,implicit_radiation) &[m
[32m+[m[32m!$omp shared(rhomat) &[m
 !$omp private(i,spsound,rhoi,p_on_rhogas,rhogas,gasfrac,uui) &[m
 !$omp private(Bxi,Byi,Bzi,psii,xi_limiteri,Bi,temperaturei,ierr,pmassi) &[m
 !$omp private(xi,yi,zi,hi) &[m
[31m-!$omp firstprivate(iactivei,iamtypei,iamgasi,iamdusti,mui,gammai,X_i,Z_i) &[m
[32m+[m[32m!$omp firstprivate(iactivei,iamtypei,iamgasi,iamdusti,mui,gammai,X_i,Z_i,rhomati) &[m
 !$omp reduction(+:ierrlist)[m
  do i=1,npart[m
     if (.not.isdead_or_accreted(xyzh(4,i))) then[m
[36m@@ -375,13 +379,14 @@[m [msubroutine cons2prim_everything(npart,xyzh,vxyzu,dvdx,rad,eos_vars,radprop,&[m
        endif[m
        if (use_krome) gammai = eos_vars(igamma,i)[m
        if (ieos==22) mui = eos_vars(imu,i)[m
[32m+[m[32m       if (ieos==26) rhomati = rhomat(i)[m
        if (maxvxyzu >= 4) then[m
           uui = vxyzu(4,i)[m
           if (uui < 0. .and. .not. ieos==23) then[m
              call warning('cons2prim','Internal energy < 0',i,'u',uui)[m
           endif[m
           call equationofstate(ieos,p_on_rhogas,spsound,rhogas,xi,yi,zi,temperaturei,eni=uui,&[m
[31m-                               gamma_local=gammai,mu_local=mui,Xlocal=X_i,Zlocal=Z_i)[m
[32m+[m[32m                               gamma_local=gammai,mu_local=mui,Xlocal=X_i,Zlocal=Z_i,rhomati=rhomati)[m
        else[m
           !isothermal[m
           call equationofstate(ieos,p_on_rhogas,spsound,rhogas,xi,yi,zi,temperaturei,mu_local=mui)[m
[1mdiff --git a/src/main/deriv.f90 b/src/main/deriv.f90[m
[1mindex d2881581a..9e22f2e6c 100644[m
[1m--- a/src/main/deriv.f90[m
[1m+++ b/src/main/deriv.f90[m
[36m@@ -37,7 +37,8 @@[m [mcontains[m
 subroutine derivs(icall,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&[m
                   Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&[m
                   dustevol,ddustevol,filfac,dustfrac,eos_vars,time,dt,dtnew,pxyzu,&[m
[31m-                  dens,metrics,apr_level)[m
[32m+[m[32m                  dens,metrics,apr_level,&[m
[32m+[m[32m                  drhomatdt)[m
  use dim,            only:mhd,fast_divcurlB,gr,periodic,do_radiation,driving,&[m
                           sink_radiation,use_dustgrowth,ind_timesteps,isothermal[m
  use io,             only:iprint,fatal,error[m
[36m@@ -86,6 +87,7 @@[m [msubroutine derivs(icall,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&[m
  real,            intent(out)   :: dtnew[m
  real,            intent(inout) :: pxyzu(:,:), dens(:)[m
  real,            intent(inout) :: metrics(:,:,:,:)[m
[32m+[m[32m real,            intent(out)   :: drhomatdt(:)[m
  integer(kind=1), intent(in)    :: apr_level(:)[m
  integer                     :: ierr,i[m
  real(kind=4)                :: t1,tcpu1,tlast,tcpulast[m
[36m@@ -188,7 +190,8 @@[m [msubroutine derivs(icall,npart,nactive,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&[m
  if (sinks_have_heating(nptmass,xyzmh_ptmass)) call ptmass_calc_enclosed_mass(nptmass,npart,xyzh)[m
  call force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&[m
             rad,drad,radprop,dustprop,dustgasprop,Vrel_disp,dustfrac,ddustevol,fext,fxyz_drag,&[m
[31m-            ipart_rhomax,dt,stressmax,eos_vars,dens,metrics,apr_level)[m
[32m+[m[32m            ipart_rhomax,dt,stressmax,eos_vars,dens,metrics,apr_level,&[m
[32m+[m[32m            drhomatdt)[m
  call do_timing('force',tlast,tcpulast)[m
 [m
  ![m
[36m@@ -246,7 +249,7 @@[m [msubroutine get_derivs_global(tused,dt_new,dt,icall)[m
  use part,         only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&[m
                         Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,filfac,&[m
                         dustfrac,ddustevol,eos_vars,pxyzu,dens,metrics,dustevol,gr,&[m
[31m-                        apr_level[m
[32m+[m[32m                        apr_level, drhomatdt[m
  use timing,       only:printused,getused[m
  use io,           only:id,master[m
  use cons2prim,    only:prim2consall[m
[36m@@ -274,7 +277,8 @@[m [msubroutine get_derivs_global(tused,dt_new,dt,icall)[m
  ! evaluate derivatives[m
  call derivs(icalli,npart,npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,Bevol,dBevol,&[m
              rad,drad,radprop,dustprop,ddustprop,dustevol,ddustevol,filfac,dustfrac,&[m
[31m-             eos_vars,time,dti,dtnew,pxyzu,dens,metrics,apr_level)[m
[32m+[m[32m             eos_vars,time,dti,dtnew,pxyzu,dens,metrics,apr_level,&[m
[32m+[m[32m             drhomatdt)[m
 [m
  call getused(t2)[m
  if (id==master .and. present(tused)) call printused(t1)[m
[1mdiff --git a/src/main/eos.f90 b/src/main/eos.f90[m
[1mindex cde1866a0..768346699 100644[m
[1m--- a/src/main/eos.f90[m
[1m+++ b/src/main/eos.f90[m
[36m@@ -132,7 +132,7 @@[m [mcontains[m
 !  (and position in the case of the isothermal disc)[m
 !+[m
 !----------------------------------------------------------------[m
[31m-subroutine equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi,eni,gamma_local,mu_local,Xlocal,Zlocal,radxi)[m
[32m+[m[32msubroutine equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi,eni,gamma_local,mu_local,Xlocal,Zlocal,radxi,rhomati)[m
  use io,            only:fatal,error,warning[m
  use part,          only:xyzmh_ptmass, nptmass[m
  use units,         only:unit_density,unit_pressure,unit_ergg,unit_velocity[m
[36m@@ -157,11 +157,13 @@[m [msubroutine equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi,eni,gam[m
  real,    intent(in),    optional :: eni[m
  real,    intent(inout), optional :: mu_local,gamma_local[m
  real,    intent(in),    optional :: Xlocal,Zlocal,radxi[m
[32m+[m[32m real,    intent(in),    optional :: rhomati[m
  integer :: ierr, i[m
  real    :: r1,r2[m
  real    :: mass_r, mass ! defined for generalised Farris prescription[m
  real    :: gammai,temperaturei,mui,imui,X_i,Z_i[m
  real    :: cgsrhoi,cgseni,cgspresi,presi,gam1,cgsspsoundi[m
[32m+[m[32m real    :: cgsrhomati[m
  real    :: uthermconst,kappaBar,kappaPart[m
  real    :: enthi,pondensi[m
  ![m
[36m@@ -545,9 +547,10 @@[m [msubroutine equationofstate(eos_type,ponrhoi,spsoundi,rhoi,xi,yi,zi,tempi,eni,gam[m
  case (26)[m
     !Give the function everything in cgs to ensure the values work correctly[m
     cgsrhoi = rhoi * unit_density[m
[32m+[m[32m    cgsrhomati = rhomati * unit_density[m
     !cgseni  = eni * unit_ergg[m
     call equationofstate_incomp(cgsrhoi,cgspresi,cgsspsoundi,gammai)[m
[31m-    ponrhoi  = real(cgspresi / (unit_pressure * rhoi))[m
[32m+[m[32m    ponrhoi  = real(cgspresi / (unit_pressure * rhomati))[m
     spsoundi = real(cgsspsoundi / unit_velocity)[m
 [m
  case default[m
[1mdiff --git a/src/main/eos_incomp.f90 b/src/main/eos_incomp.f90[m
[1mindex 7ed3f4b69..a21ea71e7 100644[m
[1m--- a/src/main/eos_incomp.f90[m
[1m+++ b/src/main/eos_incomp.f90[m
[36m@@ -60,7 +60,7 @@[m [msubroutine equationofstate_incomp(rho,pressure,spsound,gamma)[m
  ! Everything has been input in cgs[m
 [m
  rhos_inv = 1. / rhos ! rho s is already in code units, so just need to get the inverse[m
[31m- pressure = K * (rho * rhos_inv-1.)[m
[32m+[m[32m pressure = K * (rho * rhos_inv - 1.)[m
  gamma = 1.[m
  spsound = sqrt(K * rhos_inv)[m
 [m
[1mdiff --git a/src/main/force.F90 b/src/main/force.F90[m
[1mindex 6753827ed..73fa5fa30 100644[m
[1m--- a/src/main/force.F90[m
[1m+++ b/src/main/force.F90[m
[36m@@ -47,7 +47,8 @@[m [mmodule forces[m
 !   viscosity[m
 ![m
  use dim, only:maxfsum,maxxpartveciforce,maxp,ndivcurlB,&[m
[31m-               maxdusttypes,maxdustsmall,do_radiation,maxpsph[m
[32m+[m[32m               maxdusttypes,maxdustsmall,do_radiation,maxpsph,&[m
[32m+[m[32m               do_evolveddensity[m
  use mpiforce,    only:cellforce,stackforce[m
  use neighkdtree, only:leaf_is_active[m
  use kdtree,      only:inodeparts,inoderange[m
[36m@@ -193,7 +194,8 @@[m [mcontains[m
 subroutine force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&[m
                  rad,drad,radprop,dustprop,dustgasprop,Vrel_disp,&[m
                  dustfrac,ddustevol,fext,fxyz_drag,&[m
[31m-                 ipart_rhomax,dt,stressmax,eos_vars,dens,metrics,apr_level)[m
[32m+[m[32m                 ipart_rhomax,dt,stressmax,eos_vars,dens,metrics,apr_level, &[m
[32m+[m[32m                 drhomatdt)[m
 [m
  use dim,          only:maxvxyzu,mhd,mhd_nonideal,mpi,use_dust,use_apr,use_sinktree[m
  use io,           only:iprint,fatal,iverbose,id,master,real4,warning,error,nprocs[m
[36m@@ -261,6 +263,7 @@[m [msubroutine force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&[m
  real,            intent(out)   :: drad(:,:)[m
  real,            intent(inout) :: radprop(:,:)[m
  real,            intent(in)    :: dens(:), metrics(:,:,:,:)[m
[32m+[m[32m real,            intent(out)   :: drhomatdt(:)[m
  real, save :: xyzcache(4,maxcellcache)[m
 !$omp threadprivate(xyzcache)[m
  integer :: i,icell,nneigh[m
[36m@@ -430,6 +433,7 @@[m [msubroutine force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&[m
 !$omp shared(divcurlb) &[m
 !$omp shared(bevol) &[m
 !$omp shared(rad,radprop,drad) &[m
[32m+[m[32m!$omp shared(drhomatdt) &[m
 !$omp shared(eta_nimhd) &[m
 !$omp shared(alphaind) &[m
 !$omp shared(stressmax) &[m
[36m@@ -564,7 +568,8 @@[m [msubroutine force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&[m
                              dtradfacmax  ,dtcleanfacmax, &[m
 #endif[m
                              ndustres,dustresfacmax,dustresfacmean, &[m
[31m-                             rad,drad,radprop,dtrad)[m
[32m+[m[32m                             rad,drad,radprop,dtrad, &[m
[32m+[m[32m                             drhomatdt)[m
     endif[m
 [m
  enddo over_cells[m
[36m@@ -656,7 +661,8 @@[m [msubroutine force(icall,npart,xyzh,vxyzu,fxyzu,divcurlv,divcurlB,Bevol,dBevol,&[m
                                           dtradfacmax  ,dtcleanfacmax, &[m
 #endif[m
                                           ndustres,dustresfacmax,dustresfacmean, &[m
[31m-                                          rad,drad,radprop,dtrad)[m
[32m+[m[32m                                          rad,drad,radprop,dtrad, &[m
[32m+[m[32m                                          drhomatdt)[m
 [m
     enddo over_waiting[m
     !$omp enddo[m
[36m@@ -1552,7 +1558,7 @@[m [msubroutine compute_forces(i,iamgasi,iamdusti,xpartveci,hi,hi1,hi21,hi41,gradhi,g[m
          sxxj = 0.; sxyj = 0.; sxzj = 0.; syyj = 0.; syzj = 0.; szzj = 0.[m
        endif[m
 [m
[31m-       ifgas: if (iamgasi .and. iamgasj) then[m
[32m+[m[32m       ifgas: if (iamgasi .and. iamgasj) then[m[41m [m
           ![m
           !--artificial viscosity term[m
           ![m
[36m@@ -2443,6 +2449,7 @@[m [msubroutine start_cell(cell,iphase,xyzh,vxyzu,gradh,divcurlv,divcurlB,dvdx,Bevol,[m
           cell%apr(cell%npcell)                        = 1[m
        endif[m
 [m
[32m+[m
        alphai = alpha[m
        if (maxalpha==maxp) then[m
           alphai = alphaind(1,i)[m
[36m@@ -2679,7 +2686,8 @@[m [msubroutine finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,dv[m
                                          dtradfacmax  ,dtcleanfacmax, &[m
 #endif[m
                                          ndustres,dustresfacmax,dustresfacmean,&[m
[31m-                                         rad,drad,radprop,dtrad)[m
[32m+[m[32m                                         rad,drad,radprop,dtrad, &[m
[32m+[m[32m                                         drhomatdt)[m
 [m
  use io,             only:fatal,warning[m
  use dim,            only:mhd,mhd_nonideal,track_lum,use_dust,maxdvdx,use_dustgrowth,gr,use_krome,driving,isothermal,&[m
[36m@@ -2754,6 +2762,7 @@[m [msubroutine finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,dv[m
  real,    intent(in)    :: rad(:,:),radprop(:,:)[m
  real,    intent(out)   :: drad(:,:)[m
  real,    intent(inout) :: dtrad[m
[32m+[m[32m real,    intent(out)   :: drhomatdt(:)[m
  real    :: c_code,dtradi,radlambdai,radkappai[m
  real    :: xpartveci(maxxpartveciforce),fsum(maxfsum)[m
  real    :: rhoi,rho1i,rhogasi,hi,hi1,pmassi,tempi,gammai[m
[36m@@ -3013,6 +3022,7 @@[m [msubroutine finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,dv[m
     isgas: if (iamgasi) then[m
        divvi = -drhodti*rho1i[m
        divcurlv(1,i) = real(divvi,kind=kind(divcurlv)) ! store divv from forces[m
[32m+[m[32m       if(do_evolveddensity) drhomatdt(i) = drhodti[m
 [m
        if (.not.isothermal .or. track_lum) then[m
           if (maxdvdx == maxp .and. realviscosity) then[m
[36m@@ -3218,7 +3228,8 @@[m [msubroutine finish_cell_and_store_results(icall,cell,fxyzu,xyzh,vxyzu,poten,dt,dv[m
           !-get dust relative velocity in the kernel, sum all dimensions quadratically[m
           Vrel_disp(i) = sqrt(fsum(ivreldispxi)**2 + fsum(ivreldispyi)**2 + fsum(ivreldispzi)**2)[m
        endif[m
[31m-[m
[32m+[m[41m      [m
[32m+[m[32m       if(do_evolveddensity) drhomatdt(i) = 0.[m
        if (maxvxyzu > 4) fxyzu(4,i) = 0.[m
        ! timestep based on Courant condition for non-gas particles[m
        vsigdtc = vsigmax[m
[1mdiff --git a/src/main/initial.F90 b/src/main/initial.F90[m
[1mindex 5e507b8a6..c50fab95e 100644[m
[1m--- a/src/main/initial.F90[m
[1m+++ b/src/main/initial.F90[m
[36m@@ -101,7 +101,7 @@[m [mend subroutine initialise[m
 !+[m
 !----------------------------------------------------------------[m
 subroutine startrun(infile,logfile,evfile,dumpfile,noread)[m
[31m- use dim,              only:maxp,maxalpha,nalpha,mpi,ind_timesteps,inject_parts,gr[m
[32m+[m[32m use dim,              only:maxp,maxalpha,nalpha,mpi,ind_timesteps,inject_parts,gr,do_evolveddensity[m
  use io,               only:iprint,flush_warnings,fatal,id,master[m
  use boundary_dyn,     only:dynamic_bdy,init_dynamic_bdy[m
  use centreofmass,     only:get_centreofmass,print_particle_extent[m
[36m@@ -113,7 +113,9 @@[m [msubroutine startrun(infile,logfile,evfile,dumpfile,noread)[m
  use mpiutils,         only:reduceall_mpi[m
  use part,             only:npart,npartoftype,alphaind,ntot,update_npartoftypetot,&[m
                             iamtype,igas,nptmass,xyzmh_ptmass,vxyz_ptmass,&[m
[31m-                            xyzh,vxyzu,ibin,ibin_old,ibin_wake,ibelong[m
[32m+[m[32m                            xyzh,vxyzu,ibin,ibin_old,ibin_wake,ibelong,&[m
[32m+[m[32m                            rhomat,drhomatdt,&[m
[32m+[m[32m                            isdead_or_accreted[m
  use partinject,       only:update_injected_particles[m
  use readwrite_dumps,  only:dt_read_in[m
  use timing,           only:get_timings[m
[36m@@ -121,6 +123,7 @@[m [msubroutine startrun(infile,logfile,evfile,dumpfile,noread)[m
  use timestep_ind,     only:ibinnow,init_ibin,istepfrac,nbinmax[m
  use writeheader,      only:write_header[m
  use metric,           only:update_metric[m
[32m+[m[32m use granular_variables,         only:rhos[m
  character(len=*), intent(inout) :: infile[m
  character(len=*), intent(out)   :: logfile,evfile,dumpfile[m
  logical,          intent(in), optional :: noread[m
[36m@@ -205,6 +208,19 @@[m [msubroutine startrun(infile,logfile,evfile,dumpfile,noread)[m
  call get_derivs_initial(time,dumpfile,ntot,dtnew_first,ierr)[m
  if (ierr /= 0) call fatal('initial','error in derivatives setup')[m
 [m
[32m+[m[32m ! set initial value of the evolved density to be the density in granular[m
[32m+[m[32m if (do_evolveddensity) then[m
[32m+[m[32m    do i=1,npart[m
[32m+[m[32m       !might not need the isdead check[m
[32m+[m[32m       if (isdead_or_accreted(xyzh(4,i))) then[m
[32m+[m[32m          rhomat(i) = 0.[m
[32m+[m[32m       else[m
[32m+[m[32m          rhomat(i) = rhos / unit_density[m
[32m+[m[32m       endif[m
[32m+[m[32m       drhomatdt(i) = 0.[m
[32m+[m[32m    enddo[m
[32m+[m[32m endif[m
[32m+[m
  ! set initial timestep[m
  if (.not.ind_timesteps) then[m
     dt = min(dtnew_first,dtinject)[m
[1mdiff --git a/src/main/part.F90 b/src/main/part.F90[m
[1mindex b7b333cc2..c40e96649 100644[m
[1m--- a/src/main/part.F90[m
[1m+++ b/src/main/part.F90[m
[36m@@ -34,7 +34,8 @@[m [mmodule part[m
                maxTdust,store_dust_temperature,use_krome,maxp_krome, &[m
                do_radiation,gr,maxgr,maxgran,n_nden_phantom,do_nucleation,&[m
                inucleation,itau_alloc,itauL_alloc,inject_parts,use_apr,apr_maxlevel,&[m
[31m-               maxp_apr,maxptmassgr,use_sinktree,nvel_ptmass[m
[32m+[m[32m               maxp_apr,maxptmassgr,use_sinktree,nvel_ptmass,&[m
[32m+[m[32m               do_evolveddensity,max_evolved_density[m
  use dtypekdtree, only:kdnode[m
 #ifdef KROME[m
  use krome_user, only:krome_nmols[m
[36m@@ -126,6 +127,13 @@[m [mmodule part[m
  integer, parameter :: iCO       = 5 ! CO abundance[m
  real, allocatable :: abundance(:,:)[m
 ![m
[32m+[m[32m!--evolved density[m
[32m+[m[32m![m
[32m+[m[32m real, allocatable :: rhomat(:)    ! evolved quantity[m
[32m+[m[32m real, allocatable :: drhomatdt(:)   ! its time derivative[m
[32m+[m[32m real, allocatable :: rhomatpred(:)    ! predicted value (for the predictor-corrector step)[m
[32m+[m[32m character(len=*), parameter :: rhomat_label = 'rhomat'[m
[32m+[m[32m![m
 !--KROME chemistry[m
 ![m
 #ifdef KROME[m
[36m@@ -556,6 +564,9 @@[m [msubroutine allocate_part[m
  call allocate_array("shortsinktree", shortsinktree, maxptmass, maxptmass)[m
  call allocate_array("gtgrad", gtgrad, 3, maxptmass)[m
  call allocate_array('noverlap', noverlap, maxp)[m
[32m+[m[32m call allocate_array('rhomat', rhomat, max_evolved_density)[m
[32m+[m[32m call allocate_array('drhomatdt', drhomatdt, max_evolved_density)[m
[32m+[m[32m call allocate_array('rhomatpred', rhomatpred, max_evolved_density)[m
 [m
 end subroutine allocate_part[m
 [m
[36m@@ -648,6 +659,7 @@[m [msubroutine deallocate_part[m
  if (allocated(shortsinktree))deallocate(shortsinktree)[m
  if (allocated(gtgrad))       deallocate(gtgrad)[m
  if (allocated(noverlap))     deallocate(noverlap)[m
[32m+[m[32m if (allocated(rhomat))     deallocate(rhomat, drhomatdt, rhomatpred)[m
 [m
 end subroutine deallocate_part[m
 [m
[36m@@ -729,6 +741,12 @@[m [msubroutine init_part[m
     twas(:)       = 0.[m
  endif[m
  if (use_sinktree) iphase(maxpsph+1:maxp) = isink[m
[32m+[m[32m !-- Evolved density for μI rheological model (incompressible eos)[m
[32m+[m[32m if(do_evolveddensity) then[m
[32m+[m[32m    rhomat(:) = 0.[m
[32m+[m[32m    drhomatdt(:) = 0.[m
[32m+[m[32m    rhomatpred(:) = 0.[m
[32m+[m[32m endif[m
 [m
  ideadhead = 0[m
 ![m
[36m@@ -1305,6 +1323,9 @@[m [msubroutine copy_particle(src,dst,new_part)[m
     iorig(dst) = iorig(src) ! we are moving the particle within the list; maintain ID[m
  endif[m
  if (inject_parts) iseed_sink(dst) = iseed_sink(src)[m
[32m+[m[32m if(do_evolveddensity) then[m
[32m+[m[32m   rhomat(dst) = rhomat(src)[m
[32m+[m[32m endif[m
 [m
 end subroutine copy_particle[m
 [m
[36m@@ -1418,6 +1439,11 @@[m [msubroutine copy_particle_all(src,dst,new_part)[m
     iorig(dst) = iorig(src) ! we are moving the particle within the list; maintain ID[m
  endif[m
  if (inject_parts) iseed_sink(dst) = iseed_sink(src)[m
[32m+[m[32m if (do_evolveddensity) then[m
[32m+[m[32m    rhomat(dst) = rhomat(src)[m
[32m+[m[32m    drhomatdt(dst) = drhomatdt(src)[m
[32m+[m[32m    rhomatpred(dst) = rhomatpred(src)[m
[32m+[m[32m endif[m
 [m
 end subroutine copy_particle_all[m
 [m
[36m@@ -1517,6 +1543,12 @@[m [msubroutine combine_two_particles(keep,discard)[m
  if (itau_alloc == 1) tau(keep) = 0.5*(tau(keep) + tau(discard))[m
  if (itauL_alloc == 1) tau_lucy(keep) = 0.5*(tau_lucy(keep) + tau_lucy(discard))[m
 [m
[32m+[m[32m if(do_evolveddensity) then[m
[32m+[m[32m    rhomat(keep) = 0.5 * (rhomat(keep) + rhomat(discard))[m
[32m+[m[32m    drhomatdt(keep) = 0.5 * (drhomatdt(keep) + drhomatdt(discard))[m
[32m+[m[32m    rhomatpred(keep) = 0.5 * (rhomatpred(keep) + rhomatpred(discard))[m
[32m+[m[32m endif[m
[32m+[m
  if (use_krome) then[m
     T_gas_cool(keep)       = 0.5*(T_gas_cool(keep) + T_gas_cool(discard))[m
  endif[m
[36m@@ -1736,6 +1768,11 @@[m [msubroutine fill_sendbuf(i,xtemp,nbuf)[m
        call fill_buffer(xtemp,dt_in(i),nbuf)[m
        call fill_buffer(xtemp,twas(i),nbuf)[m
     endif[m
[32m+[m[32m    if(do_evolveddensity) then[m
[32m+[m[32m       call fill_buffer(xtemp,rhomat(i), nbuf)[m
[32m+[m[32m       call fill_buffer(xtemp,drhomatdt(i), nbuf)[m
[32m+[m[32m       call fill_buffer(xtemp,rhomatpred(i), nbuf)[m
[32m+[m[32m    endif[m
     call fill_buffer(xtemp,iorig(i),nbuf)[m
     if (inject_parts) call fill_buffer(xtemp,iseed_sink(i),nbuf)[m
     if (use_apr) call fill_buffer(xtemp,apr_level(i),nbuf)[m
[36m@@ -1823,6 +1860,11 @@[m [msubroutine unfill_buffer(ipart,xbuf)[m
     dt_in(ipart)        = real(unfill_buf(xbuf,j),kind=kind(dt_in))[m
     twas(ipart)         = unfill_buf(xbuf,j)[m
  endif[m
[32m+[m[32m if(do_evolveddensity) then[m
[32m+[m[32m    rhomat(ipart) = unfill_buf(xbuf,j)[m
[32m+[m[32m    drhomatdt(ipart) = unfill_buf(xbuf,j)[m
[32m+[m[32m    rhomatpred(ipart) = unfill_buf(xbuf,j)[m
[32m+[m[32m endif[m
  iorig(ipart)           = nint(unfill_buf(xbuf,j),kind=8)[m
  if (inject_parts) iseed_sink(ipart) = nint(unfill_buf(xbuf,j),kind=8)[m
  if (use_apr) apr_level(ipart) = nint(unfill_buf(xbuf,j),kind=kind(apr_level))[m
[1mdiff --git a/src/main/readwrite_dumps.f90 b/src/main/readwrite_dumps.f90[m
[1mindex 602858c8d..b1a142b55 100644[m
[1m--- a/src/main/readwrite_dumps.f90[m
[1m+++ b/src/main/readwrite_dumps.f90[m
[36m@@ -29,6 +29,7 @@[m [mmodule readwrite_dumps[m
 [m
  logical, public :: opened_full_dump[m
  logical, public :: dt_read_in[m
[32m+[m[32m logical, public :: got_evolve_dens = .false.[m[41m [m
 [m
  integer, parameter :: maxphead = 256         ! max items in header[m
  integer, parameter :: is_small_dump = 1978[m
[36m@@ -46,7 +47,7 @@[m [msubroutine write_fulldump(t,dumpfile,ntotal,iorder,sphNG)[m
  use dim,   only:maxp,maxvxyzu,maxalpha,ndivcurlv,ndivcurlB,maxgrav,gravity,use_dust,&[m
                    track_lum,use_dustgrowth,store_dust_temperature,gr,do_nucleation,&[m
                    ind_timesteps,mhd_nonideal,use_krome,h2chemistry,update_muGamma,mpi,use_apr,&[m
[31m-                   inject_parts[m
[32m+[m[32m                   inject_parts,do_evolveddensity[m
  use eos,   only:ieos,eos_is_non_ideal,eos_outputs_mu,eos_outputs_gasP,eos_outputs_temp,eos_outputs_gamma[m
  use io,    only:idump,iprint,real4,id,master,error,warning,nprocs[m
  use part,  only:xyzh,xyzh_label,vxyzu,vxyzu_label,Bevol,Bevol_label,Bxyz,Bxyz_label,npart,maxtypes, &[m
[36m@@ -60,7 +61,8 @@[m [msubroutine write_fulldump(t,dumpfile,ntotal,iorder,sphNG)[m
  !dvdx,dvdx_label,&[m
                    rad,rad_label,radprop,radprop_label,do_radiation,maxirad,maxradprop,itemp,igasP,igamma,&[m
                    iorig,iseed_sink,iX,iZ,imu,nucleation,nucleation_label,n_nucleation,tau,itau_alloc,tau_lucy,itauL_alloc,&[m
[31m-                   luminosity,eta_nimhd,eta_nimhd_label,apr_level[m
[32m+[m[32m                   luminosity,eta_nimhd,eta_nimhd_label,apr_level,&[m
[32m+[m[32m                   rhomat,rhomat_label[m
  use part,  only:metrics,metricderivs,tmunus[m
  use options,    only:use_dustfrac,use_porosity,use_var_comp,icooling[m
  use dump_utils, only:tag,open_dumpfile_w,allocate_header,&[m
[36m@@ -278,6 +280,7 @@[m [msubroutine write_fulldump(t,dumpfile,ntotal,iorder,sphNG)[m
        if (do_nucleation) call write_array(1,nucleation,nucleation_label,n_nucleation,npart,k,ipass,idump,nums,nerr)[m
        if (itau_alloc == 1) call write_array(1,tau,'tau',npart,k,ipass,idump,nums,nerr)[m
        if (itauL_alloc == 1) call write_array(1,tau_lucy,'tau_lucy',npart,k,ipass,idump,nums,nerr)[m
[32m+[m[32m       if (do_evolveddensity) call write_array(1,rhomat,rhomat_label,npart,k,ipass,idump,nums,nerr)[m
        if (store_dust_temperature) call write_array(1,dust_temp,'Tdust',npart,k,ipass,idump,nums,nerr)[m
        if (do_radiation) then[m
           call write_array(1,rad,rad_label,maxirad,npart,k,ipass,idump,nums,nerr)[m
[36m@@ -510,6 +513,7 @@[m [msubroutine read_dump(dumpfile,tfile,hfactfile,idisk1,iprint,id,nprocs,ierr,heade[m
  if (id==master .and. iverbose >= 0) write(iprint,"(/,1x,a,i3)") '>>> reading setup from file: '//trim(dumpfile)//' on unit ',idisk1[m
  opened_full_dump = .true.[m
  dt_read_in = .false.[m
[32m+[m[32m got_evolve_dens =.false.[m
  ![m
  ! open dump file[m
  ![m
[36m@@ -960,7 +964,7 @@[m [msubroutine read_phantom_arrays(i1,i2,noffset,narraylengths,nums,npartread,nparto[m
  use dump_utils, only:read_array,match_tag[m
  use dim,        only:use_dust,h2chemistry,maxalpha,maxp,gravity,maxgrav,maxvxyzu,do_nucleation, &[m
                         use_dustgrowth,maxdusttypes,maxphase,gr,store_dust_temperature,&[m
[31m-                        ind_timesteps,use_krome,use_apr,mhd,inject_parts[m
[32m+[m[32m                        ind_timesteps,use_krome,use_apr,mhd,inject_parts,do_evolveddensity[m
  use part,       only:xyzh,xyzh_label,vxyzu,vxyzu_label,dustfrac,dustfrac_label,abundance,abundance_label, &[m
                         alphaind,poten,xyzmh_ptmass,xyzmh_ptmass_label,vxyz_ptmass,vxyz_ptmass_label, &[m
                         Bevol,Bxyz,Bxyz_label,nabundances,iphase,idust, &[m
[36m@@ -968,7 +972,8 @@[m [msubroutine read_phantom_arrays(i1,i2,noffset,narraylengths,nums,npartread,nparto[m
                         VrelVf,VrelVf_label,dustgasprop,dustgasprop_label,filfac,filfac_label,pxyzu,pxyzu_label,dust_temp, &[m
                         rad,rad_label,radprop,radprop_label,do_radiation,maxirad,maxradprop,ifluxx,ifluxy,ifluxz, &[m
                         nucleation,nucleation_label,n_nucleation,ikappa,tau,itau_alloc,tau_lucy,itauL_alloc,&[m
[31m-                        ithick,ilambda,iorig,iseed_sink,dt_in,krome_nmols,T_gas_cool,apr_level[m
[32m+[m[32m                        ithick,ilambda,iorig,iseed_sink,dt_in,krome_nmols,T_gas_cool,apr_level,&[m
[32m+[m[32m                        rhomat,rhomat_label[m
  use eos_stamatellos, only:ttherm_store,ueqi_store,tau_store,du_store[m
  use sphNGutils, only:mass_sphng,got_mass,set_gas_particle_mass[m
  use options,    only:use_porosity[m
[36m@@ -1089,6 +1094,9 @@[m [msubroutine read_phantom_arrays(i1,i2,noffset,narraylengths,nums,npartread,nparto[m
              if (itauL_alloc == 1) then[m
                 call read_array(tau_lucy,'tau_lucy',got_ray_tracer,ik,i1,i2,noffset,idisk1,tag,match,ierr)[m
              endif[m
[32m+[m[32m             if (do_evolveddensity) then[m
[32m+[m[32m                call read_array(rhomat,rhomat_label,got_evolve_dens,ik,i1,i2,noffset,idisk1,tag,match,ierr)[m
[32m+[m[32m             endif[m
              if (store_dust_temperature) then[m
                 call read_array(dust_temp,'Tdust',got_Tdust,ik,i1,i2,noffset,idisk1,tag,match,ierr)[m
              endif[m
[1mdiff --git a/src/main/step_leapfrog.f90 b/src/main/step_leapfrog.f90[m
[1mindex 268d2f9e9..e6eb13772 100644[m
[1m--- a/src/main/step_leapfrog.f90[m
[1m+++ b/src/main/step_leapfrog.f90[m
[36m@@ -90,7 +90,8 @@[m [mend subroutine init_step[m
 !------------------------------------------------------------[m
 subroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
  use dim,            only:maxp,maxvxyzu,maxptmass,maxalpha,nalpha,h2chemistry,&[m
[31m-                          use_dustgrowth,use_krome,gr,do_radiation,use_apr,use_sinktree[m
[32m+[m[32m                          use_dustgrowth,use_krome,gr,do_radiation,use_apr,use_sinktree,&[m
[32m+[m[32m                          do_evolveddensity[m
  use io,             only:iprint,fatal,iverbose,id,master,warning[m
  use options,        only:iexternalforce,use_dustfrac,implicit_radiation,&[m
                           use_porosity,icooling[m
[36m@@ -102,7 +103,8 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
                           dustfrac,dustevol,ddustevol,eos_vars,alphaind,nptmass,&[m
                           dustprop,ddustprop,dustproppred,pxyzu,dens,metrics,ics,&[m
                           filfac,filfacpred,mprev,filfacprev,aprmassoftype,&[m
[31m-                          fxyz_ptmass_tree[m
[32m+[m[32m                          fxyz_ptmass_tree,&[m
[32m+[m[32m                          rhomat,drhomatdt,rhomatpred[m
  use part,           only:nptmass,xyzmh_ptmass,vxyz_ptmass,fxyz_ptmass, &[m
                           dsdt_ptmass,fsink_old,ibin_wake,dptmass, &[m
                           pxyzu_ptmass,metrics_ptmass[m
[36m@@ -177,6 +179,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
  !$omp shared(dustprop,ddustprop,dustproppred,ufloor,icooling,Tfloor) &[m
  !$omp shared(mprev,filfacprev,filfac,use_porosity) &[m
  !$omp shared(ibin,ibin_old,twas,timei) &[m
[32m+[m[32m !$omp shared(rhomat, drhomatdt) &[m
  !$omp firstprivate(itype) &[m
  !$omp private(i,hdti) &[m
  !$omp reduction(+:nvfloorp)[m
[36m@@ -225,6 +228,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
        if (itype==igas) then[m
           if (mhd)          Bevol(:,i) = Bevol(:,i) + hdti*dBevol(:,i)[m
           if (do_radiation) rad(:,i)   = rad(:,i) + hdti*drad(:,i)[m
[32m+[m[32m          if (do_evolveddensity) rhomat(i) = rhomat(i) + hdti*drhomatdt(i)[m[41m [m
           if (use_dustfrac) then[m
              dustevol(:,i) = abs(dustevol(:,i) + hdti*ddustevol(:,i))[m
              if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)[m
[36m@@ -300,6 +304,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
 !$omp shared(eos_vars,ufloor,icooling,Tfloor) &[m
 !$omp shared(twas,timei) &[m
 !$omp shared(rad,drad,radpred)&[m
[32m+[m[32m!$omp shared(rhomat, drhomatdt, rhomatpred) &[m
 !$omp private(hi,rhoi,tdecay1,source,ddenom,hdti) &[m
 !$omp private(i,spsoundi,alphaloci) &[m
 !$omp firstprivate(pmassi,itype,alpha) &[m
[36m@@ -323,6 +328,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
              if (use_dustgrowth) dustproppred(:,i) = dustprop(:,i)[m
              if (use_dustfrac)   dustpred(:,i) = dustevol(:,i)[m
              if (do_radiation)   radpred(:,i) = rad(:,i)[m
[32m+[m[32m             if (do_evolveddensity) rhomatpred(i) = rhomat(i)[m
              cycle predict_sph[m
           endif[m
        endif[m
[36m@@ -371,6 +377,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
              if (use_dustgrowth) dustproppred(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)[m
           endif[m
           if (do_radiation) radpred(:,i) = rad(:,i) + hdti*drad(:,i)[m
[32m+[m[32m          if (do_evolveddensity) rhomatpred(i) = rhomat(i) + hdti*drhomatdt(i)[m
        endif[m
        ![m
        ! viscosity switch ONLY (conductivity and resistivity do not use MM97-style switches)[m
[36m@@ -417,7 +424,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
     call derivs(1,npart,nactive,xyzh,vpred,fxyzu,fext,divcurlv,&[m
                 divcurlB,Bpred,dBevol,radpred,drad,radprop,dustproppred,ddustprop,&[m
                 dustpred,ddustevol,filfacpred,dustfrac,eos_vars,timei,dtsph,dtnew,&[m
[31m-                ppred,dens,metrics,apr_level)[m
[32m+[m[32m                ppred,dens,metrics,apr_level,drhomatdt)[m
     if (do_radiation .and. implicit_radiation) then[m
        rad = radpred[m
        vxyzu(4,1:npart) = vpred(4,1:npart)[m
[36m@@ -472,6 +479,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
 !$omp shared(ibin_dts,nbinmax) &[m
 !$omp private(dti,hdti) &[m
 !$omp shared(rad,radpred,drad)&[m
[32m+[m[32m!$omp shared(rhomat,drhomatdt) &[m
 !$omp private(i,vxi,vyi,vzi) &[m
 !$omp private(pxi,pyi,pzi,p2i) &[m
 !$omp private(erri,v2i,eni) &[m
[36m@@ -507,6 +515,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
                 if (itype==igas) then[m
                    if (mhd)          Bevol(:,i) = Bevol(:,i) + dti*dBevol(:,i)[m
                    if (do_radiation) rad(:,i)   = rad(:,i)   + dti*drad(:,i)[m
[32m+[m[32m                   if (do_evolveddensity) rhomat(i) = rhomat(i) + dti*drhomatdt(i)[m
                    if (use_dustfrac) then[m
                       dustevol(:,i) = dustevol(:,i) + dti*ddustevol(:,i)[m
                       if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) + dti*ddustprop(:,i)[m
[36m@@ -542,6 +551,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
              if (itype==igas) then[m
                 if (mhd)          Bevol(:,i) = Bevol(:,i) + hdti*dBevol(:,i)[m
                 if (do_radiation) rad(:,i)   = rad(:,i)   + hdti*drad(:,i)[m
[32m+[m[32m                if (do_evolveddensity) rhomat(i) = rhomat(i) + hdti*drhomatdt(i)[m
                 if (use_dustfrac) then[m
                    dustevol(:,i) = dustevol(:,i) + hdti*ddustevol(:,i)[m
                    if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) + hdti*ddustprop(:,i)[m
[36m@@ -613,6 +623,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
                 ![m
                 if (mhd)          Bevol(:,i) = Bevol(:,i)  + hdtsph*dBevol(:,i)[m
                 if (do_radiation) rad(:,i)   = rad(:,i) + hdtsph*drad(:,i)[m
[32m+[m[32m                if (do_evolveddensity) rhomat(i) = rhomat(i) + hdtsph*drhomatdt(i)[m
                 if (use_dustfrac) then[m
                    dustevol(:,i) = dustevol(:,i) + hdtsph*ddustevol(:,i)[m
                    if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) + hdtsph*ddustprop(:,i)[m
[36m@@ -645,6 +656,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
 !$omp shared(dustprop,ddustprop,dustproppred,use_dustfrac,dustevol,dustpred,ddustevol) &[m
 !$omp shared(filfac,filfacpred,use_porosity) &[m
 !$omp shared(rad,drad,radpred,icooling,Tfloor,xyzh) &[m
[32m+[m[32m!$omp shared(rhomat, drhomatdt, rhomatpred) &[m
 !$omp firstprivate(itype) &[m
 !$omp schedule(static)[m
        until_converged: do i=1,npart[m
[36m@@ -675,6 +687,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
              if (mhd)          Bpred(:,i)  = Bevol(:,i)[m
              if (use_dustfrac) dustpred(:,i) = dustevol(:,i)[m
              if (do_radiation) radpred(:,i) = rad(:,i)[m
[32m+[m[32m             if (do_evolveddensity) rhomatpred(i) = rhomat(i)[m
              ![m
              ! shift v back to the half step[m
              ![m
[36m@@ -696,6 +709,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
                    if (use_dustgrowth) dustprop(:,i) = dustprop(:,i) - hdtsph*ddustprop(:,i)[m
                 endif[m
                 if (do_radiation) rad(:,i) = rad(:,i) - hdtsph*drad(:,i)[m
[32m+[m[32m                if (do_evolveddensity) rhomat(i) = rhomat(i) - hdtsph*drhomatdt(i)[m
              endif[m
           endif[m
        enddo until_converged[m
[36m@@ -713,7 +727,7 @@[m [msubroutine step(npart,nactive,t,dtsph,dtextforce,dtnew)[m
        if (gr) vpred = vxyzu ! Need primitive utherm as a guess in cons2prim[m
        call derivs(2,npart,nactive,xyzh,vpred,fxyzu,fext,divcurlv,divcurlB, &[m
                      Bpred,dBevol,radpred,drad,radprop,dustproppred,ddustprop,dustpred,ddustevol,filfacpred,&[m
[31m-                     dustfrac,eos_vars,timei,dtsph,dtnew,ppred,dens,metrics,apr_level)[m
[32m+[m[32m                     dustfrac,eos_vars,timei,dtsph,dtnew,ppred,dens,metrics,apr_level,drhomatdt)[m
        if (gr) vxyzu = vpred ! May need primitive variables elsewhere?[m
        if (do_radiation .and. implicit_radiation) then[m
           rad = radpred[m
