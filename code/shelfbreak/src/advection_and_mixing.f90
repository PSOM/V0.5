subroutine advection_and_mixing(m,n,dtimel,step) 
#include "cppdefs.h"
!     ---------------------------------------------                     
USE header
!     convecn(m,n,dtime) means use level m and write s,T to level n     
!     computes Cx,Cy,Cz at the cell centers.                            
!     uflx,vflx,wflx,sflx,Tflx are the fluxes defined at                
!     cell faces (using QUICK).                                         
!                                                                       

INTEGER :: n,m,step

REAL(kind=rc_kind) :: dtimel

REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: var      ! var=(s,T,u,v,w,Tr(it,:,:,:)                    
REAL(kind=rc_kind), dimension(    0:NI+1,0:NJ+1, 0:NK+1) :: uvarx    ! uvarx  is the source term, divergence of the advective fluxes 
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif   ! vardif is the source term from diabatic processes 
REAL(kind=rc_kind), dimension(    1:NI  ,1:NJ  , 1:NK  ) :: vardif2   ! vardif is the source term from diabatic processes 

INTEGER :: selvar
INTEGER :: iv_compute_kz
INTEGER :: av_comp(5+ntr)

iv_compute_kz=1

! By default, all variables will be subject to the following loop
av_comp=1

! If the simulation is a "rhoonly" simulation (rho is stored in s), T=0 and does not need to be updated.
#ifdef rhoonly
av_comp(1)=0
#endif

do selvar=1,5+ntr
  if(av_comp(selvar)==1) then

    if(selvar==1) var=T(:,:,:,m)
    if(selvar==2) var=s(:,:,:,m)
    if(selvar==3) var=u(:,:,:,m)
    if(selvar==4) var=v(:,:,:,m)
    if(selvar==5) var=w(:,:,:,m)
    do it=1,ntr
      if(selvar==5+it) var=Tr(it,:,:,:,m)
    enddo


    ! computation of the advective fluxes, using QUICK scheme
    CALL advect(var,uvarx)
                         
    ! computation of the horizontal diabatic fluxes
    if(selvar<5.5) then

      vardif=0.; vardif2=0.

!     if(ABS(user1-1.)<1e-5) then
       call mixing_horizontal(var,vardif);!PRINT*,"VISCOUS HORI" 
!       PRINT*, "MW ",MAXVAL(ABS(W))
!     endif

!     if(ABS(user1-2.)<1e-5) then
!       call mixing_isopycnal_biharmonic(var,vardif,0.1);PRINT*,"BIHARMONIC REDI"
!       vardif=vardif+vardif2; 
!     endif

!     if(ABS(user1-3.)<1e-5) then
!       call mixing_isopycnal_biharmonic(var,vardif,0.01);PRINT*,"BIHARMONIC REDI"
!       vardif=vardif+vardif2; 
!     endif

!     if(ABS(user1-4.)<1e-5) then
!       call mixing_isopycnal_biharmonic(var,vardif,0.001);PRINT*,"BIHARMONIC REDI"
!       vardif=vardif+vardif2; 
!     endif

!     if(ABS(user1-5.)<1e-5) then
!       call mixing_isopycnal(var,vardif,1.);PRINT*,"VISCOUS REDI";
!     endif
!
!     if(ABS(user1-6.)<1e-5) then
!       call mixing_isopycnal(var,vardif,0.1);PRINT*,"VISCOUS REDI";
!     endif

!     if(ABS(user1-7.)<1e-5) then
!       call mixing_horizontal(var,vardif);PRINT*,"VISCOUS HORI" 
!     endif

!     if(ABS(user1-8.)<1e-5) then
!       vardif=0.;PRINT*,"NO VISCOUS"
!     endif

!     if(ABS(user1-9.)<1e-5) then
!       vardif=0.;PRINT*,"NO VISCOUS"
!     endif


      uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
    endif
                       
    ! computation of the vertical diabatic fluxes
    if(selvar<4.5) then
      vardif=0.;
!     if(user1<7.5) then
     call mixing_vertical(var,vardif,m,step,iv_compute_kz); 
!     endif

!     if(user1>7.5) then
!       vardif=0.;
!     endif

      uvarx(1:NI,1:NJ,1:NK)=uvarx(1:NI,1:NJ,1:NK)-vardif(1:NI,1:NJ,1:NK)
      iv_compute_kz=0;
    endif

    ! final summation  
    if(selvar==1) T(1:NI,1:NJ,1:NK,n)=T(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    if(selvar==2) s(1:NI,1:NJ,1:NK,n)=s(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    if(selvar==3) cx(1:NI,1:NJ,1:NK) =u(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    if(selvar==4) cy(1:NI,1:NJ,1:NK) =v(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    if(selvar==5) cz(1:NI,1:NJ,1:NK) =w(1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    do it=1,ntr
      if(selvar==5+it) Tr(it,1:NI,1:NJ,1:NK,n) =Tr(it,1:NI,1:NJ,1:NK,0)-dtimel*Jacinv(1:NI,1:NJ,1:NK)*uvarx(1:NI,1:NJ,1:NK)
    enddo

  endif
enddo !selvar

#ifdef allow_particle
  PRINT*,"ALLOW PARTICLE"
#endif


return 
                                                                        
END

                                           
