
  subroutine vec2xy(lat_points,lon_points,numpatch,mempatch,&
		    ixy_patch,jxy_patch,wtxy_patch,itypwat,nfcon,nforc,nfldv,&
                    fcon,forcxy,fldv,fldxy_r,msk_patch)

! ----------------------------------------------------------------------
! perfrom grid-average from subgrid 1d vector
! subgrid to grid average mapping: average a subgrid input vector [fldv] 
! of length to a 2-d [lon_points] x [lat_points] output array [fldxy_r]
!
! Created by Yongjiu Dai
!--------------!-------------------------------------------------------
! 01: taux     ! wind stress: E-W [kg/m/s2]
! 02: tauy     ! wind stress: N-S [kg/m/s2]
! 03: fsena    ! sensible heat from canopy height to atmosphere [W/m2]
! 04: lfevpa   ! latent heat flux from canopy height to atmosphere [W/m2]
! 05: fevpa    ! evapotranspiration from canopy to atmosphere [mm/s]
! 06: fsenl    ! sensible heat from leaves [W/m2]
! 07: fevpl    ! evaporation+transpiration from leaves [mm/s]
! 08: etr      ! transpiration rate [mm/s]
! 09: fseng    ! sensible heat flux from ground [W/m2]
! 10: fevpg    ! evaporation heat flux from ground [mm/s]
! 11: fgrnd    ! ground heat flux [W/m2]
! 12: sabvsun  ! solar absorbed by sunlit canopy [W/m2]
! 13: sabvsha  ! solar absorbed by shaded [W/m2]
! 14: sabg     ! solar absorbed by ground  [W/m2]
! 15: olrg     ! outgoing long-wave radiation from ground+canopy [W/m2]
! 16: rnet     ! net radiation [W/m2]
! 17: xerr     ! the error of water banace [mm/s]
! 18: zerr     ! the error of energy balance [W/m2]
! 19: rsur     ! surface runoff [mm/s]
! 20: rnof     ! total runoff [mm/s]
! 21: assim    ! canopy assimilation rate [mol m-2 s-1]
! 22: respc    ! respiration (plant+soil) [mol m-2 s-1]
!--------------!-------------------------------------------------------
! 23:32: tss   ! soil temperature [K]
! 33:42: wliq  ! liquid water in soil layers [kg/m2]
! 43:52: wice  ! ice lens in soil layers [kg/m2]
!--------------!-------------------------------------------------------
! 53: tg       ! ground surface temperature [K]
! 54: tlsun    ! sunlit leaf temperature [K]
! 55: tlsha    ! shaded leaf temperature [K]
! 56: ldew     ! depth of water on foliage [mm]
! 57: scv      ! snow cover, water equivalent [mm]
! 58: snowdp   ! snow depth [meter]
! 59: fsno     ! fraction of snow cover on ground
! 60: sigf     ! fraction of veg cover, excluding snow-covered veg [-]
! 61: green    ! leaf greenness
! 62: lai      ! leaf area index
! 63: sai      ! stem area index
! 64: alb(1,1) ! averaged albedo [visible, direct]
! 65: alb(1,2) ! averaged albedo [visible, diffuse]
! 66: alb(2,1) ! averaged albedo [near-infrared, direct]
! 67: alb(2,2) ! averaged albedo [near-infrared,diffuse]
! 68: emis     ! averaged bulk surface emissivity
! 69: z0ma     ! effective roughness [m]
!--------------!-------------------------------------------------------
! 70: trad     ! radiative temperature of surface [K]
! 71: ustar    ! u* in similarity theory [m/s]
! 72: tstar    ! t* in similarity theory [kg/kg]
! 73: qstar    ! q* in similarity theory [kg/kg]
! 74: zol      ! dimensionless height (z/L) used in Monin-Obukhov theory
! 75: rib      ! bulk Richardson number in surface layer
! 76: fm       ! integral of profile function for momentum
! 77: fh       ! integral of profile function for heat
! 78: fq       ! integral of profile function for moisture
!--------------!-------------------------------------------------------
! 79: tref     ! 2 m height air temperature [kelvin]
! 80: qref     ! 2 m height air specific humidity [kg/kg]
! 81: u10m     ! 10m u-velocity [m/s]
! 82: v10m     ! 10m v-velocity [m/s]
! 83: f10m     ! integral of profile function for momentum at 10m [-]
!--------------!-------------------------------------------------------
! 84: us       ! wind in eastward direction [m/s]
! 85: vs       ! wind in northward direction [m/s]
! 86: tm       ! temperature at reference height [kelvin]
! 87: qm       ! specific humidity at reference height [kg/kg]
! 88: prc      ! convective precipitation [mm/s]
! 89: prl      ! large scale precipitation [mm/s]
! 90: pbot     ! atmospheric pressure at the surface [pa]
! 91: frl      ! atmospheric infrared (longwave) radiation [W/m2]
! 92: solar    ! downward solar radiation at surface [W/m2]
!--------------!-------------------------------------------------------

      use precision
      use phycon_module, only: vonkar, stefnc, cpair, rgas, grav
      implicit none

! arguments:
      integer, intent(in) :: lon_points            ! number of longitude points on model grid
      integer, intent(in) :: lat_points            ! number of latitude points on model grid

      integer, intent(in) :: numpatch              ! total number of patches of grids
      integer, intent(in) :: mempatch              ! =max(1,numpatch)
      integer, intent(in) :: ixy_patch(mempatch)   ! patch longitude index
      integer, intent(in) :: jxy_patch(mempatch)   ! patch latitude index
      integer, intent(in) :: itypwat(mempatch)     ! land water type
      real(r8), intent(in) :: wtxy_patch(mempatch) ! patch weight

      integer, INTENT(in) :: nfcon                 ! number of time constant variables
      integer, INTENT(in) :: nforc                 ! number of forcing variables
      integer, intent(in) :: nfldv                 ! number of output variables
      real(r8), INTENT(in) :: fcon(mempatch,nfcon) ! time constant variables
      real(r8), intent(in) :: fldv(mempatch,nfldv) ! output variables
      real(r8), intent(in) :: forcxy(lon_points,lat_points,nforc) ! xy gridded forcing
      real(r8), intent(out) :: fldxy_r(lon_points,lat_points,nfldv) ! xy gridded output
      integer, INTENT(in)  :: msk_patch(mempatch)

! local variables
      integer  :: i,j,k,L                            ! indices
      real(r8) :: swt(lon_points,lat_points)         ! sum of wt
      real(r8) rhoair,thm,th,thv,ur,displa,zldis,hu,ht,hq
      real(r8) z0m,z0h,z0q,us,vs,tm,qm,pbot,psrf
      real(r8) obu,temp1,temp2,temp12m,temp22m
      real(r8) um,thvstar,beta,zii,wc,wc2

!----------------------------------------------------------------------

      fldxy_r(:,:,:) = 0.0
      swt(:,:) = 0.0

!$omp parallel
!$omp do private(i,j,k)
      do k = 1, numpatch
         i = ixy_patch(k)
         j = jxy_patch(k)
!$omp atomic
         swt(i,j)=swt(i,j)+wtxy_patch(k)

![1-69]
         do L = 1, 69
!$omp atomic
            fldxy_r(i,j,L) = fldxy_r(i,j,L) + wtxy_patch(k)*fldv(k,L)
         enddo

![79-83]
         do L = 79, 83
!$omp atomic
            fldxy_r(i,j,L) = fldxy_r(i,j,L) + wtxy_patch(k)*fldv(k,L)
         enddo
      enddo
!$omp end do

      do k = 1, numpatch
         i = ixy_patch(k)
         j = jxy_patch(k)
         if (msk_patch(k)==0) swt(i,j) = 0
      enddo

!$omp do &
!$omp private(i,j) &
!$omp private(z0m,z0h,z0q,displa,hu,ht,hq,zldis,us,vs,tm,qm,pbot,psrf,rhoair) &
!$omp private(thm,th,thv,beta,zii,thvstar,ur,um,wc,wc2,obu,temp1,temp2,temp12m,temp22m)
      do j = 1, lat_points
![1-69]
         do i = 1, lon_points
            if(swt(i,j).gt.0.)then
               fldxy_r(i,j,1:69) = fldxy_r(i,j,1:69)/swt(i,j)
            else
               fldxy_r(i,j,1:69) = -9999.
            endif
         enddo

![70-78]
         do i = 1, lon_points       
            if(swt(i,j).gt.0.)then
               z0m = fldxy_r(i,j,69)
               z0h = fldxy_r(i,j,69)
               z0q = fldxy_r(i,j,69)
               displa = 2./3.*z0m/0.07
       
               hu = max(forcxy(i,j,16),5.+displa)
               ht = max(forcxy(i,j,17),5.+displa)
               hq = max(forcxy(i,j,18),5.+displa)
!lk               zldis = hu-displa
               zldis = hu + z0m
       
               us = forcxy(i,j,3)
               vs = forcxy(i,j,4)
               tm = forcxy(i,j,5)
               qm = forcxy(i,j,6)
               pbot = forcxy(i,j,9)
               psrf = forcxy(i,j,10)
     
               rhoair = (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm)
 
               fldxy_r(i,j,70) = (fldxy_r(i,j,15)/stefnc)**0.25 
               fldxy_r(i,j,71) = sqrt(max(1.e-6_r8,sqrt(fldxy_r(i,j,1)**2+fldxy_r(i,j,2)**2))/rhoair) 
               fldxy_r(i,j,72) = -fldxy_r(i,j,3)/(rhoair*fldxy_r(i,j,71))/cpair
               fldxy_r(i,j,73) = -fldxy_r(i,j,5)/(rhoair*fldxy_r(i,j,71))
 
               thm = tm + 0.0098*ht
               th = tm*(100000./psrf)**(rgas/cpair)
               thv = th*(1.+0.61*qm)       
 
               fldxy_r(i,j,74) = zldis*vonkar*grav&
                   * (fldxy_r(i,j,72)+0.61*th*fldxy_r(i,j,73))&
                   / (fldxy_r(i,j,71)**2*thv)
 
               if(fldxy_r(i,j,74) .ge. 0.)then   !stable
                  fldxy_r(i,j,74) = min(2._r8,max(fldxy_r(i,j,74),1.e-6_r8))
               else                              !unstable
                  fldxy_r(i,j,74) = max(-100._r8,min(fldxy_r(i,j,74),-1.e-6_r8))
               endif
 
               beta = 1.
               zii = 1000.
               thvstar=fldxy_r(i,j,72)+0.61*th*fldxy_r(i,j,73)
               ur = sqrt(us*us+vs*vs)
               if(fldxy_r(i,j,74) .ge. 0.)then
                  um = max(ur,0.1_r8)
               else
                  wc = (-grav*fldxy_r(i,j,71)*thvstar*zii/thv)**(1./3.)
                 wc2 = beta*beta*(wc*wc)
                  um = max(0.1_r8,sqrt(ur*ur+wc2))
               endif
 
               obu = zldis/fldxy_r(i,j,74)
               call moninobuk(hu,ht,hq,displa,z0m,z0h,z0q,&
                    obu,um,fldxy_r(i,j,71),temp1,temp2,temp12m,temp22m,&
                    fldxy_r(i,j,83),fldxy_r(i,j,76),fldxy_r(i,j,77),fldxy_r(i,j,78))
    
!               fldxy_r(i,j,75) = fldxy_r(i,j,74)*vonkar**3*fldxy_r(i,j,71)**2/(temp1*um**2)
               fldxy_r(i,j,75) = fldxy_r(i,j,74)/vonkar*fldxy_r(i,j,71)**2/(temp1*um**2)
               fldxy_r(i,j,75) = min(5._r8,fldxy_r(i,j,75)) 

            else
               fldxy_r(i,j,70:78) = -9999.
            endif
    
         enddo

![79-83]
         do i = 1, lon_points
            if(swt(i,j).gt.0.)then
               fldxy_r(i,j,79:83) = fldxy_r(i,j,79:83)/swt(i,j)
            else
               fldxy_r(i,j,79:83) = -9999.
            endif
         enddo

![84-92]
         do i = 1, lon_points
            if(swt(i,j).gt.0.)then
               fldxy_r(i,j,84) = forcxy(i,j,3)
               fldxy_r(i,j,85) = forcxy(i,j,4)
               fldxy_r(i,j,86) = forcxy(i,j,5)
               fldxy_r(i,j,87) = forcxy(i,j,6)
               fldxy_r(i,j,88) = forcxy(i,j,7)
               fldxy_r(i,j,89) = forcxy(i,j,8)
               fldxy_r(i,j,90) = forcxy(i,j,9)
               fldxy_r(i,j,91) = forcxy(i,j,15)
               fldxy_r(i,j,92) = forcxy(i,j,11)+forcxy(i,j,12)+forcxy(i,j,13)+forcxy(i,j,14)
            else
               fldxy_r(i,j,84:92) = -9999.
            endif
         enddo
      enddo
!$omp end do nowait
!$omp end parallel
  end subroutine vec2xy
! ----------------------------------------------------------------------
! EOP
