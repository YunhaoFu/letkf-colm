!--F90

!#########################################################################
!# File Name: letkf_main.F90
!# Created on : 2023-03-31 21:38:50
!# Author: YunhaoFu
!# Last Modified: 2023-04-13 08:29:02
!# Description: works for fortran90 with pre-compile
!#########################################################################
MODULE letkf_main

    USE precision
    USE model_state, ONLY: longxy, latixy, &
                         x_ens, da_cnt, &
                         obs_lon, obs_lat, omb, u, hx, numpatch, &
                         grid2patch_count, grid2patch_start, mask, itywat, &
                         ixy_patch, jxy_patch, wtxy_patch, ivt

    USE obs_state,   ONLY: obs_count, diag_dir, ens_size, &
                         max_y_delta, save_diag, radius, infl, lon_points, lat_points
    USE kdtree
    USE letkf_core
    USE netcdf
    USE localization
    IMPLICIT NONE

    PRIVATE
INCLUDE './letkf.h'

    ! public subroutines
    PUBLIC            ::             letkf_ini, letkf_run, letkf_final

    !------------------------------ diagnostic use ---------------------------------!
    INTEGER,  ALLOCATABLE                 ::       diag_obscount(:,:)               !
    REAL(r8), ALLOCATABLE                 ::       diag_maxdist (:,:)               !
    REAL(r8), ALLOCATABLE                 ::       diag_rloc(:,:,:)                 ! clear description for each variables
    REAL(r8), ALLOCATABLE                 ::       diag_hdxb(:,:,:,:)               ! is shown in the nc write line
    REAL(r8), ALLOCATABLE                 ::       diag_error(:,:,:)                ! can be infered as #MARK1
    REAL(r8), ALLOCATABLE                 ::       diag_dep(:,:,:)                  ! 
    REAL(r8), ALLOCATABLE                 ::       diag_lon(:,:,:)                  !
    REAL(r8), ALLOCATABLE                 ::       diag_lat(:,:,:)                  !
    !----------------------------- localization use --------------------------------!
    INTEGER , PARAMETER                   ::       sp=kind(0.0e0)
    REAL(r8)                              ::       max_search_dist
    REAL(r8)                              ::       dist
    REAL(r8)                              ::       rd
    !--------------------------------- kd tree -------------------------------------!
    TYPE(kd_root)                         ::       obs_tree                         ! kd tree for obs
    !------------------------------- loop for ij -----------------------------------!
    INTEGER                               ::       obs_ij_cnt
    INTEGER                               ::       obs_ij_idx(max_obs)
    REAL(sp)                              ::       obs_ij_dist(max_obs)
    REAL(r8),ALLOCATABLE                  ::       obs_ij_hdxb(:,:)
    REAL(r8)                              ::       obs_ij_rdiag(max_obs)
    REAL(r8)                              ::       obs_ij_dep(max_obs)
    REAL(r8)                              ::       obs_ij_lon(max_obs)
    REAL(r8)                              ::       obs_ij_lat(max_obs)
    INTEGER                               ::       obs_lg_cnt
    REAL(r8)                              ::       obs_lg_dist(max_obs)
    REAL(r8),ALLOCATABLE                  ::       obs_lg_hdxb(:,:)
    REAL(r8)                              ::       obs_lg_rdiag(max_obs)
    REAL(r8)                              ::       obs_lg_dep(max_obs)
    REAL(r8)                              ::       obs_lg_rloc(max_obs)
    REAL(r8)                              ::       obs_lg_lon(max_obs)
    REAL(r8)                              ::       obs_lg_lat(max_obs)
    REAL(r8),ALLOCATABLE                  ::       trans(:,:)
    !---------------------------------- nc file ------------------------------------!
    INTEGER                               ::       ncid, vid
    INTEGER                               ::       i, j, nx, ny, ns, np, mobs, mnp, npatch
    !----------------------------- local components --------------------------------!
    REAL(r8), ALLOCATABLE                 ::       bkg(:,:,:)
    REAL(r8), ALLOCATABLE                 ::       bkg_mean(:,:)
    INTEGER                               ::       nvar, ens
    INTEGER                               ::       idx
    CHARACTER(LEN=1024)                   ::       temp_name
    !CHARACTER(LEN=30)                     ::       temp_name_

    CONTAINS

        SUBROUTINE letkf_ini

            ALLOCATE(obs_ij_hdxb(ens_size,max_obs))
            ALLOCATE(obs_lg_hdxb(ens_size,max_obs))
            ALLOCATE(trans(ens_size,ens_size))

            DO ns=1,obs_count
                hx(1:ens_size,ns) = hx(1:ens_size,ns) - SUM(hx(1:ens_size,ns)) / REAL(ens_size,KIND=r8)
            ENDDO

            CALL kd_init(obs_tree, REAL(obs_lon(1:obs_count),KIND=sp), REAL(obs_lat(1:obs_count),KIND=sp))
            
            ALLOCATE(bkg_mean(da_cnt,numpatch))
            ALLOCATE(bkg     (da_cnt,numpatch,ens_size))

            DO nvar=1,da_cnt
                DO np=1,numpatch
                    bkg_mean(nvar,np) = SUM(x_ens(nvar,np,1:ens_size)) / REAL(ens_size,KIND=r8)
                    bkg(nvar,np,1:ens_size) = x_ens(nvar,np,1:ens_size) - bkg_mean(nvar,np)
                ENDDO
            ENDDO

            x_ens(:,:,:) = bkg(:,:,:)

            IF(save_diag) THEN
                ALLOCATE(diag_obscount       (lon_points,lat_points))
                ALLOCATE(diag_maxdist        (lon_points,lat_points)) 
                ALLOCATE(diag_rloc           (lon_points,lat_points,max_obs))
                ALLOCATE(diag_hdxb  (ens_size,lon_points,lat_points,max_obs))
                ALLOCATE(diag_error          (lon_points,lat_points,max_obs))
                ALLOCATE(diag_dep            (lon_points,lat_points,max_obs))
                ALLOCATE(diag_lon            (lon_points,lat_points,max_obs))
                ALLOCATE(diag_lat            (lon_points,lat_points,max_obs))
                diag_obscount(:,:) = 0
            ENDIF

        ENDSUBROUTINE letkf_ini

        SUBROUTINE letkf_run

            PRINT *, "============================================================"
            PRINT *, "                    LETKF start run                         "

                !TODO use constant localization radius in regional DA
                !dist = radius(1)
                !max_search_dist = dist * stdev2max
                !IF(save_diag)  diag_maxdist(:,:) =  max_search_dist

            lon_loop: DO i=1,lon_points
                lat_loop: DO j=1,lat_points

                    dist=get_dist(latixy(i,j),radius)
                    max_search_dist = dist * stdev2max
                    IF(save_diag)  diag_maxdist(i,j) =  max_search_dist

                    IF(mask(i,j) == 0) CYCLE
                    IF(grid2patch_count(i,j) .LE. 0) CYCLE

                    ! get search radius
                    CALL kd_search_radius(obs_tree, REAL(longxy(i,j),KIND=sp), REAL(latixy(i,j),KIND=sp),&
                       REAL(max_search_dist,KIND=sp), obs_ij_idx, obs_ij_dist, obs_ij_cnt, .FALSE.) 

                    IF(obs_ij_cnt > max_obs) STOP 13
                    IF(obs_ij_cnt .LE. 0) CYCLE

                    ! prepare obs when there is one "nearing" the ij grid point
                    DO ns=1,obs_ij_cnt
                        idx                = obs_ij_idx(ns)
                        obs_ij_hdxb (:,ns) = hx(:,idx)
                        obs_ij_rdiag  (ns) = u(idx)
                        obs_ij_dep    (ns) = omb(idx)
                        obs_ij_lon    (ns) = obs_lon(idx)
                        obs_ij_lat    (ns) = obs_lat(idx)
                    END DO

                    obs_lg_cnt = 0
                    obs_loop: DO ns=1,obs_ij_cnt
                        rd = letkf_loc_gc(obs_ij_dist(ns),dist)
                        IF (rd > 1.0e-4) THEN
                            obs_lg_cnt                = obs_lg_cnt + 1
                            obs_lg_rloc  (obs_lg_cnt) = rd
                            obs_lg_hdxb(:,obs_lg_cnt) = obs_ij_hdxb (:,ns)
                            obs_lg_rdiag (obs_lg_cnt) = obs_ij_rdiag  (ns)
                            obs_lg_dep   (obs_lg_cnt) = obs_ij_dep    (ns)
                            obs_lg_lon   (obs_lg_cnt) = obs_ij_lon    (ns)
                            obs_lg_lat   (obs_lg_cnt) = obs_ij_lat    (ns)
                        END IF
                    END DO obs_loop

                    IF (save_diag) THEN
                        diag_obscount(i,j)                =   obs_lg_cnt
                        diag_rloc    (i,j,1:obs_lg_cnt)   =   obs_lg_rloc  (1:obs_lg_cnt)
                        diag_hdxb  (:,i,j,1:obs_lg_cnt)   =   obs_lg_hdxb(:,1:obs_lg_cnt)
                        diag_error   (i,j,1:obs_lg_cnt)   =   obs_lg_rdiag (1:obs_lg_cnt)
                        diag_dep     (i,j,1:obs_lg_cnt)   =   obs_lg_dep   (1:obs_lg_cnt)
                        diag_lon     (i,j,1:obs_lg_cnt)   =   obs_lg_lon   (1:obs_lg_cnt)
                        diag_lat     (i,j,1:obs_lg_cnt)   =   obs_lg_lat   (1:obs_lg_cnt)
                    END IF

                    IF (obs_lg_cnt > 0) THEN
                        !calculate trans matrix in ij grid point format
                        CALL letkf_core_solve(ens_size, obs_lg_cnt, obs_lg_hdxb(:,:obs_lg_cnt),&
                            obs_lg_rdiag(:obs_lg_cnt), obs_lg_rloc(:obs_lg_cnt), &
                            obs_lg_dep(:obs_lg_cnt), infl, trans)

                        !apply the trans calculated above into different var within the same lat lon, i.e., novertical localization
                        DO nvar=1,da_cnt
                            DO np=grid2patch_start(i,j),grid2patch_start(i,j)+grid2patch_count(i,j)-1
                                ! TODO only concern itywat <=3: soil => 0; urban and built-up => 1; wetland => 2; land ice => 3
                                IF(itywat(np) > 3) CYCLE
                                CALL dgemm('n', 'n', 1, ens_size, ens_size, 1.0_r8, bkg(nvar,np,1:ens_size), &
                                    1, trans, ens_size, 0.0_r8, x_ens(nvar,np,1:ens_size), 1)
                            ENDDO
                        ENDDO
                    END IF
            
                ENDDO lat_loop
            
            ENDDO lon_loop

            DO ens=1,ens_size
                  x_ens(:,:,ens) &
                = x_ens(:,:,ens) + bkg_mean(:,:)
            ENDDO

            PRINT *, "                     LETKF run over                         "
            PRINT *, "============================================================"

            CALL kd_free(obs_tree)

        ENDSUBROUTINE letkf_run

        SUBROUTINE letkf_final

            DEALLOCATE(bkg)
            DEALLOCATE(bkg_mean)

            DEALLOCATE(obs_ij_hdxb)
            DEALLOCATE(obs_lg_hdxb)
            DEALLOCATE(trans)

            IF(save_diag) THEN
                DEALLOCATE(diag_obscount)
                DEALLOCATE(diag_maxdist )
                DEALLOCATE(diag_rloc    )
                DEALLOCATE(diag_hdxb    )
                DEALLOCATE(diag_error   )
                DEALLOCATE(diag_dep     )
                DEALLOCATE(diag_lon     )
                DEALLOCATE(diag_lat     )
            ENDIF

        ENDSUBROUTINE letkf_final

        SUBROUTINE check(status, str)
          INTEGER, INTENT(in) :: status
          CHARACTER(*), OPTIONAL, INTENT(in) :: str
          CHARACTER(:), ALLOCATABLE :: str2
          IF(status /= nf90_noerr) THEN
             str2=TRIM(nf90_strerror(status))

             IF(PRESENT(str)) THEN
                str2=str2//": "//TRIM(str)
             END IF

             STOP 2
          END IF
        END SUBROUTINE check

ENDMODULE letkf_main
