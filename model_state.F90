!--F90

!#########################################################################
!# File Name: model_state.F90
!# Created on : 2023-03-22 09:02:38
!# Author: YunhaoFu
!# Last Modified: 2023-04-13 10:16:11
!# Description: works for fortran90 with pre-compile
!#########################################################################

MODULE model_state
    USE precision
    USE netcdf
    USE obs_state, ONLY: idate, o_lat_rever, &
    obs_raw, obs_count_raw, obs_count_r, obs_count, &
    ens_size, max_y_delta, diag_dir, lon_points, lat_points, edges, edgen, edgee, edgew

    IMPLICIT NONE

    PRIVATE
INCLUDE "./m_parameter.h"

    PUBLIC          ::   gridpoint, ens_read, read_invar,&
                         innovation_check, ens_write, &
                         obs_allocate, obs_deallocate, &
                         model_allocate, model_deallocate

    PUBLIC          ::      latixy, longxy, &
                            obs_lon, obs_lat, &
                            omb, y, u, &
                            x_ens, da_cnt , &
                            numpatch, grid2patch_count, grid2patch_start, &
                            mask, itywat, ivt, ixy_patch, jxy_patch, wtxy_patch, &
                            hx
    !--------------------------------------------------!
    !                    temp use                      !
    INTEGER                       ::   idate_(3)
    INTEGER                       ::   ens, nx, ny, nvar, npatch
    INTEGER                       ::   i, j , np, nl, ns, ns1, ns2, mnp
    INTEGER                       ::   ncid, varid, dimid, vid
    REAL                          ::   dx, dy
    CHARACTER(LEN=15)             ::   ens_str
    REAL(r8)                      ::   tmp1, tmp2
    !--------------------------------------------------!
    !              restart file descrption             !
    INTEGER                       ::   numpatch
    INTEGER, ALLOCATABLE          ::   ixy_patch(:)
    INTEGER, ALLOCATABLE          ::   jxy_patch(:)
    REAL(r8),ALLOCATABLE          ::   wtxy_patch(:)
    REAL(r8),ALLOCATABLE          ::   sumwt(:,:)
    INTEGER ,ALLOCATABLE          ::   mask(:,:)
    INTEGER, ALLOCATABLE          ::   itywat(:)
    INTEGER, ALLOCATABLE          ::   ivt(:)
    integer, allocatable          ::   grid2patch(:)
    integer, allocatable          ::   grid2patch_start(:,:)
    integer, allocatable          ::   grid2patch_count(:,:)
    !--------------------------------------------------!
    !               model state patch                  !
    INTEGER                       ::   cnt_s, cnt_e
    REAL(r8),ALLOCATABLE          ::   var_patch(:,:)
    REAL(r8),ALLOCATABLE          ::   ts_patch (:,:,:)
    REAL(r8),ALLOCATABLE          ::   wl_patch (:,:,:)
    REAL(r8),ALLOCATABLE          ::   wi_patch (:,:,:)
    REAL(r8),ALLOCATABLE          ::   tr_patch (:)
    REAL(r8),ALLOCATABLE          ::   x_ens    (:,:,:)
    INTEGER ,ALLOCATABLE          ::   lb_patch (:,:)
    REAL(r8),ALLOCATABLE          ::   sc_patch (:,:)
    !                model state grid                  !
    REAL(r8),PARAMETER            ::   fv = -999.0
    REAL(r8),ALLOCATABLE          ::   longxy(:,:)
    REAL(r8),ALLOCATABLE          ::   latixy(:,:)
    REAL(r8),ALLOCATABLE          ::   hx_ens(:,:,:)
    !--------------------------------------------------!
    !                obs state temp                    !
    REAL(r8),ALLOCATABLE          ::   obs_lon_r(:)
    REAL(r8),ALLOCATABLE          ::   obs_lat_r(:)
    REAL(r8),ALLOCATABLE          ::         y_r(:)
    REAL(r8),ALLOCATABLE          ::         u_r(:)
    REAL(r8),ALLOCATABLE          ::        dt_r(:)
    REAL(r8),ALLOCATABLE          ::        hx_r(:,:)
    REAL(r8),ALLOCATABLE          ::       omb_r(:)
    INTEGER ,ALLOCATABLE          ::   lon_cnt_r(:)
    INTEGER ,ALLOCATABLE          ::   lat_cnt_r(:)

    REAL(r8),ALLOCATABLE          ::   obs_lon(:)
    REAL(r8),ALLOCATABLE          ::   obs_lat(:)
    REAL(r8),ALLOCATABLE          ::         y(:)
    REAL(r8),ALLOCATABLE          ::         u(:)
    REAL(r8),ALLOCATABLE          ::        dt(:)
    REAL(r8),ALLOCATABLE          ::      hx(:,:)
    REAL(r8),ALLOCATABLE          ::       omb(:)
    INTEGER ,ALLOCATABLE          ::   lon_cnt(:)
    INTEGER ,ALLOCATABLE          ::   lat_cnt(:)

    CONTAINS

        SUBROUTINE obs_allocate
            ALLOCATE(obs_lon_r    (obs_count_raw))
            ALLOCATE(obs_lat_r    (obs_count_raw))
            ALLOCATE(      y_r    (obs_count_raw))
            ALLOCATE(      u_r    (obs_count_raw))
            ALLOCATE(     dt_r    (obs_count_raw))
            ALLOCATE(     hx_r(ens_size,obs_count_raw))
            ALLOCATE(    omb_r    (obs_count_raw))
            ALLOCATE(    lon_cnt_r(obs_count_raw))
            ALLOCATE(    lat_cnt_r(obs_count_raw))

            ALLOCATE(obs_lon    (obs_count_raw))
            ALLOCATE(obs_lat    (obs_count_raw))
            ALLOCATE(      y    (obs_count_raw))
            ALLOCATE(      u    (obs_count_raw))
            ALLOCATE(     dt    (obs_count_raw))
            ALLOCATE(hx(ens_size,obs_count_raw))
            ALLOCATE(    omb    (obs_count_raw))
            ALLOCATE(    lon_cnt(obs_count_raw))
            ALLOCATE(    lat_cnt(obs_count_raw))

        ENDSUBROUTINE obs_allocate

        SUBROUTINE model_allocate
            ALLOCATE(sumwt(lon_points,lat_points))
            ALLOCATE(mask(lon_points,lat_points))
            allocate(grid2patch (numpatch))
            allocate(grid2patch_start(lon_points, lat_points))
            allocate(grid2patch_count(lon_points, lat_points))
            ALLOCATE(longxy(lon_points,lat_points))
            ALLOCATE(latixy(lon_points,lat_points))
            ALLOCATE(hx_ens(ens_size,lon_points,lat_points))
        ENDSUBROUTINE model_allocate

        SUBROUTINE model_deallocate
            DEALLOCATE(sumwt)
            DEALLOCATE(mask)
            deallocate(grid2patch)
            deallocate(grid2patch_start)
            deallocate(grid2patch_count)
            DEALLOCATE(longxy)
            DEALLOCATE(latixy)
            DEALLOCATE(hx_ens)
        ENDSUBROUTINE model_deallocate

        SUBROUTINE gridpoint
            IF(o_lat_rever /= m_lat_rever) THEN
                PRINT *, "obs preprocess for model grid point(1 -> reverse; 0 -> no reverse): ", o_lat_rever
                PRINT *, "                   model grid point(1 -> reverse; 0 -> no reverse): ", m_lat_rever  
                STOP 8
            ENDIF
            dx = (edgee-edgew)/(lon_points)
            dy = (edgen-edges)/(lat_points)

            PRINT *, "edgee   edgew   edgen   edges   lon_points   lat_points"
            PRINT *, edgee  ,edgew  ,edgen  ,edges  ,lon_points  ,lat_points
            PRINT *, "dx = ", dx
            PRINT *, "dy = ", dy

            IF(m_lat_rever == 0) THEN
                PRINT *, "lat no reverse"
                DO j = 1, lat_points
                    DO i = 1, lon_points
                    latixy(i,j) = edges + (2*j-1)*dy/2.
                    longxy(i,j) = edgew + (2*i-1)*dx/2.
                    ENDDO
                ENDDO
            ELSE
                PRINT *, "lat do reverse"
                DO j = 1, lat_points
                    DO i = 1, lon_points
                    latixy(i,j) = edgen - (2*j-1)*dy/2.
                    longxy(i,j) = edgew + (2*i-1)*dx/2.
                    ENDDO
                ENDDO
            ENDIF
        ENDSUBROUTINE gridpoint

        SUBROUTINE read_invar
            character(len=256) :: mfile  ! member file
            integer, parameter :: ens_tmp = 1
            logical            :: first_save = .true.
            integer            :: i_save, j_save

            write(mfile, "('bkg.', I3.3)") ens_tmp

            CALL check(NF90_OPEN(TRIM(mfile), nf90_nowrite, ncid))
            CALL check(NF90_INQ_DIMID(ncid,"nlon",dimid))
            CALL check(NF90_INQUIRE_DIMENSION(ncid,dimid,LEN=nx))
            IF(nx /= lon_points) STOP 9
            CALL check(NF90_INQ_DIMID(ncid,"nlat",dimid))
            CALL check(NF90_INQUIRE_DIMENSION(ncid,dimid,LEN=ny))
            IF(ny /= lat_points) STOP 9
            CALL check(NF90_INQ_DIMID(ncid,"numpatch",dimid))
            CALL check(NF90_INQUIRE_DIMENSION(ncid,dimid,LEN=numpatch))

            PRINT *, "----------------------------------------------------"
            PRINT *, "Dimension: lon            lat            numpatch"
            PRINT *, lon_points, lat_points, numpatch
            PRINT *, "----------------------------------------------------"

            ALLOCATE(ixy_patch  (numpatch))
            ALLOCATE(jxy_patch  (numpatch))
            ALLOCATE(wtxy_patch (numpatch))
            ALLOCATE(itywat     (numpatch))
            ALLOCATE(   ivt     (numpatch))
            CALL check(NF90_INQ_VARID(ncid,'ixy_patch',varid))
            CALL check(NF90_GET_VAR(ncid,varid,ixy_patch))
            CALL check(NF90_INQ_VARID(ncid,'jxy_patch',varid))
            CALL check(NF90_GET_VAR(ncid,varid,jxy_patch))
            CALL check(NF90_INQ_VARID(ncid,'wtxy_patch',varid))
            CALL check(NF90_GET_VAR(ncid,varid,wtxy_patch))
            CALL check(NF90_INQ_VARID(ncid,'mask',varid))
            CALL check(NF90_GET_VAR(ncid,varid,mask))
            CALL check(NF90_INQ_VARID(ncid,'fcon',varid))
            CALL check(NF90_GET_VAR(ncid,varid,itywat,&
                start=(/1,invar_int_idx(1)/),count=(/numpatch,1/)))
            CALL check(NF90_GET_VAR(ncid,varid,ivt,&
                start=(/1,invar_int_idx(2)/),count=(/numpatch,1/)))

            CALL check(NF90_CLOSE(ncid))

            grid2patch_count = 0
            do np=1, numpatch
                i = ixy_patch(np)
                j = jxy_patch(np)
                grid2patch_count(i,j) = grid2patch_count(i,j) + 1
                if(first_save .or. (i /= i_save .or. j /= j_save)) then
                    first_save = .false.
                    i_save = i
                    j_save = j
                    grid2patch_start(i,j) = np
                endif
            enddo

            sumwt = 0.0_r8
            do np = 1, numpatch
               i = ixy_patch(np)
               j = jxy_patch(np)
               sumwt(i,j)=sumwt(i,j)+wtxy_patch(np)
            enddo
            do np = 1, numpatch
               i = ixy_patch(np)
               j = jxy_patch(np)
               if (mask(i,j)==0) sumwt(i,j) = 0.0_r8
            enddo

            PRINT *, "longxy : ", MINVAL(longxy)        , MAXVAL(longxy)
            PRINT *, "latixy : ", MINVAL(latixy)        , MAXVAL(latixy)
            PRINT *, "ixy    : ", MINVAL(ixy_patch)     , MAXVAL(ixy_patch)
            PRINT *, "jxy    : ", MINVAL(jxy_patch)     , MAXVAL(jxy_patch)
            PRINT *, "wxy    : ", MINVAL(wtxy_patch)    , MAXVAL(wtxy_patch)
            PRINT *, "sumwt  : ", MINVAL(sumwt)         , MAXVAL(sumwt)
            PRINT *, "ity    : ", MINVAL(itywat)        , MAXVAL(itywat)
            PRINT *, "ivt    : ", MINVAL(ivt)           , MAXVAL(ivt)
            PRINT *, "g2p    : ", MINVAL(grid2patch)    , MAXVAL(grid2patch)
            PRINT *, "g2p_cnt: ", MINVAL(grid2patch_count), MAXVAL(grid2patch_count)
            PRINT *, "g2p_sat: ", MINVAL(grid2patch_start), MAXVAL(grid2patch_start)

        ENDSUBROUTINE read_invar

        SUBROUTINE ens_read
            character(len=256) :: mfile  ! member file
            logical            :: lexist

            ALLOCATE(var_patch  (numpatch          ,SUM(var_length)))
            ALLOCATE(ts_patch   (numpatch,maxsnl+1:nl_soil,ens_size))
            ALLOCATE(wl_patch   (numpatch,maxsnl+1:nl_soil,ens_size))
            ALLOCATE(wi_patch   (numpatch,maxsnl+1:nl_soil,ens_size))
            ALLOCATE(x_ens      (da_cnt,numpatch          ,ens_size))
            ALLOCATE(tr_patch   (numpatch                          ))
            ALLOCATE(lb_patch   (numpatch                 ,ens_size))
            ALLOCATE(sc_patch   (numpatch                 ,ens_size))
            hx_ens   =  0.0_r8
            lb_patch =  0

            DO ens=1,ens_size
                write(mfile, "('bkg.', I3.3)") ens
                PRINT *, 'Read bkg: ', TRIM(mfile)
                INQUIRE(FILE=TRIM(mfile),EXIST=lexist)
                IF(.NOT. lexist) THEN
                    PRINT *, "Not found: ",TRIM(mfile)
                    STOP 10
                ENDIF
                CALL bkg_read(mfile)
            ENDDO
            PRINT *, 'Read bkg over'
        ENDSUBROUTINE ens_read

        SUBROUTINE bkg_read(mfile)
            character(len=*), intent(in) :: mfile

            CALL check(NF90_OPEN(TRIM(mfile), nf90_nowrite, ncid))
            CALL check(NF90_GET_ATT(ncid,nf90_global,'idate',idate_))
            PRINT *, "idate    in rst: ", idate_
            IF(idate_(1) /= idate(1) .OR. &
               idate_(2) /= idate(2) .OR. &
               idate_(3) /= idate(3) ) THEN
                PRINT *, "WARNING:"
                PRINT *, "idate   in NML: ", idate
            ENDIF

            CALL check(NF90_INQ_VARID(ncid,'fvar',varid))
            DO nvar=1, var_count
                cnt_s = SUM(var_length(1:nvar-1))
                cnt_e = SUM(var_length(1:nvar))
                PRINT*, 'put ',nvar,' var in: ', cnt_s+1, ' to ', cnt_e
                CALL check(NF90_GET_VAR(ncid,varid,&
                var_patch(1:numpatch,cnt_s+1:cnt_e),&
                    start=(/1,var_sidx(nvar)/),count=(/numpatch,var_length(nvar)/)))
            ENDDO
            CALL check(NF90_CLOSE(ncid))

            ! read Y
            tr_patch(1:numpatch                     ) = var_patch(1:numpatch,50)
            PRINT *, "trad: ", MINVAL(tr_patch(1:numpatch                     )), MAXVAL(tr_patch(1:numpatch                     ))

            ! transform Y from patch to grid
            do np = 1, numpatch
               i = ixy_patch(np)
               j = jxy_patch(np)
               hx_ens(ens,i,j) = hx_ens(ens,i,j) + wtxy_patch(np) * tr_patch(np)
            enddo
            DO j=1,lat_points
                DO i=1,lon_points
                    IF(sumwt(i,j) .GT. 0.) THEN
                        hx_ens(ens,i,j) = hx_ens(ens,i,j) / sumwt(i,j)
                    ELSE
                        hx_ens(ens,i,j) = fv
                    ENDIF
                ENDDO
            ENDDO

            ! read X
            ts_patch(1:numpatch,maxsnl+1:nl_soil,ens) = var_patch(1:numpatch, 1:15)
            wl_patch(1:numpatch,maxsnl+1:nl_soil,ens) = var_patch(1:numpatch,16:30)
            wi_patch(1:numpatch,maxsnl+1:nl_soil,ens) = var_patch(1:numpatch,31:45)
            x_ens (1,1:numpatch                 ,ens) = var_patch(1:numpatch,46)
            x_ens (3,1:numpatch                 ,ens) = var_patch(1:numpatch,47)
            x_ens (4,1:numpatch                 ,ens) = var_patch(1:numpatch,48)
            sc_patch(1:numpatch                 ,ens) = var_patch(1:numpatch,49)

            DO np=1,numpatch
                ! TODO only concern itywat <=3: soil => 0; urban and built-up => 1; wetland => 2; land ice => 3
                IF(itywat(np) > 3) CYCLE
                DO nl=maxsnl+1,0
                    IF(wl_patch(np,nl,ens)+wi_patch(np,nl,ens) > 0.) lb_patch(np,ens) = lb_patch(np,ens) - 1
                ENDDO
                IF(ts_patch(np,lb_patch(np,ens)+1,ens) /= x_ens(1,np,ens)) THEN
                    PRINT *, "itywat type: ", itywat(np)
                    PRINT *, "ivt    type: ", ivt(np)
                    PRINT *, "t soil at  : ", lb_patch(np,ens)+1," :", ts_patch(np,lb_patch(np,ens)+1,ens)
                    PRINT *, "t_grnd     : ", x_ens(1,np,ens)
                    PRINT *, "t soil whole level: ", ts_patch(np,:,ens)
                    STOP 12
                ENDIF
                x_ens(2,np,ens) = ts_patch(np,lb_patch(np,ens)+0,ens)
                x_ens(5,np,ens) = wl_patch(np,lb_patch(np,ens)+1,ens)
                x_ens(6,np,ens) = wl_patch(np,lb_patch(np,ens)+0,ens)
                x_ens(7,np,ens) = wi_patch(np,lb_patch(np,ens)+1,ens)
                x_ens(8,np,ens) = wi_patch(np,lb_patch(np,ens)+0,ens)
            ENDDO

            PRINT *, "tss[1] : ", MINVAL(x_ens (1,1:numpatch,ens)), MAXVAL(x_ens (1,1:numpatch,ens))
            PRINT *, "tss[2] : ", MINVAL(x_ens (2,1:numpatch,ens)), MAXVAL(x_ens (2,1:numpatch,ens))
            PRINT *, "tsun   : ", MINVAL(x_ens (3,1:numpatch,ens)), MAXVAL(x_ens (3,1:numpatch,ens))
            PRINT *, "tsha   : ", MINVAL(x_ens (4,1:numpatch,ens)), MAXVAL(x_ens (4,1:numpatch,ens))
            PRINT *, "wliq[1]: ", MINVAL(x_ens (5,1:numpatch,ens)), MAXVAL(x_ens (5,1:numpatch,ens))
            PRINT *, "wliq[2]: ", MINVAL(x_ens (6,1:numpatch,ens)), MAXVAL(x_ens (6,1:numpatch,ens))
            PRINT *, "wice[1]: ", MINVAL(x_ens (7,1:numpatch,ens)), MAXVAL(x_ens (7,1:numpatch,ens))
            PRINT *, "wice[2]: ", MINVAL(x_ens (8,1:numpatch,ens)), MAXVAL(x_ens (8,1:numpatch,ens))

        ENDSUBROUTINE bkg_read

        SUBROUTINE ens_write
            character(len=256) :: mfile  ! member file
            logical            :: lexist


            DO ens=1,ens_size
                write(mfile, "('ana.', I3.3)") ens
                PRINT *, 'Write ana: ', TRIM(mfile)
                INQUIRE(FILE=TRIM(mfile),EXIST=lexist)
                IF(.NOT. lexist) THEN
                    PRINT *, "Not found: ",TRIM(mfile)
                    STOP 11
                ENDIF
                CALL ana_write(mfile)
            ENDDO
            PRINT *, 'Write ana over'

            DEALLOCATE(ixy_patch  )
            DEALLOCATE(jxy_patch  )
            DEALLOCATE(wtxy_patch )
            DEALLOCATE(itywat     )
            DEALLOCATE(ivt        )
            DEALLOCATE(var_patch  )
            DEALLOCATE(ts_patch   )
            DEALLOCATE(wl_patch   )
            DEALLOCATE(wi_patch   )
            DEALLOCATE(x_ens      )
            DEALLOCATE(tr_patch   )
            DEALLOCATE(lb_patch   )
            DEALLOCATE(sc_patch   )

        ENDSUBROUTINE ens_write

        SUBROUTINE ana_write(mfile)
           character(len=*), intent(in) :: mfile
           real(r8)                     :: scv

            CALL check(NF90_OPEN(TRIM(mfile), nf90_write, ncid))

            CALL check(NF90_GET_ATT(ncid,nf90_global,'idate',idate_))
            PRINT *, "idate    in rst: ", idate_
            IF(idate_(1) /= idate(1) .OR. &
               idate_(2) /= idate(2) .OR. &
               idate_(3) /= idate(3) ) THEN
                PRINT *, "WARNING:"
                PRINT *, "idate   in NML: ", idate
            ENDIF

            
            DO np=1,numpatch
                ! TODO only concern itywat <=3: soil => 0; urban and built-up => 1; wetland => 2; land ice => 3
                IF(itywat(np) > 3) CYCLE
                scv = 0.0_r8
                DO nl=maxsnl+1,0
                    scv = scv + wl_patch(np,nl,ens) + wi_patch(np,nl,ens)
                ENDDO
                sc_patch(np,ens) = scv
                ts_patch(np,lb_patch(np,ens)+1,ens) = x_ens(1,np,ens)
                ts_patch(np,lb_patch(np,ens)+0,ens) = x_ens(2,np,ens)
                wl_patch(np,lb_patch(np,ens)+1,ens) = x_ens(5,np,ens)
                wl_patch(np,lb_patch(np,ens)+0,ens) = x_ens(6,np,ens)
                wi_patch(np,lb_patch(np,ens)+1,ens) = x_ens(7,np,ens)
                wi_patch(np,lb_patch(np,ens)+0,ens) = x_ens(8,np,ens)
                !IF(ts_patch(np,lb_patch(np,ens)+1,ens) /= x_ens(1,np,ens)) THEN
                !    PRINT *, "itywat type: ", itywat(np)
                !    PRINT *, "ivt    type: ", ivt(np)
                !    PRINT *, "t soil at  : ", lb_patch(np,ens)+1," :", ts_patch(np,lb_patch(np,ens)+1,ens)
                !    PRINT *, "t_grnd     : ", x_ens(1,np,ens)
                !    PRINT *, "t soil whole level: ", ts_patch(np,:,ens)
                !    STOP 14
                !ELSE
                !    IF(ts_patch(np,lb_patch(np,ens)+1,ens) - x_ens(1,np,ens) /= 0.0_r8) STOP 30
                !    ts_patch(np,lb_patch(np,ens)+1,ens) = x_ens(1,np,ens)
                !ENDIF
            ENDDO

            var_patch(1:numpatch,1:15)  = ts_patch(1:numpatch,maxsnl+1:nl_soil,ens)
            var_patch(1:numpatch,16:30) = wl_patch(1:numpatch,maxsnl+1:nl_soil,ens)
            var_patch(1:numpatch,31:45) = wi_patch(1:numpatch,maxsnl+1:nl_soil,ens)
            var_patch(1:numpatch,46)    = x_ens (1,1:numpatch                 ,ens)
            var_patch(1:numpatch,47)    = x_ens (3,1:numpatch                 ,ens)
            var_patch(1:numpatch,48)    = x_ens (4,1:numpatch                 ,ens)
            var_patch(1:numpatch,49)    = sc_patch(1:numpatch                 ,ens)

            PRINT *, "tss[1] : ", MINVAL(x_ens (1,1:numpatch,ens)), MAXVAL(x_ens (1,1:numpatch,ens))
            PRINT *, "tss[2] : ", MINVAL(x_ens (2,1:numpatch,ens)), MAXVAL(x_ens (2,1:numpatch,ens))
            PRINT *, "tsun   : ", MINVAL(x_ens (3,1:numpatch,ens)), MAXVAL(x_ens (3,1:numpatch,ens))
            PRINT *, "tsha   : ", MINVAL(x_ens (4,1:numpatch,ens)), MAXVAL(x_ens (4,1:numpatch,ens))
            PRINT *, "wliq[1]: ", MINVAL(x_ens (5,1:numpatch,ens)), MAXVAL(x_ens (5,1:numpatch,ens))
            PRINT *, "wliq[2]: ", MINVAL(x_ens (6,1:numpatch,ens)), MAXVAL(x_ens (6,1:numpatch,ens))
            PRINT *, "wice[1]: ", MINVAL(x_ens (7,1:numpatch,ens)), MAXVAL(x_ens (7,1:numpatch,ens))
            PRINT *, "wice[2]: ", MINVAL(x_ens (8,1:numpatch,ens)), MAXVAL(x_ens (8,1:numpatch,ens))
            
            CALL check(NF90_INQ_VARID(ncid,'fvar',varid))
            DO nvar=1, var_count
                IF(var_write(nvar) == 1) THEN
                    cnt_s = SUM(var_length(1:nvar-1))
                    cnt_e = SUM(var_length(1:nvar))
                    PRINT*, 'write ',nvar,' var in: ', cnt_s+1, ' to ', cnt_e
                    CALL check(NF90_PUT_VAR(ncid,varid,&
                    var_patch(1:numpatch,cnt_s+1:cnt_e),&
                        start=(/1,var_sidx(nvar)/),count=(/numpatch,var_length(nvar)/)))
                ENDIF
            ENDDO

            CALL check(NF90_CLOSE(ncid))

        ENDSUBROUTINE ana_write

        SUBROUTINE innovation_check
            obs_count   = 0
            obs_count_r = 0
            DO ns=1,obs_count_raw
            
                ! fillvalue for hx check
                IF(hx_ens(1,obs_raw(ns)%model_lon_count,obs_raw(ns)%model_lat_count) == fv) CYCLE

                tmp1 = SUM(hx_ens(1:ens_size,obs_raw(ns)%model_lon_count,obs_raw(ns)%model_lat_count)) / REAL(ens_size,KIND=r8)
                tmp2 = obs_raw(ns)%y - tmp1

                obs_count_r = obs_count_r + 1
                obs_lon_r        (obs_count_r)   =   obs_raw(ns)%obs_lon
                obs_lat_r        (obs_count_r)   =   obs_raw(ns)%obs_lat
                lon_cnt_r        (obs_count_r)   =   obs_raw(ns)%model_lon_count
                lat_cnt_r        (obs_count_r)   =   obs_raw(ns)%model_lat_count
                y_r              (obs_count_r)   =   obs_raw(ns)%y
                u_r              (obs_count_r)   =   obs_raw(ns)%u
                dt_r             (obs_count_r)   =   obs_raw(ns)%dt
                hx_r  (1:ens_size,obs_count_r)   =   hx_ens(1:ens_size,obs_raw(ns)%model_lon_count,obs_raw(ns)%model_lat_count)
                omb_r            (obs_count_r)   =   tmp2

                ! innovation check
                IF(ABS(tmp2) > max_y_delta) CYCLE

                obs_count = obs_count + 1
                obs_lon        (obs_count)   =   obs_raw(ns)%obs_lon
                obs_lat        (obs_count)   =   obs_raw(ns)%obs_lat
                lon_cnt        (obs_count)   =   obs_raw(ns)%model_lon_count
                lat_cnt        (obs_count)   =   obs_raw(ns)%model_lat_count
                y              (obs_count)   =   obs_raw(ns)%y
                u              (obs_count)   =   obs_raw(ns)%u
                dt             (obs_count)   =   obs_raw(ns)%dt
                hx  (1:ens_size,obs_count)   =   hx_ens(1:ens_size,obs_raw(ns)%model_lon_count,obs_raw(ns)%model_lat_count)
                omb            (obs_count)   =   tmp2
            ENDDO

            PRINT *, "obs count raw :    ",obs_count_raw
            PRINT *, "fillvalue for hx : ",obs_count_raw - obs_count_r 
            PRINT *, "----------------------------------------------------------"
            PRINT *, "obs count_r:", obs_count_r
            PRINT *, "obs lon_r  :  ", MINVAL(      obs_lon_r(1:obs_count)), MAXVAL(obs_lon_r(1:obs_count))
            PRINT *, "obs lat_r  :  ", MINVAL(      obs_lat_r(1:obs_count)), MAXVAL(obs_lat_r(1:obs_count))
            PRINT *, "obs y_r    :  ", MINVAL(            y_r(1:obs_count)), MAXVAL(      y_r(1:obs_count))
            PRINT *, "obs u_r    :  ", MINVAL(            u_r(1:obs_count)), MAXVAL(      u_r(1:obs_count))
            PRINT *, "obs dtime_r:  ", MINVAL(           dt_r(1:obs_count)), MAXVAL(     dt_r(1:obs_count))
            PRINT *, "hx_r       :  ", MINVAL(hx_r(1:ens_size,1:obs_count)), MAXVAL(hx_r(1:ens_size,1:obs_count))
            PRINT *, "omb_r      :  ", MINVAL(          omb_r(1:obs_count)), MAXVAL(    omb_r(1:obs_count))
            PRINT *, "----------------------------------------------------------"
            PRINT *, "obs count:  ", obs_count
            PRINT *, "obs lon  :  ", MINVAL(      obs_lon(1:obs_count)), MAXVAL(obs_lon(1:obs_count))
            PRINT *, "obs lat  :  ", MINVAL(      obs_lat(1:obs_count)), MAXVAL(obs_lat(1:obs_count))
            PRINT *, "obs y    :  ", MINVAL(            y(1:obs_count)), MAXVAL(      y(1:obs_count))
            PRINT *, "obs u    :  ", MINVAL(            u(1:obs_count)), MAXVAL(      u(1:obs_count))
            PRINT *, "obs dtime:  ", MINVAL(           dt(1:obs_count)), MAXVAL(     dt(1:obs_count))
            PRINT *, "hx       :  ", MINVAL(hx(1:ens_size,1:obs_count)), MAXVAL(hx(1:ens_size,1:obs_count))
            PRINT *, "omb      :  ", MINVAL(          omb(1:obs_count)), MAXVAL(    omb(1:obs_count))

        ENDSUBROUTINE innovation_check

        SUBROUTINE obs_deallocate
            DEALLOCATE(obs_lon_r)
            DEALLOCATE(obs_lat_r)
            DEALLOCATE(      y_r)
            DEALLOCATE(      u_r)
            DEALLOCATE(     dt_r)
            DEALLOCATE(     hx_r)
            DEALLOCATE(    omb_r)
            DEALLOCATE( lon_cnt_r)
            DEALLOCATE( lat_cnt_r)

            DEALLOCATE(obs_lon)
            DEALLOCATE(obs_lat)
            DEALLOCATE(      y)
            DEALLOCATE(      u)
            DEALLOCATE(     dt)
            DEALLOCATE(     hx)
            DEALLOCATE(    omb)
            DEALLOCATE( lon_cnt)
            DEALLOCATE( lat_cnt)

        ENDSUBROUTINE obs_deallocate

        SUBROUTINE check(status, str)
          INTEGER, INTENT(in) :: status
          CHARACTER(*), OPTIONAL, INTENT(in) :: str
          CHARACTER(:), ALLOCATABLE :: str2
          IF(status /= nf90_noerr) THEN
             str2=TRIM(nf90_strerror(status))

             IF(PRESENT(str)) THEN
                str2=str2//": "//TRIM(str)
             END IF

             STOP 3
          END IF
        END SUBROUTINE check

ENDMODULE model_state
