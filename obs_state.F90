!--F90

!#########################################################################
!# File Name: obs_state.F90
!# Created on : 2023-03-22 20:05:33
!# Author: YunhaoFu
!# Last Modified: 2023-05-22 13:15:07
!# Description: works for fortran90 with pre-compile
!#########################################################################

MODULE obs_state
    USE precision
    USE netcdf

    IMPLICIT NONE
    PRIVATE

INCLUDE './o_parameter.h'

    PUBLIC  ::  obs_ini, obs_final

    PUBLIC  ::  idate, o_lat_rever, obs_raw, &
                obs_count_raw, obs_count_r, obs_count, &
                ens_size, max_y_delta, save_diag, diag_dir, &
                radius, infl, lon_points, lat_points, edges, edgen, edgee, edgew

    !------------------ NML file ----------------------!
    INTEGER                       ::   idate(3)        ! input idate
    INTEGER                       ::   lon_points      ! lon points in model
    INTEGER                       ::   lat_points      ! lat points in model
    REAL(r8)                      ::   edges           ! south edge in model  
    REAL(r8)                      ::   edgen           ! north edge in model
    REAL(r8)                      ::   edgew           ! west  edge in model
    REAL(r8)                      ::   edgee           ! east  edge in model
    INTEGER                       ::   ens_size        ! ensemble size
    REAL(r8)                      ::   max_y_delta     ! innovation threshold
    LOGICAL                       ::   save_diag       ! save diagnostic or not
    CHARACTER(LEN=1024)           ::   diag_dir        ! diagnostic output absolute dir
    REAL(r8)                      ::   radius(2)       ! localization radius [m]
    REAL(r8)                      ::   infl            ! inflation
    !--------------------------------------------------!
    !                   temp use                       !
    INTEGER                       ::   ncid
    INTEGER                       ::   dimid
    INTEGER                       ::   varid
    INTEGER,  ALLOCATABLE         ::   tmp_int_1(:)
    REAL(r8), ALLOCATABLE         ::   tmp_real_1(:)
    INTEGER                       ::   ns
    !--------------------------------------------------!
    !                 obs state                        !
    INTEGER                       ::   obs_count_raw   ! obs count pass dtime windows, flag, fillvalue, etc.
    INTEGER                       ::   obs_count_r     ! obs count pass above + hx which is not fillvalue
    INTEGER                       ::   obs_count       ! obs count pass above + innovation threshold
    TYPE             ::       obs_register
        REAL(r8)     ::       obs_lon
        REAL(r8)     ::       obs_lat
        INTEGER      ::       model_lon_count
        INTEGER      ::       model_lat_count
        REAL(r8)     ::       y
        REAL(r8)     ::       u
        REAL(r8)     ::       dt
    ENDTYPE
    TYPE(obs_register), ALLOCATABLE :: obs_raw(:)

    CONTAINS

        SUBROUTINE obs_ini
            NAMELIST /description/  idate    , &
                                    lon_points,&
                                    lat_points,&
                                    edges    , &
                                    edgen    , &
                                    edgew    , &
                                    edgee    , &
                                    ens_size , &
                                    max_y_delta, &
                                    save_diag, &
                                    diag_dir , &
                                    radius   , &
                                    infl

            READ(5,nml=description)
            CLOSE(5)
            PRINT *, "idate in NML: ",idate

            CALL obs_read()
            !DO ns=1,obs_count_raw
            !    PRINT*, "lon: ", obs_raw(ns)%obs_lon, " count: ", obs_raw(ns)%model_lon_count
            !    PRINT*, "lat: ", obs_raw(ns)%obs_lat, " count: ", obs_raw(ns)%model_lat_count
            !ENDDO
            PRINT *, "Read obs nc over"
        ENDSUBROUTINE obs_ini

        SUBROUTINE obs_read
            CALL check(NF90_OPEN("obs.nc", nf90_nowrite, ncid))
            CALL check(NF90_INQ_DIMID(ncid,"station",dimid))
            CALL check(NF90_INQUIRE_DIMENSION(ncid,dimid,LEN=obs_count_raw))
            PRINT *, "obs raw count: ",obs_count_raw
            ALLOCATE(obs_raw(obs_count_raw))
            ALLOCATE(tmp_int_1(obs_count_raw))
            ALLOCATE(tmp_real_1(obs_count_raw))
            CALL check(NF90_INQ_VARID(ncid,'lon',varid))
            CALL check(NF90_GET_VAR(ncid,varid,tmp_real_1))
            PRINT *, "obs lon: ", MINVAL(tmp_real_1), MAXVAL(tmp_real_1)
            DO ns=1,obs_count_raw
                obs_raw(ns)%obs_lon = tmp_real_1(ns)
            ENDDO
            CALL check(NF90_INQ_VARID(ncid,'lat',varid))
            CALL check(NF90_GET_VAR(ncid,varid,tmp_real_1))
            PRINT *, "obs lat: ", MINVAL(tmp_real_1), MAXVAL(tmp_real_1)
            DO ns=1,obs_count_raw
                obs_raw(ns)%obs_lat = tmp_real_1(ns)
            ENDDO
            CALL check(NF90_INQ_VARID(ncid,'lon_count',varid))
            CALL check(NF90_GET_VAR(ncid,varid,tmp_int_1))
            PRINT *, "obs lon count: ", MINVAL(tmp_int_1), MAXVAL(tmp_int_1)
            DO ns=1,obs_count_raw
                obs_raw(ns)%model_lon_count = tmp_int_1(ns)
            ENDDO
            CALL check(NF90_INQ_VARID(ncid,'lat_count',varid))
            CALL check(NF90_GET_VAR(ncid,varid,tmp_int_1))
            PRINT *, "obs lat count: ", MINVAL(tmp_int_1), MAXVAL(tmp_int_1)
            DO ns=1,obs_count_raw
                obs_raw(ns)%model_lat_count = tmp_int_1(ns)
            ENDDO
            CALL check(NF90_INQ_VARID(ncid,'lst',varid))
            CALL check(NF90_GET_VAR(ncid,varid,tmp_real_1))
            PRINT *, "obs: ", MINVAL(tmp_real_1), MAXVAL(tmp_real_1)
            DO ns=1,obs_count_raw
                obs_raw(ns)%y = tmp_real_1(ns)
            ENDDO
            CALL check(NF90_INQ_VARID(ncid,'lst_error',varid))
            CALL check(NF90_GET_VAR(ncid,varid,tmp_real_1))
            PRINT *, "obs error: ", MINVAL(tmp_real_1), MAXVAL(tmp_real_1)
            DO ns=1,obs_count_raw
                obs_raw(ns)%u = tmp_real_1(ns)**2.0_r8
            ENDDO
            CALL check(NF90_INQ_VARID(ncid,'dtime',varid))
            CALL check(NF90_GET_VAR(ncid,varid,tmp_real_1))
            PRINT *, "obs dtime: ", MINVAL(tmp_real_1), MAXVAL(tmp_real_1)
            DO ns=1,obs_count_raw
                obs_raw(ns)%dt = tmp_real_1(ns)
            ENDDO
            CALL check(NF90_CLOSE(ncid))
        ENDSUBROUTINE obs_read

        SUBROUTINE obs_final
                DEALLOCATE(tmp_int_1)
                DEALLOCATE(tmp_real_1)
                DEALLOCATE(obs_raw)
                PRINT *, "obs count: ", obs_count
                IF(obs_count .LE. 0) THEN
                    PRINT *, ''
                    PRINT *, "            No observations pass the innovation check"
                    PRINT *, "=================================================================="
                    PRINT *, "             FYH-LETKF over:  STATE = no  obs, no da    "
                    PRINT *, "=================================================================="
                ELSE
                    PRINT *, ''
                    PRINT *, "             Has observations pass the innovation check"
                    PRINT *, "=================================================================="
                    PRINT *, "             FYH-LETKF over:  STATE = has  obs, da done    "
                    PRINT *, "=================================================================="
                ENDIF
        ENDSUBROUTINE obs_final

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

ENDMODULE obs_state
