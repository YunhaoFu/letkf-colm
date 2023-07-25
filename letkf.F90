!--F90

!#########################################################################
!# File Name: letkf.F90
!# Created on : 2023-03-22 09:01:15
!# Author: YunhaoFu
!# Last Modified: 2023-04-06 20:30:55
!# Description: works for fortran90 with pre-compile
!#########################################################################

PROGRAM letkf
    USE obs_state
    USE model_state
    USE letkf_main
    IMPLICIT NONE

    CALL obs_ini()

        CALL obs_allocate()
        CALL model_allocate()
        CALL gridpoint()
        CALL read_invar()
        CALL ens_read()
        CALL innovation_check()
        IF(obs_count .LE. 0) THEN
            CALL obs_deallocate()
            CALL model_deallocate()
            CALL obs_final()
        ELSE
            CALL letkf_ini()
            CALL letkf_run()
            CALL letkf_final()

            CALL ens_write()
            CALL obs_deallocate()
            CALL model_deallocate()
            CALL obs_final()
        ENDIF

ENDPROGRAM letkf

