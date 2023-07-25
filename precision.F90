!--F90

!#########################################################################
!# File Name: precision.F90
!# Created on : 2023-03-22 20:06:01
!# Author: YunhaoFu
!# Last Modified: 2023-03-22 20:06:12
!# Description: works for fortran90 with pre-compile
!#########################################################################

MODULE precision
    IMPLICIT NONE
    PRIVATE

    PUBLIC :: r8
    !TODO SELECTED_REAL_KIND(8) is equal to SELECTED_REAL_KIND(12)
    INTEGER, PARAMETER:: r8 = SELECTED_REAL_KIND(8)

END MODULE precision
