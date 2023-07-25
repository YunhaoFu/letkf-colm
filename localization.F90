MODULE localization
    USE precision
    IMPLICIT NONE

    PRIVATE

    ! public subroutines
    PUBLIC    ::     get_dist, letkf_loc_gc

    INTEGER, PARAMETER  ::   sp=kind(0.0e0)

    CONTAINS

        FUNCTION get_dist(lat,r) RESULT(dist)

            REAL(r8), INTENT(IN)    :: lat
            REAL(r8), INTENT(IN)    :: r(2)
            REAL(sp)                :: dist
            dist = ( (ABS(lat)-0.0)*r(2)+(90.0-ABS(lat))*r(1) ) / 90.0

        END FUNCTION get_dist

        PURE FUNCTION letkf_loc_gc(z, L) RESULT(res)
            REAL(sp), INTENT(in) :: z  !< value to localize
            REAL(r8), INTENT(in) :: L  !< the equivalent to the Gaussian standard deviation
            REAL(r8) :: res
            REAL(r8) :: c
            REAL(r8) :: abs_z, z_c

            c = L / SQRT(0.3_r8)
            abs_z = ABS(z)
            z_c = abs_z / c

            IF (abs_z >= 2*c) THEN
               res = 0.0
            ELSEIF (abs_z > c) THEN
               res = &
                      0.08333_r8 * z_c**5 &
                    - 0.50000_r8 * z_c**4 &
                    + 0.62500_r8 * z_c**3 &
                    + 1.66667_r8 * z_c**2 &
                    - 5.00000_r8 * z_c &
                    + 4_r8 &
                    - 0.66667_r8 * c/abs_z
            ELSE
               res = &
                     -0.25000_r8 * z_c**5 &
                    + 0.50000_r8 * z_c**4 &
                    + 0.62500_r8 * z_c**3 &
                    - 1.66667_r8 * z_c**2 &
                    + 1_r8
            END IF
          END FUNCTION letkf_loc_gc

ENDMODULE localization
