    !--------------------------------------------------!
    !                   parameter.h                    !
    INTEGER, PARAMETER            ::   m_lat_rever= 0  ! 1 -> reverse; 0 -> no reverse
    INTEGER, PARAMETER            ::   maxsnl     = -5
    INTEGER, PARAMETER            ::   nl_soil    = 10
    INTEGER, PARAMETER            ::   maxpatch   = 25
    INTEGER, PARAMETER            ::   da_cnt     = 8
    INTEGER          , PARAMETER  ::   invar_int_count = 2
    INTEGER          , PARAMETER  ::   invar_int_idx(invar_int_count)&
                                = (/ 3 ,&             ! Land water type
                                     4/)              ! Land cover type of classification
    INTEGER          , PARAMETER  ::   var_count = 8
    INTEGER          , PARAMETER  ::   var_sidx(var_count)&
                                = (/ 31 ,&           !  Soil temperature
                                     46 ,&           !  Liquid water in layers
                                     61 ,&           !  Ice lens in layers
                                     76 ,&           !  Ground surface temperature
                                     77 ,&           !  Sunlit leaf temperature
                                     78 ,&           !  Shaded leaf temperature
                                     81 ,&           !  Snow cover, water equivalent
                                     113/)           !  Radiative temperature of surface
    INTEGER          , PARAMETER  ::   var_length(var_count)&
                                = (/ 15 ,&
                                     15 ,&
                                     15 ,&
                                     1  ,&
                                     1  ,&
                                     1  ,&
                                     1  ,&
                                     1  /)
    INTEGER          , PARAMETER  ::   var_write(var_count)&  ! 0 -> only read; 1 -> read and write
                                = (/ 1 ,&
                                     1 ,&
                                     1 ,&
                                     1 ,&
                                     1 ,&
                                     1 ,&
                                     1 ,&
                                     0 /)
