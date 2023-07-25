### .txt

###########################################################################
### File Name: CMakeLists.txt
### Created on: 2023-03-22 20:09:03
### Author: YunhaoFu
### Last Modified: 2023-07-25 09:25:31
### Description: works for CMake or txt only
###########################################################################

cmake_minimum_required (VERSION 3.2)
enable_language(Fortran)

project (FYH-LETKF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
  "Build type, options are: Debug Release RelWithDebInfo" FORCE)
endif()

# set(CMAKE_Fortran_COMPILER,"/opt/intel/oneapi/compiler/2022.1.0/linux/bin/intel64/ifort")
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# MPI
#find_package(MPI REQUIRED)

# NetCDF
set(NETCDF_F90 "YES")
find_package(NetCDF REQUIRED)

# LAPACK / BLAS
find_package(LAPACK REQUIRED)

set( letkf_lib_srcs
    precision.F90
    model_state.F90
    obs_state.F90
    letkf_main.F90
    letkf_core.F90
    localization.F90
    kdtree.F90
)

add_library( letkf STATIC ${letkf_lib_srcs} )

# LAPACK
target_link_libraries(letkf  PRIVATE ${LAPACK_LIBRARIES})
#NetCDF
include_directories(${NETCDF_INCLUDE_DIRS})
target_link_libraries(letkf  PRIVATE  ${NETCDF_LIBRARIES})

add_executable(letkfdriver letkf.F90)

set_target_properties(letkfdriver PROPERTIES LINKER_LANGUAGE Fortran)
target_link_libraries(letkfdriver letkf )
