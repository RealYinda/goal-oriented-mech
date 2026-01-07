cmake_minimum_required(VERSION 3.3)

project(ThermalSolver CXX C Fortran)

set(NDIM 3)
set(CREL opt)

set(libfind_ROOT "/opt/jaumin1.10")
set(JAUMIN_ROOT  "/opt/jaumin1.10")
set(JPSOL_ROOT   "/opt/jaumin1.10")
list(INSERT CMAKE_MODULE_PATH 0 ${libfind_ROOT}/share/cmake/Modules/scns)
list(APPEND CMAKE_MODULE_PATH   ${JAUMIN_ROOT}/share/cmake/Modules/scns)
list(APPEND CMAKE_MODULE_PATH   ${JPSOL_ROOT}/share/cmake/Modules/scns)

set(JPSOL_BACKEND JAUMIN)
set(JAUMIN_NDIM   ${NDIM})
set(JAUMIN_CREL   ${CREL})
set(JPSOL_NDIM    ${NDIM})
set(JPSOL_CREL    ${CREL})


# 设置 CGAL 的安装路径
set(CGAL_DIR /home/Dev/CGAL/cgal-releases-CGAL-4.13)
find_package(CGAL REQUIRED)
include(LibFindMacros)
libfind_check_external_packages(deps NAMES JPSOL ${JPSOL_BACKEND} mpi openmp REQUIRED)
include_directories(${deps_INCLUDE_DIRS})
add_definitions(${deps_DEFINITIONS})
# 指定Eigen路径
set(EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3")
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(control element fem interface integrator material shapefunction ./ )
file(GLOB_RECURSE src_files ./*.C ./*.cpp ./*.c ./*.f ./*.F90
                            ./fortran/${NDIM}d/*.f
                            ./fortran/${NDIM}d/*.F90)

add_executable(main${NDIM}d ${src_files})
target_link_libraries(main${NDIM}d ${JPSOL_LIBRARIES} ${JAUMIN_LIBRARIES} CGAL::CGAL)
