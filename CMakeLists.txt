cmake_minimum_required(VERSION 2.8)

project(orbopt Fortran)

# HOWTO compile the sources from the source code directory:
#
# $ mkdir build && cd build
# $ ccmake ..
# $ make -j4
#
# done! you may want to set the compile type to Debug or Release in the ccmake GUI
#
# in order to set the fortran/c/c++ compiler use in the build directory:
# $ FC=ifort cmake ..
#
# to build with parallel MKL:
# ---------------------------
# Intel
# $ FC=ifort cmake -DENABLE_THREADED_MKL=ON..
# GNU
# $ FC=gfortran cmake -DENABLE_THREADED_MKL=ON..

option(ENABLE_AUTO_BLAS       "Enable CMake to autodetect BLAS"     ON)
option(ENABLE_AUTO_LAPACK     "Enable CMake to autodetect LAPACK"   ON)
option(ENABLE_PROFILE         "Compile with profiling option    "   OFF)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/math
   )

if(NOT CMAKE_Fortran_MODULE_DIRECTORY)
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/modules)
endif()

include_directories(
${CMAKE_Fortran_MODULE_DIRECTORY}
)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Release
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

include(ConfigCompilerFlags)

######################################################################
# Version information
######################################################################
set(ORBOPT_YEAR 2016)
set(ORBOPT_VERSION_MAJOR 0)
set(ORBOPT_VERSION_MINOR 2)
set(ORBOPT_VERSION_BUILD "")

if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
    find_package(Git)
    if(GIT_FOUND)
        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-list --max-count=1 HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_REVISION
            ERROR_QUIET
            )
        if(NOT ${GIT_REVISION} STREQUAL "")
            string(STRIP ${GIT_REVISION} GIT_REVISION)
        endif()

        execute_process(
            COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_BRANCH
            ERROR_QUIET
            )
        if(NOT ${GIT_BRANCH} STREQUAL "")
            string(STRIP ${GIT_BRANCH} GIT_BRANCH)
        endif()

    endif()
else(EXISTS ${PROJECT_SOURCE_DIR}/.git)
    set(GIT_REVISION "N/A")
    set(GIT_BRANCH   "N/A")
endif(EXISTS ${PROJECT_SOURCE_DIR}/.git)

set(ORBOPT_VERSION_BUILD "${GIT_REVISION} (${GIT_BRANCH})")

set(ORBOPT_VERSION "${ORBOPT_VERSION_MAJOR}.${ORBOPT_VERSION_MINOR}")
if(ORBOPT_VERSION_BUILD)
  set(ORBOPT_GIT_VERSION "${ORBOPT_VERSION_BUILD}")
else(ORBOPT_VERSION_BUILD)
  set(ORBOPT_GIT_VERSION "")
endif(ORBOPT_VERSION_BUILD)
set(ORBOPT_VERSION_STRING "ORBOPT - version: ${ORBOPT_VERSION}")
MESSAGE(STATUS "ORBOPT - version   :  ${ORBOPT_VERSION}")
MESSAGE(STATUS "ORBOPT - git commit:  ${GIT_REVISION} (branch: ${GIT_BRANCH})")

######################################################################
# Configure files
######################################################################

configure_file(${PROJECT_SOURCE_DIR}/src/orbopt_version.f90.in ${PROJECT_SOURCE_DIR}/src/orbopt_version.f90)
######################################################################
# Install
######################################################################

install(FILES ${PROJECT_SOURCE_DIR}/src/orbopt_version.f90
        DESTINATION include COMPONENT headers)

set(ORBOPT_SOURCE_FILES
   ${PROJECT_SOURCE_DIR}/src/hostname.F90
   ${PROJECT_SOURCE_DIR}/src/orbopt_header.F90
   ${PROJECT_SOURCE_DIR}/src/Closed_Shell.f90
   ${PROJECT_SOURCE_DIR}/src/ConDMRGLR.f90
   ${PROJECT_SOURCE_DIR}/src/Coupled_solver.f90
   ${PROJECT_SOURCE_DIR}/src/Coupling.f90
   ${PROJECT_SOURCE_DIR}/src/Coupling_terms.f90
   ${PROJECT_SOURCE_DIR}/src/Coupling_update.f90
   ${PROJECT_SOURCE_DIR}/src/Date_time.f90
   ${PROJECT_SOURCE_DIR}/src/Davidson.f90
   ${PROJECT_SOURCE_DIR}/src/Derivaties.f90
   ${PROJECT_SOURCE_DIR}/src/Derivaties_deltaR.f90
   ${PROJECT_SOURCE_DIR}/src/Energy_gen.f90
   ${PROJECT_SOURCE_DIR}/src/Environment.f90
   ${PROJECT_SOURCE_DIR}/src/Fock_gen.f90
   ${PROJECT_SOURCE_DIR}/src/Global_control.f90
   ${PROJECT_SOURCE_DIR}/src/Gmat_gen.f90
   ${PROJECT_SOURCE_DIR}/src/Grouptable.f90
   ${PROJECT_SOURCE_DIR}/src/HFenergy.f90
   ${PROJECT_SOURCE_DIR}/src/Initial_value.f90
   ${PROJECT_SOURCE_DIR}/src/Input.f90
   ${PROJECT_SOURCE_DIR}/src/Integrals.f90
   ${PROJECT_SOURCE_DIR}/src/Iteration.f90
   ${PROJECT_SOURCE_DIR}/src/Iteration_ext.f90
   ${PROJECT_SOURCE_DIR}/src/LagrangeRDMs.f90
   ${PROJECT_SOURCE_DIR}/src/MPSci_update.f90
   ${PROJECT_SOURCE_DIR}/src/Mathlib.f90
   ${PROJECT_SOURCE_DIR}/src/Matrix.f90
   ${PROJECT_SOURCE_DIR}/src/NonLinear.f90
   ${PROJECT_SOURCE_DIR}/src/Operator_update.f90
   ${PROJECT_SOURCE_DIR}/src/Operator_update_CPNR.f90
   ${PROJECT_SOURCE_DIR}/src/Operator_update_LR.f90
   ${PROJECT_SOURCE_DIR}/src/Orbital_Hessian.f90
   ${PROJECT_SOURCE_DIR}/src/PmultU.f90
   ${PROJECT_SOURCE_DIR}/src/PreDMRGLR.f90
   ${PROJECT_SOURCE_DIR}/src/PreOneSCF.f90
   ${PROJECT_SOURCE_DIR}/src/Precaii.f90
   ${PROJECT_SOURCE_DIR}/src/Precabb.f90
   ${PROJECT_SOURCE_DIR}/src/Print_info.f90
   ${PROJECT_SOURCE_DIR}/src/Redundancies.f90
   ${PROJECT_SOURCE_DIR}/src/Run_DMRG.f90
   ${PROJECT_SOURCE_DIR}/src/Scf_save.f90
   ${PROJECT_SOURCE_DIR}/src/Solver_Augmented_Hessian.f90
   ${PROJECT_SOURCE_DIR}/src/Solver_Hessian.f90
   ${PROJECT_SOURCE_DIR}/src/Solver_LR.f90
   ${PROJECT_SOURCE_DIR}/src/Solver_LR_partition.f90
   ${PROJECT_SOURCE_DIR}/src/Solver_WMK.f90
   ${PROJECT_SOURCE_DIR}/src/Solver_WMK2.f90
   ${PROJECT_SOURCE_DIR}/src/Transform_INT.f90
   ${PROJECT_SOURCE_DIR}/src/Transform_Nindex.f90
   ${PROJECT_SOURCE_DIR}/src/Y_gen.f90
   ${PROJECT_SOURCE_DIR}/src/orbopt_version.f90
)

set(ORBOPT_MAIN
   ${PROJECT_SOURCE_DIR}/src/MAIN.f90
)

#math detection from DIRAC...
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)

add_library(ORBOPT ${ORBOPT_SOURCE_FILES})

macro(compile_standalone _executable _source LIST_OF_EXECUTABLES)
    add_executable(
        ${_executable}
        ${_source}
    )
    set_property(TARGET ${_executable} PROPERTY LINKER_LANGUAGE Fortran)
    set(LIST_OF_EXECUTABLES
        ${LIST_OF_EXECUTABLES}
        ${_executable}
    )
endmacro()

# set standalone executables
compile_standalone(orb_opt "${ORBOPT_MAIN}" "${LIST_OF_EXECUTABLES}")

foreach(
    EXECUTABLE
    ${LIST_OF_EXECUTABLES}
    )
    if(${CMAKE_SYSTEM_NAME} STREQUAL "AIX")
        SET_TARGET_PROPERTIES(${EXECUTABLE} PROPERTIES LINK_FLAGS "-Wl,-bbigtoc")
    endif()
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        SET_TARGET_PROPERTIES(${EXECUTABLE} PROPERTIES LINK_FLAGS "-Wl,-E")
    endif()
    target_link_libraries(
        ${EXECUTABLE}
        ORBOPT
        ${LAPACK_LIBRARIES}
        ${BLAS_LIBRARIES}
        ${LINK_LIBS}
        )
endforeach()

# copy input files to binary location (to be consistent with old Makefile version)
message ("-- Copying input data for ORBOPT")
configure_file(${PROJECT_SOURCE_DIR}/input-files/INTEGRAL_ORDER_1_IRREPS.grp ${PROJECT_BINARY_DIR} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/input-files/INTEGRAL_ORDER_2_IRREPS.grp ${PROJECT_BINARY_DIR} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/input-files/INTEGRAL_ORDER_4_IRREPS.grp ${PROJECT_BINARY_DIR} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/input-files/INTEGRAL_ORDER_8_IRREPS.grp ${PROJECT_BINARY_DIR} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/input-files/MASORB.orb                  ${PROJECT_BINARY_DIR} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/input-files/infile.inp                  ${PROJECT_BINARY_DIR} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/rdmsave_su2.py                      ${PROJECT_BINARY_DIR} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/pickup.py                           ${PROJECT_BINARY_DIR} COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/src/pickup_eff.py                       ${PROJECT_BINARY_DIR} COPYONLY)

