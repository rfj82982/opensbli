cmake_minimum_required(VERSION 3.18)
project(APP C CXX)

# Sets the installation location (if it is not set), something which usually we
#Â won't do. This is used when we do `cmake --install .`
if(NOT APP_INSTALL_DIR)
  set(APP_INSTALL_DIR "$ENV{HOME}/OPS-APPS")
  message(
    WARNING
      "We will install the applications to ${APP_INSTALL_DIR}! Please use -DAPP_INSTALL_DIR to specify!"
  )
endif()

# if APPs are compiled independently then setup compilers, dependencies and find
# the OPS libraries on the system
if(${CMAKE_PROJECT_NAME} STREQUAL APP)
  # if show the compiling process in detail
  option(OPS_VERBOSE_WARNING "Turn on verbose warning messages" OFF)
  option(OPS_TEST "Turn on tests for Apps" OFF)
  if(NOT OPS_VERBOSE_WARNING)
    message(
      "We show concise compiling information by defautl! Use -DOPS_VERBOSE_WARNING=ON to switch on."
    )
  endif()
  set(CMAKE_VERBOSE_MAKEFILE ${OPS_VERBOSE_WARNING})
  set(OPS_APP_SRC ${CMAKE_SOURCE_DIR})
  # Use the Release mode by default
  if(NOT CMAKE_BUILD_TYPE)
    message(
      "We choose the Release mode! Please uses -DCMAKE_BUILD_TYPE to change.")
    set(CMAKE_BUILD_TYPE Release)
  endif()
  if(NOT LIBTRID_PATH)
    message(
      WARNING
        "The install path of the scalar tridiagonal solver library (LIBTRID_PATH) is not specified and the apps which need it will not compiled!"
    )
  endif()

  # Configure Compilers C
  set(CMAKE_C_STANDARD 99)
  if(${CMAKE_C_COMPILER_ID} STREQUAL GNU)
    set(CMAKE_C_FLAGS "-fPIC -Wall")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
  elseif(${CMAKE_C_COMPILER_ID} STREQUAL Clang)
    set(CMAKE_C_FLAGS "-fPIC -DUNIX -Wall")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_C_FLAGS_DEBUG "-g -O0")
  elseif(${CMAKE_C_COMPILER_ID} STREQUAL MSVC)
    set(CMAKE_C_FLAGS "-W1")
    set(CMAKE_C_FLAGS_RELEASE "/O2 /Ot ")
    set(CMAKE_C_FLAGS_DEBUG "/Z7 /DEBUG")
  elseif(${CMAKE_C_COMPILER_ID} STREQUAL Intel)
    set(CMAKE_C_FLAGS "-xHost")
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
  else()
    message(
      WARNING
        "There is no default flag for the chosen C compiler! Please use -DCFLAG=XXX to supply!"
    )
  endif()
  # Configure Compilers C++
  set(CMAKE_CXX_STANDARD 11)
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    set(CMAKE_CXX_FLAGS "-fPIC -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    set(CMAKE_CXX_FLAGS "-fPIC -DUNIX -Wall")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    set(CMAKE_CXX_FLAGS "-W1")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /Ot ")
    set(CMAKE_CXX_FLAGS_DEBUG "/Z7 /DEBUG")
  elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Intel)
    set(CMAKE_CXX_FLAGS "-xHost")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_DEBUG " -O0 -g -DOPS_DEBUG ")
  else()
    message(
      WARNING
        "There is no default flag for the chosen C++ compiler! Please use -DCXXFLAG=XXX to supply!"
    )
  endif()

  if(CFLAG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAG}")
  else()
    message(
      "Extra C compiler options by using CFLAG, e.g., -DCFLAG=\"-ffloat-store\" for GNU!"
    )
  endif()
  if(CXXFLAG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAG}")
  else()
    message(
      "Extra C++ compiler options by using CXXFLAG, e.g., -DCXXFLAG=\"-ffloat-store\" for GNU!"
    )
  endif()

  # Trying to find the OPS library
  set(SEARCH_PATHS "${OPS_INSTALL_DIR}" "/usr/local" "/usr")
  # Tyring to find the include directory
  set(OPS_HEAD
      ops_checkpointing.h
      ops_hdf5.h
      ops_lib_core.h
      ops_opencl_reduction.h
      ops_tridiag.h
      ops_cuda_reduction.h
      ops_instance.h
      ops_device_rt_support.h
      ops_util.h
      ops_device_rt_support.h
      ops_internal1.h
      ops_macros.h
      ops_seq.h
      queue.h
      ops_exceptions.h
      ops_internal2.h
      ops_mpi_core.h
      ops_seq_v2.h)
  set(OPS_INCLUDE_DIR "")
  foreach(head ${OPS_HEAD})
    unset(head_inc CACHE)
    find_path(
      head_inc
      NAMES ${head}
      PATHS ${SEARCH_PATHS}
      PATH_SUFFIXES include)
    if(head_inc)
      list(APPEND OPS_INCLUDE_DIR ${head_inc})
    else()
      message(FATAL_ERROR "We cannot find ${head}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES OPS_INCLUDE_DIR)
  include_directories("${OPS_INCLUDE_DIR}")
  # Trying to find the library directory
  if(LIBTRID_PATH)
    set(OPS_LIB
        ops_cuda
        ops_hdf5_seq
        ops_mpi_cuda
        ops_opencl
        ops_hdf5_mpi
        ops_mpi
        ops_mpi_opencl
        ops_seq
        ops_trid_seq
        ops_trid_cuda
        ops_trid_mpi
        ops_trid_mpi_cuda)
  else()
    set(OPS_LIB
        ops_cuda
        ops_hdf5_seq
        ops_mpi_cuda
        ops_opencl
        ops_hdf5_mpi
        ops_mpi
        ops_mpi_opencl
        ops_seq)
  endif()
  foreach(library ${OPS_LIB})
    find_library(
      ${library}_FOUND
      NAMES ${library}
      PATHS ${SEARCH_PATHS}
      PATH_SUFFIXES lib)
    if(${library}_FOUND)
      unset(dir CACHE)
      get_filename_component(dir ${${library}_FOUND} DIRECTORY)
      list(APPEND OPS_LIBRARY_DIR ${dir})
    else()
      message("We cannot find ${library}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES OPS_LIBRARY_DIR)
  link_directories("${OPS_LIBRARY_DIR}")
  # set the OPS_INSTALL_PATH for testing
  set(OPS_INSTALL_PATH "${OPS_INSTALL_DIR}")

  # Try to find the required dependency
  find_package(MPI QUIET)
  find_package(HDF5 COMPONENTS C HL QUIET)
  find_package(CUDAToolkit QUIET)
  find_package(OpenACC QUIET)
  find_package(OpenCL QUIET)
  find_package(OpenMP QUIET)
  find_package(Python3 REQUIRED)
  if(NOT Python3_FOUND)
    message(
      FATAL_ERROR
        "We cannot find Python3 and the Python translator needs Python3! Please use -DPython3_EXECUTABLE to specify the path."
    )
  else()
    find_path(
      opsc
      NAMES "ops.py"
      PATHS ${SEARCH_PATHS}
      PATH_SUFFIXES bin/ops_translator/c)
    if(NOT opsc)
      message(FATAL_ERROR "We cannot find the python C/C++ translator!")
    endif()
  endif()
  set(OPS_C_TRANSLATOR "${opsc}/ops.py")
  # Configure the "include" dir for compiling
  if(NOT HDF5_FOUND)
    message(
      WARNING
        "We cannot find the parallel HDF5 library. The HDF5 IO routines won't work! Please use -DHDF5_ROOT to specify the path!"
    )
  endif()
  if(NOT MPI_FOUND)
    message(
      WARNING "We cannot find a MPI environment! The MPI codes won't compile!")
  endif()

  if(LIBTRID_PATH)
    include_directories(${LIBTRID_PATH}/include)
  endif()

  if(CUDAToolkit_FOUND)
    set(CMAKE_CUDA_COMPILER ${CUDAToolkit_NVCC_EXECUTABLE})
    if(GPU_ARCH)
      set(CMAKE_CUDA_ARCHITECTURES
          70
          CACHE STRING "CUDA architectures")
    else()
      message(WARNING "Please the GPU architecture using -DGPU_ARCH=XXX!")
    endif()
    enable_language(CUDA)
  else()
    message(WARNING "We cannot find a CUDA kit! The CUDA codes won't compile!")
  endif()

  if(NOT OpenCL_FOUND)
    message(
      WARNING "We cannot find a OpenCL kit! The OpenCL codes won't compile!")
  endif()

  if(USE_NCCL)
    find_package(NCCL)
    if(NOT NCCL_FOUND)
      message(
        WARNING
          "We cannot find NCCL! Tridiagonal library wrapper will be compiled with NCCL"
      )
    endif()
  endif()

  set(HDF5_SEQ ${HDF5_FOUND})
  if(ops_cuda_FOUND)
    set(CUDA ${CUDAToolkit_FOUND})
  endif()
  if(ops_opencl_FOUND)
    set(OPENCL ${OpenCL_FOUND})
  endif()
  if(${MPI_FOUND})
    set(MPI TRUE)
    set(HDF5_MPI ${HDF5_FOUND})
    if(ops_cuda_FOUND)
      set(CUDA_MPI ${CUDAToolkit_FOUND})
    endif()
    if(ops_opencl_FOUND)
      set(OPENCL_MPI ${OpenCL_FOUND})
    endif()
  endif()

endif()


# Prepare the macro for compiling apps
# Name: App name
# Odd: Key words for source files are included in other source files and need no explicit compilation
# Others: Key words for source codes need no code generation using ops.py
# Extra: Key words for ource codes that cannot be treated by the macro
# Trid: Whether the tridiagonal library is required for this application
macro(BUILD_OPS_C_SAMPLE Name Odd Others Extra Trid)
  # Copy all source and head files into tmp, e.g. build/tmp. If tmp doesn't
  # exist then CMake will create it
  file(GLOB CPP "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
  list(FILTER CPP EXCLUDE REGEX "ops.cpp")
  file(COPY ${CPP} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tmp/")
  file(GLOB HEADS "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  file(COPY ${HEADS} DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tmp/")
  set(TMP_SOURCE_DIR "${CMAKE_CURRENT_BINARY_DIR}/tmp")
  # Get the files to be compileds:
  # DEV: non-optimised codes (full directory and filename), which I think are going to be translated
  # OTHERS: source codes which don't need code generation (full directory and filename)
  # OPS: codes generated by ops.py (full directory and filename)
  file(GLOB DEV "${TMP_SOURCE_DIR}/*.cpp")
  list(FILTER DEV EXCLUDE REGEX ${Odd})
  foreach(oth ${Others})
    list(FILTER DEV EXCLUDE REGEX ${oth})
  endforeach()
  foreach(ext ${Extra})
    list(FILTER DEV EXCLUDE REGEX ${ext})
  endforeach()
  list(GET DEV 0 Kernel)
  get_filename_component(KernerName ${Kernel} NAME_WE)
  # execute process to translate C and CPP source, e.g. $OPS_TRANSLATOR/c/ops.py
  execute_process(COMMAND ${OPS_C_TRANSLATOR} ${DEV} WORKING_DIRECTORY ${TMP_SOURCE_DIR})
  file(GLOB OPS "${TMP_SOURCE_DIR}/*ops*.cpp")
  list(FILTER OPS EXCLUDE REGEX ${Odd})
  file(GLOB OTHERS "${TMP_SOURCE_DIR}/*.cpp")
  foreach(OpsFile ${OPS})
    list(REMOVE_ITEM OTHERS ${OpsFile})
    string(REPLACE "_ops" "" NoOPS ${OpsFile})
    list(REMOVE_ITEM OTHERS ${NoOPS})
  endforeach()
  list(FILTER OTHERS EXCLUDE REGEX ${Odd})
  foreach(ext ${Extra})
    list(FILTER OTHERS EXCLUDE REGEX ${ext})
  endforeach()

  # Copy input parameters
  file(GLOB_RECURSE INPUT "${CMAKE_CURRENT_SOURCE_DIR}/*.in")


  # Compile sequential app
  file(COPY ${INPUT} DESTINATION "${TMP_SOURCE_DIR}")
  if(EXISTS "${TMP_SOURCE_DIR}/user_types.h")
    install(FILES "${TMP_SOURCE_DIR}/user_types.h"
            DESTINATION ${APP_INSTALL_DIR}/${Name})
  endif()
  install(FILES ${INPUT} DESTINATION ${APP_INSTALL_DIR}/${Name})
  add_executable(
    ${Name}_seq ${OPS} ${OTHERS}
                "${TMP_SOURCE_DIR}/MPI_OpenMP/${KernerName}_cpu_kernels.cpp")
  target_include_directories(${Name}_seq PRIVATE ${TMP_SOURCE_DIR})
  if(HDF5_SEQ)
    target_link_libraries(${Name}_seq ops_hdf5_seq hdf5::hdf5 hdf5::hdf5_hl
                          MPI::MPI_CXX)
  endif()
  if(${Trid} STREQUAL "YES")
    set(TRID Yes)
    target_link_libraries(${Name}_seq ops_trid_seq -L${LIBTRID_PATH}/lib
                          -ltridcpu)
  else()
    set(TRID No)
  endif()
  target_link_libraries(${Name}_seq ops_seq OpenMP::OpenMP_CXX)
  install(TARGETS ${Name}_seq DESTINATION ${APP_INSTALL_DIR}/${Name})

  # Compile CUDA app, if the CUDA OPS library was found
  if(CUDA)
    add_executable(
      ${Name}_cuda ${OPS} ${OTHERS}
                   "${TMP_SOURCE_DIR}/CUDA/${KernerName}_kernels.cu")
    target_include_directories(${Name}_cuda PRIVATE ${TMP_SOURCE_DIR})
    target_link_libraries(${Name}_cuda ops_cuda CUDA::cudart_static
                          OpenMP::OpenMP_CXX)
    target_compile_options(${Name}_cuda PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                                                --fmad=false>)
    if(HDF5_SEQ)
      target_link_libraries(${Name}_cuda ops_hdf5_seq hdf5::hdf5 hdf5::hdf5_hl
                            MPI::MPI_CXX)
    endif()
    if(TRID)
      target_link_libraries(
        ${Name}_cuda
        ops_trid_cuda
        -L${LIBTRID_PATH}/lib
        -ltridcuda
        CUDA::cublas_static
        CUDA::cublasLt_static
        CUDA::culibos)
    endif()
    install(TARGETS ${Name}_cuda DESTINATION ${APP_INSTALL_DIR}/${Name})
  endif()

  # Tridiagonal library does not currently support OpenCL
  # OpenCL needs the kenrel codes for first running
  install(DIRECTORY "${TMP_SOURCE_DIR}/OpenCL"
          DESTINATION ${APP_INSTALL_DIR}/${Name})
  if(OPENCL AND NOT TRID)
    add_executable(
      ${Name}_opencl
      ${OPS} ${OTHERS}
      "${TMP_SOURCE_DIR}/OpenCL/${KernerName}_opencl_kernels.cpp")
    target_include_directories(${Name}_opencl PRIVATE ${TMP_SOURCE_DIR}
                                                      ${OpenCL_INCLUDE_DIRS})
    target_link_libraries(${Name}_opencl ops_opencl ${OpenCL_LIBRARIES}
                          OpenMP::OpenMP_CXX)
    if(HDF5_SEQ)
      target_link_libraries(${Name}_opencl ops_hdf5_seq hdf5::hdf5
                            hdf5::hdf5_hl MPI::MPI_CXX)
    endif()
    install(TARGETS ${Name}_opencl DESTINATION ${APP_INSTALL_DIR}/${Name})
  endif()

  # Compile HIP app
  if(OPS_HIP)
    if(HIP AND NOT TRID)
      add_executable(
        ${Name}_hip ${OPS} ${OTHERS}
                    "${TMP_SOURCE_DIR}/HIP/${KernerName}_kernels.cpp")
      message("${TMP_SOURCE_DIR}/HIP/${KernerName}_hip_kernel.cpp")
      target_include_directories(${Name}_hip PRIVATE ${TMP_SOURCE_DIR})
      target_link_libraries(${Name}_hip ops_hip hip::device OpenMP::OpenMP_CXX)
    endif()
  endif()

  # Compile MPI app, including CUDA MPI if that is possible
  if(MPI)
    add_executable(
      ${Name}_mpi ${OPS} ${OTHERS}
                  "${TMP_SOURCE_DIR}/MPI_OpenMP/${KernerName}_cpu_kernels.cpp")
    target_include_directories(${Name}_mpi PRIVATE ${TMP_SOURCE_DIR})
    target_compile_definitions(${Name}_mpi PRIVATE "-DOPS_MPI")
    if(HDF5_MPI)
      target_link_libraries(${Name}_mpi PRIVATE ops_hdf5_mpi hdf5::hdf5
                                                hdf5::hdf5_hl MPI::MPI_CXX)
    endif()
    if(TRID)
      target_link_libraries(
        ${Name}_mpi PRIVATE ops_trid_mpi -L${LIBTRID_PATH}/lib -ltridcpu_mpi
                            MPI::MPI_CXX)
    endif()
    target_link_libraries(${Name}_mpi PRIVATE ops_mpi MPI::MPI_CXX
                                              OpenMP::OpenMP_CXX)
    install(TARGETS ${Name}_mpi DESTINATION ${APP_INSTALL_DIR}/${Name})
    # Compile CUDA MPI
    if(CUDA_MPI)
      add_executable(
        ${Name}_mpi_cuda ${OPS} ${OTHERS}
                         "${TMP_SOURCE_DIR}/CUDA/${KernerName}_kernels.cu")
      target_include_directories(${Name}_mpi_cuda PRIVATE ${TMP_SOURCE_DIR})
      target_compile_definitions(${Name}_mpi_cuda PRIVATE "-DOPS_MPI")
      target_compile_options(${Name}_mpi_cuda
                             PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: --fmad=false >)
      if(HDF5_MPI)
        target_link_libraries(
          ${Name}_mpi_cuda PRIVATE ops_hdf5_mpi hdf5::hdf5 hdf5::hdf5_hl
                                   MPI::MPI_CXX)
      endif()
      if(TRID)
        target_link_libraries(
          ${Name}_mpi_cuda PRIVATE ops_trid_mpi_cuda -L${LIBTRID_PATH}/lib
                                   -ltridcuda_mpi MPI::MPI_CXX)
        if(NCCL_FOUND)
          target_link_libraries(${Name}_mpi_cuda PRIVATE ${NCCL_LIBRARY})
        endif()
      endif()
      target_link_libraries(
        ${Name}_mpi_cuda PRIVATE ops_mpi_cuda CUDA::cudart_static MPI::MPI_CXX
                                 OpenMP::OpenMP_CXX)
      install(TARGETS ${Name}_mpi_cuda DESTINATION ${APP_INSTALL_DIR}/${Name})
    endif()

    # Compile OpenCL MPI app
    if(OPENCL_MPI AND NOT TRID)
      add_executable(
        ${Name}_mpi_opencl
        ${OPS} ${OTHERS}
        "${TMP_SOURCE_DIR}/OpenCL/${KernerName}_opencl_kernels.cpp")
      target_include_directories(
        ${Name}_mpi_opencl PRIVATE ${TMP_SOURCE_DIR} ${OpenCL_INCLUDE_DIRS})
      target_compile_definitions(${Name}_mpi_opencl PRIVATE "-DOPS_MPI")
      if(HDF5_MPI)
        target_link_libraries(
          ${Name}_mpi_opencl PRIVATE ops_hdf5_mpi hdf5::hdf5 hdf5::hdf5_hl
                                     MPI::MPI_CXX)
      endif()
      target_link_libraries(
        ${Name}_mpi_opencl PRIVATE ops_mpi_opencl ${OpenCL_LIBRARIES}
                                   MPI::MPI_CXX OpenMP::OpenMP_CXX)
      install(TARGETS ${Name}_mpi_opencl DESTINATION ${APP_INSTALL_DIR}/${Name})
    endif()
  endif()
endmacro()

build_ops_c_sample(OpenSBLI "NONE" "NONE" "NONE" "NO" "NO")