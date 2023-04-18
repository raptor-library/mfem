cmake_minimum_required(VERSION 3.10)

# set the project name
project(spack-example-mfem LANGUAGES CXX C)

# Compiler defs
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic -Wall -Wextra")
include(CMakePrintHelpers)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# CUDA
include(CheckLanguage)
check_language(CUDA)
if(DEFINED CMAKE_CUDA_COMPILER AND NOT ${CMAKE_CUDA_COMPILER} MATCHES "NOTFOUND")
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    set(CMAKE_CUDA_ARCHITECTURES "70")
    set(CMAKE_CUDA_FLAGS "--expt-extended-lambda")
endif()


# CMake modules
message(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/config/cmake/modules)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(ZLIB REQUIRED)
find_package(SUITESPARSE REQUIRED)
find_package(METIS REQUIRED)
find_package(HYPRE REQUIRED)
find_package(MFEM REQUIRED)
find_package(NAVIER REQUIRED)
find_package(GTest REQUIRED)
find_package(CEED REQUIRED)
find_package(SUNDIALS REQUIRED)

#Some options used in executables and libraries
function(set_mfem_target_options target source_file)
target_compile_options(${target}
    PRIVATE ${LAPACK_LINKER_FLAGS}
    PRIVATE ${BLAS_LINKER_FLAGS}
    PRIVATE ${MPI_LINK_FLAGS}
    )

target_link_libraries(${target}
    PRIVATE NAVIER::NAVIER
    PRIVATE MFEM::MFEM
    PRIVATE SUITESPARSE::SUITESPARSE
    PRIVATE SUITESPARSE::UMFPACK
    PRIVATE SUITESPARSE::KLU
    PRIVATE SUNDIALS::nvecserial
    PRIVATE SUNDIALS::nvecparallel
    PRIVATE METIS::METIS
    PRIVATE HYPRE::HYPRE
    PRIVATE ${LAPACK_LIBRARIES}
    PRIVATE ${BLAS_LIBRARIES}
    PRIVATE ${MPI_C_LIBRARIES}
    PRIVATE ZLIB::ZLIB
    PRIVATE GTest::gtest_main
    PRIVATE CEED::CEED
    )


target_include_directories(${target}
    PRIVATE ${MPI_C_INCLUDE_DIRS}
    )

if(DEFINED CMAKE_CUDA_COMPILER AND NOT ${CMAKE_CUDA_COMPILER} MATCHES "NOTFOUND")
    set_source_files_properties(${source_file} PROPERTIES LANGUAGE CUDA)
    target_link_libraries(${target}
        PRIVATE cuda
        PRIVATE cusparse
        PRIVATE cublasLt #new target for cuda > 10.1
        )
endif()

endfunction()

#For executables using mfem
function(add_mfem_executable executable_name source_file)

#Add the executable
add_executable(${executable_name} ${source_file})

set_mfem_target_options(${executable_name} ${source_file})

install(TARGETS ${executable_name}
    DESTINATION ${PROJECT_SOURCE_DIR}/bin
    )

endfunction()

#For libraries using mfem
function(add_mfem_library library_name source_file)

# add the library
add_library(${library_name} ${source_file})

set_mfem_target_options(${library_name} ${source_file})

install(TARGETS ${library_name}
    EXPORT fluids-mfem
    DESTINATION ${PROJECT_SOURCE_DIR}/lib
    )

endfunction()

# Testing
enable_testing()
include(GoogleTest)

add_subdirectory(src)
add_subdirectory(test)
