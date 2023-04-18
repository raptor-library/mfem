#[=======================================================================[.rst:
FindNAVIER
-------

Finds the Navier miniapp library.

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``NAVIER_FOUND``
  True if the system has the NAVIER library.
``NAVIER_INCLUDE_DIRS``
  Include directories needed to use NAVIER.
``NAVIER_LIBRARIES``
  Libraries needed to link to NAVIER.
#]=======================================================================]

find_package(NAVIER QUIET NO_MODULE)
if(NAVIER_FOUND)
    message("NAVIER ${NAVIER_FIND_VERSION} found.")
    if(NAVIER_FIND_COMPONENTS)
        message("NAVIER components found:")
        message("${NAVIER_FIND_COMPONENTS}")
    endif()
    return()
endif()

set(NAVIER_FOUND TRUE)

## Find headers and libraries
execute_process(COMMAND bash -c "spack location -i mfem" OUTPUT_VARIABLE MFEM_LOC OUTPUT_STRIP_TRAILING_WHITESPACE)
find_library(NAVIER_LIBRARIES navier PATHS "${MFEM_LOC}/share/mfem/miniapps/navier")
find_path(NAVIER_INCLUDE_DIRS navier_solver.hpp PATHS "${MFEM_LOC}/share/mfem/miniapps/navier/lib")
if(NOT NAVIER_LIBRARIES OR NOT NAVIER_INCLUDE_DIRS)
    set(NAVIER_FOUND FALSE)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NAVIER DEFAULT_MSG
                                  NAVIER_FOUND
                                  NAVIER_LIBRARIES
                                  NAVIER_INCLUDE_DIRS)

if (NAVIER_FOUND)
    add_library(NAVIER::NAVIER UNKNOWN IMPORTED)
    set_target_properties(NAVIER::NAVIER PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${NAVIER_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${NAVIER_LIBRARIES}
        )
endif()

mark_as_advanced(NAVIER_LIBRARIES NAVIER_INCLUDE_DIRS NAVIER_FOUND)
