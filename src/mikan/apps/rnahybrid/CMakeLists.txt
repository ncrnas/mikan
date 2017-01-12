project (mk_rnahybrid)

cmake_minimum_required (VERSION 2.8.2)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add lib directories
add_subdirectory(../../lib/rnahybrid2_mfe rh2)
add_subdirectory(../../lib/mikan mk)

# Add include directories
include_directories (${SEQAN_INCLUDE_DIRS}
                     ${RH2_HEADERS}
                     ${MIKAN_HEADERS})

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})

# Add executable applications
add_executable (mk-rnahybrid rnahybrid2_mfe.cpp)

# Add libraries
add_library(rh2 ${RH2_SOURCE})
add_library(mk ${MIKAN_SOURCE})

# Add dependencies found by find_package (SeqAn)
if (MSVC)
  set (LIBM_NAME "")
else ()
  set (LIBM_NAME "m")
endif ()

target_link_libraries (mk-rnahybrid rh2 mk ${SEQAN_LIBRARIES} ${LIBM_NAME})

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

include(${CMAKE_SOURCE_DIR}/test/lib/googletest/cmake/internal_utils.cmake)
config_compiler_and_linker()  # Defined in internal_utils.cmake.

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Install mk-rnahybrid in ${PREFIX}/bin directory
install (TARGETS mk-rnahybrid
         DESTINATION bin)

install (FILES ${CMAKE_SOURCE_DIR}/../LICENSE
        ${CMAKE_SOURCE_DIR}/../README.md
        DESTINATION share/doc/mikan/mk_rnahybrid)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

set (CPACK_PACKAGE_NAME "mk-rnahybrid")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "mk-rnahybrid: mikan version of the RNAhybrid algorithm")
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Takaya Saito <takaya.saito@outlook.com>")
set (CPACK_PACKAGE_VENDOR "Takaya Saito <takaya.saito@outlook.com>")
set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/../LICENSE)
set (CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/../README.md)
set (CPACK_PACKAGE_VERSION "${MIKAN_VERSION}")
set (CPACK_PACKAGE_VERSION_MAJOR "${MIKAN_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${MIKAN_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${MIKAN_VERSION_PATCH}")

include (CPack)
