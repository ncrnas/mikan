# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /mikan/apps/mk/CMakeLists.txt
#
# CMakeLists.txt file for mk apps
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (mk_miranda)

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
add_subdirectory(../../lib/miranda3 mr3as)

# Add include directories
include_directories (${SEQAN_INCLUDE_DIRS}
                     ${MR3_HEADERS})

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})

# Add executable applications
add_executable (mk-miranda mkmiranda/miranda3.cpp)

# Add libraries
add_library(mr3as ${MR3_SOURCE})

# Add dependencies found by find_package (SeqAn)
if (MSVC)
  set (LIBM_NAME "")
else ()
  set (LIBM_NAME "m")
endif ()

target_link_libraries (mk-miranda mr3as ${SEQAN_LIBRARIES} ${LIBM_NAME})

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

include(${CMAKE_SOURCE_DIR}/test/lib/googletest/cmake/internal_utils.cmake)
config_compiler_and_linker()  # Defined in internal_utils.cmake.

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Install mk-miranda in ${PREFIX}/bin directory
install (TARGETS mk-miranda
         DESTINATION bin)

install (FILES ${CMAKE_SOURCE_DIR}/../LICENSE
        ${CMAKE_SOURCE_DIR}/../README.md
        DESTINATION share/doc/mikan/mkmiranda)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

set (CPACK_PACKAGE_NAME "mk_miranda")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "mk-miranda: mikan version of the miRanda algorithm")
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Takaya Saito <takaya.saito@outlook.com>")
set (CPACK_PACKAGE_VENDOR "Takaya Saito <takaya.saito@outlook.com>")
set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/../LICENSE)
set (CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/../README.md)
set (CPACK_PACKAGE_VERSION "${MIKAN_VERSION}")
set (CPACK_PACKAGE_VERSION_MAJOR "${MIKAN_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${MIKAN_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${MIKAN_VERSION_PATCH}")

include (CPack)
