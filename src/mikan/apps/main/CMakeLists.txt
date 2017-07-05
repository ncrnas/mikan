project (mk_main)

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
add_subdirectory(../../lib/miranda3 mr3as)
add_subdirectory(../../lib/pita_ddg pita)
add_subdirectory(../../lib/rnahybrid2_mfe rh2)
add_subdirectory(../../lib/targetminer1 tm1)
add_subdirectory(../../lib/targetscan5_cs ts5)
add_subdirectory(../../lib/two_step_svm tssvm)
add_subdirectory(../../lib/vienna_rna vrna)
add_subdirectory(../../lib/mikan mk)

# Add include directories
include_directories (${SEQAN_INCLUDE_DIRS}
                     ${MR3_HEADERS}
                     ${PITA_HEADERS}
                     ${RH2_HEADERS}
                     ${TM1_HEADERS}
                     ${TS5_HEADERS}
                     ${TSSVM_HEADERS}
                     ${VR16_HEADERS}
                     ${MIKAN_HEADERS})

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})

# Add executable applications
add_executable (mikan mikan.cpp)

# Add libraries
if(NOT TARGET mr3as)
  add_library(mr3as ${MR3_SOURCE})
endif()
if(NOT TARGET pita)
  add_library(pita ${PITA_SOURCE})
endif()
if(NOT TARGET rh2)
  add_library(rh2 ${RH2_SOURCE})
endif()
if(NOT TARGET tm1)
  add_library(tm1 ${TM1_SOURCE})
endif()
if(NOT TARGET ts5)
  add_library(ts5 ${TS5_SOURCE})
endif()
if(NOT TARGET tssvm)
  add_library(tssvm ${TSSVM_SOURCE})
endif()
if(NOT TARGET vrna)
  add_library(vrna ${VR16_SOURCE})
endif()
if(NOT TARGET mk)
  add_library(mk ${MIKAN_SOURCE})
endif()

# Add dependencies found by find_package (SeqAn)
if (MSVC)
  set (LIBM_NAME "")
else ()
  set (LIBM_NAME "m")
endif ()

target_link_libraries (mikan mr3as pita rh2 tm1 ts5 tssvm vrna mk ${SEQAN_LIBRARIES} ${LIBM_NAME})

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

include(${CMAKE_SOURCE_DIR}/test/lib/googletest/cmake/internal_utils.cmake)
config_compiler_and_linker()  # Defined in internal_utils.cmake.

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Install mk-miranda in ${PREFIX}/bin directory
install (TARGETS mikan
         DESTINATION bin)

install (FILES ${CMAKE_SOURCE_DIR}/../LICENSE
        ${CMAKE_SOURCE_DIR}/../README.md
        DESTINATION share/doc/mikan)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

set (CPACK_PACKAGE_NAME "mikan")
set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "mikan: ensemble approach to predict miRNA targes")
set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Takaya Saito <takaya.saito@outlook.com>")
set (CPACK_PACKAGE_VENDOR "Takaya Saito <takaya.saito@outlook.com>")
set (CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/../LICENSE)
set (CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/../README.md)
set (CPACK_PACKAGE_VERSION "${MIKAN_VERSION}")
set (CPACK_PACKAGE_VERSION_MAJOR "${MIKAN_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${MIKAN_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${MIKAN_VERSION_PATCH}")

include (CPack)
