# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /mikan/apps/mk/CMakeLists.txt
#
# CMakeLists.txt file for mk apps
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (mikan_apps_mk)
message (STATUS "Configuring mikan/apps/mk")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add lib directories
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/rnahybrid2_mfe rh2)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/targetscan5_cs ts5)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/two_step_svm tssvm)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/pita_ddg pita)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/miranda3 mr3as)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/targetminer1 tm1)

# Add include directories
include_directories (${SEQAN_INCLUDE_DIRS} 
                     ${RH2_HEADERS}
                     ${TS5_HEADERS}
                     ${TSSVM_HEADERS}
                     ${PITA_HEADERS}
                     ${MR3_HEADERS}
                     ${TM1_HEADERS})

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})

# Add executable applications
add_executable (rnahybrid2_mfe rnahybrid2_mfe.cpp)
add_executable (targetscan5_cs targetscan5_cs.cpp)
add_executable (two_step_svm two_step_svm.cpp)
add_executable (pita_ddg pita_ddg.cpp)
add_executable (miranda3 miranda3.cpp)
add_executable (targetminer1 targetminer1.cpp)

# Add libraries
add_library(rh2 ${RH2_SOURCE})
add_library(ts5 ${TS5_SOURCE})
add_library(tssvm ${TSSVM_SOURCE})
add_library(pita ${PITA_SOURCE})
add_library(mr3as ${MR3_SOURCE})
add_library(tm1 ${TM1_SOURCE})

# Add dependencies found by find_package (SeqAn)
target_link_libraries (rnahybrid2_mfe ${SEQAN_LIBRARIES} rh2 m)
target_link_libraries (targetscan5_cs ${SEQAN_LIBRARIES} ts5 m)
target_link_libraries (two_step_svm ${SEQAN_LIBRARIES} tssvm m)
target_link_libraries (pita_ddg ${SEQAN_LIBRARIES} pita m)
target_link_libraries (miranda3 ${SEQAN_LIBRARIES} mr3as m)
target_link_libraries (targetminer1 ${SEQAN_LIBRARIES} tm1 m)

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (rnahybrid2_mfe)
  seqan_setup_install_vars (targetscan5_cs)
  seqan_setup_install_vars (two_step_svm)
  seqan_setup_install_vars (pita_ddg)
  seqan_setup_install_vars (targetminer1)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install mk programs in ${PREFIX}/bin directory
install (TARGETS rnahybrid2_mfe
         DESTINATION bin)
install (TARGETS targetscan5_cs
         DESTINATION bin)
install (TARGETS two_step_svm
         DESTINATION bin)
install (TARGETS pita_ddg
         DESTINATION bin)
install (TARGETS miranda3
         DESTINATION bin)
install (TARGETS targetminer1
         DESTINATION bin)
                  
# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/mk for SeqAn release builds
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (empt)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:mk")
  set (CPACK_PACKAGE_NAME "mk")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "empt")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (mk "mk")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:mk")

