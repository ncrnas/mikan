project(mk_mikanconf)

cmake_minimum_required(VERSION 2.8.2)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies
set(SEQAN_FIND_DEPENDENCIES NONE)
find_package(SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add lib directory
add_subdirectory(../../lib/mikan_ens mkens)

# Add include directory
include_directories(${SEQAN_INCLUDE_DIRS} ${MKE_HEADERS})

# Add definitions set by find_package (SeqAn)
add_definitions(${SEQAN_DEFINITIONS})

# Add executable application
add_executable(mkconf ../../lib/mikan_ens/mke_config.cpp mikan_config.cpp)

# Add dependencies found by find_package (SeqAn)
target_link_libraries(mkconf ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

config_compiler_and_linker()  # Defined in internal_utils.cmake.

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Install mk-miranda in ${PREFIX}/bin directory
install(TARGETS mkconf
        DESTINATION bin)

install(FILES ${CMAKE_SOURCE_DIR}/../LICENSE
        ${CMAKE_SOURCE_DIR}/../README.md
        DESTINATION share/doc/mikan)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

set(CPACK_PACKAGE_NAME "mkconf")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "mkconf: print mikan config template to stdout")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Takaya Saito <takaya.saito@outlook.com>")
set(CPACK_PACKAGE_VENDOR "Takaya Saito <takaya.saito@outlook.com>")
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_SOURCE_DIR}/../LICENSE)
set(CPACK_RESOURCE_FILE_README ${CMAKE_SOURCE_DIR}/../README.md)
set(CPACK_PACKAGE_VERSION "${MIKAN_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${MIKAN_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${MIKAN_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${MIKAN_VERSION_PATCH}")

include(CPack)
