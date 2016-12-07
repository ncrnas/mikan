project (test_mk_rh2)

cmake_minimum_required (VERSION 2.8.2)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search dependencies
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

add_subdirectory(${CMAKE_SOURCE_DIR}/mikan/lib/rnahybrid2_mfe rh2)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../lib/mikan mikantest)

# Add include directories
include_directories (${SEQAN_INCLUDE_DIRS} ${gtest_SOURCE_DIR}/include
        ${gtest_SOURCE_DIR} ${RH2_HEADERS} ${MKTEST_HEADERS})

# Add definitions set by find_package (SeqAn)
add_definitions (${SEQAN_DEFINITIONS})
if (MSVC)
    add_definitions (/D TEST_DATA_PATH=${CMAKE_CURRENT_LIST_DIR}/../data/)
else()
    add_definitions (-DTEST_DATA_PATH=${CMAKE_CURRENT_LIST_DIR}/../data/)
endif()

# Add executable applications
#add_executable (gtest_mk_pita seed_01_nmer.cpp seed_02_gu.cpp seed_03_mm.cpp seed_04_mmgu.cpp seed_05_bt.cpp
#        site_01_nmer_1.cpp site_01_nmer_2.cpp site_02_gu_1.cpp site_02_gu_2.cpp site_03_mm_1.cpp site_03_mm_2.cpp
#        site_04_mmgu_1.cpp site_04_mmgu_2.cpp site_05_bt_1.cpp)
add_executable (gtest_mk_rh2 seed_01_nmer.cpp seed_02_gu.cpp site_01_nmer_1.cpp site_01_nmer_2.cpp
        site_02_gu_1.cpp site_02_gu_2.cpp )

# Add libraries
if(NOT TARGET mikantest)
    add_library(mikantest ${MKTEST_SOURCE})
endif()

# Add dependencies found by find_package (SeqAn)
target_link_libraries (gtest_mk_rh2 gtest gtest_main)
target_link_libraries (gtest_mk_rh2 ${SEQAN_LIBRARIES} rh2 mikantest)

# Add CXX flags found by find_package (SeqAn)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Define helper functions and macros used by Google Test.
include(${CMAKE_CURRENT_SOURCE_DIR}/../lib/googletest/cmake/internal_utils.cmake)
config_compiler_and_linker()  # Defined in internal_utils.cmake.

# Add test
add_test (NAME test_mk_rh2 COMMAND gtest_mk_rh2)
