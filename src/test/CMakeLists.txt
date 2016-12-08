project (test)

cmake_minimum_required (VERSION 2.8.2)

# Add include paths
set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/seqan/extras/include
        ${CMAKE_SOURCE_DIR}/seqan/core/include
        ${CMAKE_SOURCE_DIR})

add_subdirectory (lib/googletest)
add_subdirectory (mk_common)
add_subdirectory (mk_miranda)
add_subdirectory (mk_pita)
add_subdirectory (mk_rh2)
add_subdirectory (mk_tm1)
add_subdirectory (mk_ts5)
