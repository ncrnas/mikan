cmake_minimum_required (VERSION 2.8.2)

# Add sub directories.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../vienna_rna vrna)

set(MR3_SOURCE
   ${VR16_SOURCE}
   ${CMAKE_CURRENT_SOURCE_DIR}/mr3_core.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/mr3_option.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/mr3_seed_site.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/mr3_score.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/mr3_align.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/mr3_site_cluster.cpp
   PARENT_SCOPE
)
set(MR3_HEADERS
   ${VR16_HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   PARENT_SCOPE
)
