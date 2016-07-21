cmake_minimum_required (VERSION 2.8.2)

# Add sub directories.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../vienna_rna vrna)

set(PITA_SOURCE
   ${VR16_SOURCE}
   ${CMAKE_CURRENT_SOURCE_DIR}/pita_core.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/pita_option.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/pita_seed_site.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/pita_score.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/pita_site_cluster.cpp
   PARENT_SCOPE
)
set(PITA_HEADERS
   ${VR16_HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   PARENT_SCOPE
)
