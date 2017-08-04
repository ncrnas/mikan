cmake_minimum_required(VERSION 2.8.2)

# Add sub directory
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mfe)

set(RH2_SOURCE
        ${RH2_MFE_SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_site_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_site_filter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_rna_score.cpp
        PARENT_SCOPE
        )
set(RH2_HEADERS
        ${RH2_MFE_HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )
