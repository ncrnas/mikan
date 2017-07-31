cmake_minimum_required(VERSION 2.8.2)

# Add sub directories.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/RNAhybrid)

set(RH2_SOURCE
        ${HB_SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/rh2_site_filter.cpp
        PARENT_SCOPE
        )
set(RH2_HEADERS
        ${HB_HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )
