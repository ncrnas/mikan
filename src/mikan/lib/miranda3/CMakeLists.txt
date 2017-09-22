cmake_minimum_required(VERSION 2.8.2)

set(MR3_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/mr3_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mr3_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mr3_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mr3_site_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mr3_align.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mr3_site_filter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mr3_rna_score.cpp
        PARENT_SCOPE
        )
set(MR3_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/dp
        PARENT_SCOPE
        )
