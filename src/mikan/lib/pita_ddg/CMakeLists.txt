cmake_minimum_required(VERSION 2.8.2)

set(PITA_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/pita_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pita_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pita_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pita_site_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pita_site_filter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/pita_rna_score.cpp
        PARENT_SCOPE
        )
set(PITA_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )
