cmake_minimum_required(VERSION 2.8.2)

set(MIKAN_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_nd_array.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_sequence.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_input.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_seed_seq.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_site_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_rna_with_sites.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_site_filter.cpp
        PARENT_SCOPE
        )
set(MIKAN_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )
