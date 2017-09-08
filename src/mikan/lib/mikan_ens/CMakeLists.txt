cmake_minimum_required(VERSION 2.8.2)

set(MKE_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/mke_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mke_config.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mke_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mke_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mke_site_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mke_rna_sites.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mke_rna_score.cpp
        PARENT_SCOPE
        )
set(MKE_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )

