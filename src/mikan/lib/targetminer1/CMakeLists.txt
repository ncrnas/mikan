cmake_minimum_required(VERSION 2.8.2)

set(TM1_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_site_feature.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_site_cluster.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_mrna_feature.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_mrna_svm.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_align.cpp
        PARENT_SCOPE
        )
set(TM1_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )

