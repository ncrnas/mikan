cmake_minimum_required(VERSION 2.8.2)

# Add sub directory
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/svm)

set(TM1_SOURCE
        ${TM1_SVM_SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_site_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_site_filter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_rna_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_align.cpp
        PARENT_SCOPE
        )
set(TM1_HEADERS
        ${TM1_SVM_HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )

