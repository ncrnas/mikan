cmake_minimum_required(VERSION 2.8.2)

# Add sub directories
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/svm)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/svm_init)

set(TSSVM_SOURCE
        ${TSSVM_SVM_SOURCE}
        ${TSSVM_SVM_INIT_SOURCE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_seed_site.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_align.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_filter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_score.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_rna_score.cpp
        PARENT_SCOPE
        )
set(TSSVM_HEADERS
        ${TSSVM_SVM_HEADERS}
        ${TSSVM_SVM_INIT_HEADERS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )

