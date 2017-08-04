cmake_minimum_required(VERSION 2.8.2)

set(TSSVM_SVM_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_feature.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_mrna_feature.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_mrna_svm.cpp
        PARENT_SCOPE
        )
set(TSSVM_SVM_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )