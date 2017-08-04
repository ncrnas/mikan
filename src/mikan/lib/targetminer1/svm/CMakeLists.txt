cmake_minimum_required(VERSION 2.8.2)

set(TM1_SVM_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_mrna_feature.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_mrna_svm.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tm1_site_feature.cpp
        PARENT_SCOPE
        )
set(TM1_SVM_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )