cmake_minimum_required(VERSION 2.8.2)

set(TSSVM_SVM_INIT_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_alpha.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u00.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u01.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u02.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u03.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u04.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u05.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u06.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u07.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u08.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u09.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u10.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u11.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u12.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u13.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u14.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tssvm_site_svm_sv_u15.cpp
        PARENT_SCOPE
        )
set(TSSVM_SVM_INIT_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )