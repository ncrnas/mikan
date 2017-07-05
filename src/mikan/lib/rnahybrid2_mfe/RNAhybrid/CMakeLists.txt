cmake_minimum_required(VERSION 2.8.2)

set(HB_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/hybrid_core.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/hybrid_energy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/hybrid_table.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/hybrid_backtrace.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/hybrid_prettyprint.cpp
        PARENT_SCOPE
        )
set(HB_HEADERS
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        PARENT_SCOPE
        )