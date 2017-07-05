cmake_minimum_required (VERSION 2.8.2)

set(MIKAN_SOURCE
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_memory.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_option.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_sequence.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/mk_core.cpp
   PARENT_SCOPE
)
set(MIKAN_HEADERS
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   PARENT_SCOPE
)
