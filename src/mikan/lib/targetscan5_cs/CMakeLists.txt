cmake_minimum_required (VERSION 2.8.2)

set(TS5_SOURCE
   ${CMAKE_CURRENT_SOURCE_DIR}/ts5_core.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/ts5_option.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/ts5_seed_site.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/ts5_feature.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/ts5_score.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/ts5_align.cpp
   PARENT_SCOPE
)
set(TS5_HEADERS
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   PARENT_SCOPE
)

