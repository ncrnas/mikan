cmake_minimum_required (VERSION 2.8.2)

set(VR16_SOURCE
   ${CMAKE_CURRENT_SOURCE_DIR}/vr16_ddg_core.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/vr16_duplex.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/vr16_fold.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/vr16_part_func.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/vr16_energy.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/vr16_pair_mat.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/vr16_params.cpp
   PARENT_SCOPE
)
set(VR16_HEADERS
   ${CMAKE_CURRENT_SOURCE_DIR}/include
   PARENT_SCOPE
)