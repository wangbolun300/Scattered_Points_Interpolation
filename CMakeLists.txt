cmake_minimum_required(VERSION 3.1)
project(Sparse_Interp)



################################################################################

# project-options


# libigl
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            OFF)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)

################################################################################

### Configuration
set(SPARSE_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(SPARSE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sparse_interp)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)
include(Sparse_InterpUtils)

#add_library(${PROJECT_NAME}
#  src/main.cpp
#
#)

add_subdirectory(sparse_interp)
add_library(${PROJECT_NAME}  
${SPARSE_SOURCES}
)
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core igl::triangle)
############################
# B spline library
include(B_spline/SOURCE.cmake)
prepend_current_path(SPLINE_SOURCES)
add_library(SPLINE_LIB ${SPLINE_SOURCES})
SI_copy_headers(${SPLINE_SOURCES})
SI_set_source_group(${SPLINE_SOURCES})
target_sources(SPLINE_LIB PRIVATE ${SPLINE_SOURCES})
target_include_directories(SPLINE_LIB PUBLIC B_spline/)


add_executable(${PROJECT_NAME}_bin
    src/main.cpp
  )
  # target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_bin PUBLIC SPLINE_LIB igl::core igl::opengl_glfw)
