cmake_minimum_required(VERSION 3.12)

project(Advection C)

# ---- Configure compiler ----------------------------------------------------
if(CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release)
endif()

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  add_compile_options(-w -O3)
elseif (${CMAKE_BUILD_TYPE} MATCHES Debug)
  add_compile_options(
    -g
    -O0
    -Wall
    -Wextra
    -Wpedantic
    -Wno-newline-eof
    -Wno-missing-braces
    -Wno-char-subscripts
    -Wno-unused-variable
    -Wno-unused-function
    -Wno-unused-parameter
    )
else()
  message(WARNING "Unhandled build type: ${CMAKE_BUILD_TYPE}")
endif()

# ---- Variable documentation ------------------------------------------------
set(ISCD_COMMONS_LIBRARY "" CACHE STRING "Path to ISCD::Commons library")
set(ISCD_COMMONS_INCLUDE_DIR "" CACHE STRING "Path to ISCD::Commons include directory")

# ---- ISCD::Commons ---------------------------------------------------------
find_library(Commons
  NAMES
  Commons
  HINTS
  ${ISCD_COMMONS_LIBRARY}
  )

# ---- Math ------------------------------------------------------------------
find_library(M_LIB m)

# ---- Executable ------------------------------------------------------------
file(GLOB_RECURSE Advection_SRCS sources/*.c)
file(GLOB_RECURSE Advection_HEADERS sources/*.h)

add_executable(Advection ${Advection_SRCS} ${Advection_HEADERS})

target_include_directories(Advection
  PUBLIC
  ${ISCD_COMMONS_INCLUDE_DIR}
  )
target_link_libraries(Advection
  PUBLIC
  pthread
  ${M_LIB}
  ${Commons}
  )
