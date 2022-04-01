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

# ---- ISCD::Commons ---------------------------------------------------------
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/third-party/ISCD/Commons
  ${PROJECT_BINARY_DIR}/third-party/ISCD/Commons)
target_include_directories(Commons
  INTERFACE
  ${PROJECT_SOURCE_DIR}/third-party/ISCD/Commons/sources)
add_library(ISCD::Commons ALIAS Commons)

# ---- Executable ------------------------------------------------------------
file(GLOB_RECURSE Advection_SRCS sources/*.c)
file(GLOB_RECURSE Advection_HEADERS sources/*.h)

add_executable(Advection ${Advection_SRCS} ${Advection_HEADERS})

target_link_libraries(Advection
  PRIVATE
  m
  pthread
  ISCD::Commons)
