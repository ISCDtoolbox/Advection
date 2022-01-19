cmake_minimum_required(VERSION 3.12)
project(Advection)

file(GLOB_RECURSE sources sources/*.c)
file(GLOB_RECURSE headers sources/*.h)

add_library(libAdvection ${sources})

find_library(Commons
  NAMES
  Commons
  HINTS
  ${ISCD_COMMONS_LIBRARY}
  )

find_library(M_LIB m)

target_include_directories(libAdvection PUBLIC ${ISCD_COMMONS_INCLUDE_DIR})
target_link_libraries(libAdvection
  PUBLIC
  ${Commons}
  pthread
  ${M_LIB}
  )

# Executable
add_executable(Advection sources/advect.c)
target_link_libraries(Advection libAdvection ${Commons})

# This install is really borked, so I do not want to use it
# INSTALL(               TARGETS Advection LIBRARY  DESTINATION "$ENV{HOME}/lib")
# INSTALL(               FILES ${header_files} DESTINATION "$ENV{HOME}/include")
# INSTALL(               TARGETS advect RUNTIME DESTINATION "$ENV{HOME}/bin")
