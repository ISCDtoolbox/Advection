# Install script for directory: /Users/dapogny/Advection

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/dapogny/lib/libAdvection.dylib")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/dapogny/lib" TYPE SHARED_LIBRARY FILES "/Users/dapogny/Advection/build/libAdvection.dylib")
  if(EXISTS "$ENV{DESTDIR}/Users/dapogny/lib/libAdvection.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/dapogny/lib/libAdvection.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/dapogny/lib"
      "$ENV{DESTDIR}/Users/dapogny/lib/libAdvection.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/dapogny/lib"
      "$ENV{DESTDIR}/Users/dapogny/lib/libAdvection.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/dapogny/lib/libAdvection.dylib")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/dapogny/include/ad_calls.h;/Users/dapogny/include/advect.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/dapogny/include" TYPE FILE FILES
    "/Users/dapogny/Advection/sources/ad_calls.h"
    "/Users/dapogny/Advection/sources/advect.h"
    )
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/dapogny/bin/advect")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/Users/dapogny/bin" TYPE EXECUTABLE FILES "/Users/dapogny/Advection/build/advect")
  if(EXISTS "$ENV{DESTDIR}/Users/dapogny/bin/advect" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/dapogny/bin/advect")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/dapogny/lib"
      "$ENV{DESTDIR}/Users/dapogny/bin/advect")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/dapogny/Advection/build"
      "$ENV{DESTDIR}/Users/dapogny/bin/advect")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/dapogny/lib"
      "$ENV{DESTDIR}/Users/dapogny/bin/advect")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" "$ENV{DESTDIR}/Users/dapogny/bin/advect")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/dapogny/Advection/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
