# Install script for directory: /Users/xmt/code/github/spurt

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/xmt/XAVIER.default_install_dir")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE DIRECTORY FILES "/Users/xmt/code/github/spurt/bin" USE_SOURCE_PERMISSIONS)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE DIRECTORY FILES "/Users/xmt/code/github/spurt/scripts" USE_SOURCE_PERMISSIONS)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/utils" TYPE DIRECTORY FILES "/Users/xmt/code/github/spurt/utils/GLSLShaderCode" USE_SOURCE_PERMISSIONS)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/xmt/code/github/spurt/build/boundary-aware-rectgrid/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/third_party/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/crease/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/data/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/format/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/flow/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/geometry/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/graphics/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/image/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/learning/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/maps/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/math/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/misc/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/reconstruction/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/tensor/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/topology/cmake_install.cmake")
  include("/Users/xmt/code/github/spurt/build/utils/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/xmt/code/github/spurt/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
