cmake_minimum_required(VERSION 3.1)
project(navtoolbox)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(TARGET navtoolbox)
find_package(Eigen3 NO_MODULE)
add_definitions(${EIGEN3_DEFINITIONS})
include_directories(${EIGEN3_INCLUDE_DIRS})


file(GLOB HEADER_FILES ./include/*.h)
file(GLOB SOURCE_FILES ./src/*.cpp)

add_library(${TARGET} STATIC
            ${SOURCE_FILES}
            ${HEADER_FILES})
target_include_directories(${TARGET} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_link_libraries(${TARGET} Eigen3::Eigen3)

set_property(TARGET navtoolbox PROPERTY CXX_STANDARD 11)
set_property(TARGET navtoolbox PROPERTY CXX_STANDARD_REQUIRED ON)
