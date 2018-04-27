cmake_minimum_required(VERSION 3.0.0)

project(kPers)

# Set default build type to release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# warning flags
set(WARNING_FLAGS "-Wredundant-decls -Wcast-align -Wmissing-declarations -Wmissing-include-dirs -Wswitch-enum -Wswitch-default -Wextra -Wall -Werror -Winvalid-pch -Wredundant-decls -Wformat=2 -Wmissing-format-attribute -Wformat-nonliteral")
set(ADD_FLAGS_DEBUG "${WARNING_FLAGS} -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_SANITIZE_VECTOR")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${ADD_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${ADD_FLAGS_DEBUG}")

# optimisation flags
set(ADD_FLAGS_RELEASE "-O3 -flto -march=native -mtune=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ADD_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${ADD_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${ADD_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${ADD_FLAGS_RELEASE}")

# Seqan section
find_package(ZLIB)
find_package(BZip2)
find_package(OpenMP)
set(SeqAn_FIND_QUIETLY TRUE)
find_package(SeqAn)

include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

include_directories(inc)

file(
    GLOB_RECURSE
    src_files
    src/*
    )

add_executable(${PROJECT_NAME} ${src_files})
target_link_libraries(${PROJECT_NAME} ${SEQAN_LIBRARIES})
