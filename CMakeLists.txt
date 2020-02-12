cmake_minimum_required(VERSION 3.1)
project(VectorRender LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
set(Boost_USE_STATIC_LIBS ON)
find_package(Eigen REQUIRED)
find_package(Boost REQUIRED)
find_package(CGAL REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")  # enable assert
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")  # enable assert
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG}")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fopenmp -Wno-int-in-bool-context -Wno-sign-compare")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}")
    set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")
endif()

if(MSVC)
    ADD_DEFINITIONS(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /wd4267 /wd4244 /wd4018 /wd4800") 
endif()

include(${CGAL_USE_FILE})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(${CGAL_INCLUDE_DIRS})

set(
    VectorRender_SRC
    src/main.cc
    src/types.h
    src/plane_param.cc
    src/plane_param.h
    src/split_data.cc
    src/split_data.h
)

add_executable(
    vector_render
    ${VectorRender_SRC}
)

target_link_libraries(
    vector_render
    ${CGAL_LIBRARIES}
)
