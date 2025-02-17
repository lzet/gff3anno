cmake_minimum_required(VERSION 3.5)

project(gffparser LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(MAKE_TEST "make test utility" OFF)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(Threads REQUIRED)
add_library(
    ${PROJECT_NAME} STATIC
    include/gffparser.h src/gffparser.cpp
)
target_link_libraries(${PROJECT_NAME} Threads::Threads)

if(MAKE_TEST)
    add_executable(${PROJECT_NAME}_test src/test.cpp)
    target_link_libraries(${PROJECT_NAME}_test ${PROJECT_NAME})
endif(MAKE_TEST)
