cmake_minimum_required(VERSION 3.16)

project(gff3anno LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(
    3rdparty/gffparser
)
add_subdirectory(
    3rdparty/zlib-1.3.1
)

include_directories(
    3rdparty/gffparser/include
    3rdparty/getopts
    3rdparty/bxzstr-1.2.3/include
)

add_executable(${PROJECT_NAME}
    3rdparty/getopts/getopts.cpp
    main.cpp
)
target_link_libraries(${PROJECT_NAME}
    stdc++fs
    gffparser
    zlibstatic
)

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
