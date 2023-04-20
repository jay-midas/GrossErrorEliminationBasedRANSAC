cmake_minimum_required(VERSION 3.24)

project(MAGSAC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_subdirectory(model)
add_subdirectory(magsac)
add_subdirectory(data)
add_subdirectory(test)