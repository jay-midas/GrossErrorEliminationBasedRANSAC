cmake_minimum_required(VERSION 3.24)

project(magsac)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB HEADERS "*.h")
file(GLOB SOURCES "*.cc")

add_library(${PROJECT_NAME} ${HEADERS} ${SOURCES})

# add_executable(${PROJECT_NAME} ${HEADERS} ${SOURCES})

set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(${PROJECT_NAME} INTERFACE "..")
target_link_libraries(${PROJECT_NAME} PUBLIC model data)
