cmake_minimum_required(VERSION 3.4)

project(test)

add_executable(${PROJECT_NAME} main.cc)
target_compile_options(${PROJECT_NAME} PUBLIC -g)

target_link_libraries(${PROJECT_NAME} model data magsac)