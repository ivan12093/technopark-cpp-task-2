cmake_minimum_required(VERSION 3.18)
project(technopark_cpp_task_2)

set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES "src/*.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")

add_library(${PROJECT_NAME}_lib STATIC ${SOURCES})
target_include_directories(${PROJECT_NAME}_lib PUBLIC include)

add_executable(technopark_cpp_task_2
        src/main.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC include)

target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_lib)

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

if (BUILD_COVERAGE)
    add_coverage(${PROJECT_NAME})
    add_coverage(${PROJECT_NAME}_lib)
endif()

