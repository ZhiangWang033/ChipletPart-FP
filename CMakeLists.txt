cmake_minimum_required(VERSION 3.16.2)
project(chiplet-fp)

find_package(Threads REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)


include_directories(
    src
    ${Boost_INCLUDE_DIRS}
)

file(GLOB SOURCES "src/*.cpp")

add_executable(exe ${SOURCES})

TARGET_LINK_LIBRARIES(exe
    ${CMAKE_THREAD_LIBS_INIT}
    ${Boost_LIBRARIES}
    OpenMP::OpenMP_CXX
)

