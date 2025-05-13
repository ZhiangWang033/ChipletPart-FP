cmake_minimum_required(VERSION 3.16.2)
project(floorlet)

find_package(Threads)

find_package(Boost)


include_directories(src
    ${Boost_INCLUDE_DIRS}
    )

file(GLOB SOURCES "src/*.cpp")

add_executable(exe ${SOURCES})

TARGET_LINK_LIBRARIES(exe
    ${CMAKE_THREAD_LIBS_INIT}
    ${Boost_LIBRARIES}
    )

