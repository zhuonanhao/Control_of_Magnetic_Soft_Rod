cmake_minimum_required(VERSION 3.11)
project(magnetic-tube)

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wno-deprecated-declarations")
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(Threads REQUIRED)
find_package(Eigen3 3.4 REQUIRED)
find_package(OpenMP)


include_directories(
        ${EIGEN_INCLUDE_DIRS}
        ${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIRS})



add_executable(simDER
        src/main.cpp
        src/world.cpp

        src/rod_mechanics/inner_forces/elasticStretchingForce.cpp
        src/rod_mechanics/inner_forces/elasticBendingForce.cpp
        src/rod_mechanics/inner_forces/elasticTwistingForce.cpp
        src/rod_mechanics/inner_forces/inertialForce.cpp

        src/rod_mechanics/external_forces/externalGravityForce.cpp
        src/rod_mechanics/external_forces/dampingForce.cpp
        # src/dampingForce.cpp
        src/rod_mechanics/external_forces/externalContactForce.cpp
        src/rod_mechanics/control_forces/externalMagneticForce.cpp

        src/rod_mechanics/elasticRod.cpp
        src/time_steppers/timeStepper.cpp
        src/initialization/setInput.cpp
        )

target_include_directories(simDER PRIVATE
    src
)

target_link_libraries(simDER PRIVATE
    ${OPENGL_LIBRARIES}
    ${GLUT_LIBRARIES}
    Eigen3::Eigen
    ${LAPACK_LIBRARIES}
    ${BLAS_LIBRARIES}

    lapack

    # gfortran
    Threads::Threads
    m
)

if (OpenMP_CXX_FOUND)
    target_link_libraries(simDER PUBLIC OpenMP::OpenMP_CXX)
endif ()
