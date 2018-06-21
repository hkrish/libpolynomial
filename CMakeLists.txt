
cmake_minimum_required(VERSION 3.5)

# set(LLVM_PATH "/usr/local/opt/llvm")

# set(CMAKE_C_COMPILER "${LLVM_PATH}/bin/clang")
# set(CMAKE_CXX_COMPILER "${LLVM_PATH}/bin/clang++")
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.11")

project(libPolynomial5 LANGUAGES C CXX)

set(CMAKE_C_FLAGS_RELEASE "-O3 -Wall -Wpedantic")
set(CMAKE_C_FLAGS_DEBUG "-g -Wall")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")

# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse -msse2")

include_directories(
  )

link_directories(
  )

set(libsrc
  src/Polynomial5.h
  src/Polynomial5.c
  )

add_library(polynomial5static STATIC ${libsrc})

add_library(polynomial5 SHARED ${libsrc})
set_property(TARGET polynomial5 PROPERTY POSITION_INDEPENDENT_CODE 1)

add_executable(
  tests
  examples/tests.c
  )

target_link_libraries(
  tests
  polynomial5
  "m"
  )

install(TARGETS polynomial5 polynomial5static DESTINATION lib)
install(FILES src/Polynomial5.h DESTINATION include)
