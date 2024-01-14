include_guard()

if (TARGET armadillo::armadillo)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    armadillo
    GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
    GIT_TAG ef4736878b2d1b2bbef1a0f22e18e776f171feaf
    GIT_SHALLOW TRUE
)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not use shared lib.")
FetchContent_MakeAvailable(armadillo)

add_library(armadillo::armadillo ALIAS armadillo)
