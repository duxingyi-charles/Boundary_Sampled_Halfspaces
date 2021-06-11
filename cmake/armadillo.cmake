include_guard()

if (TARGET armadillo::armadillo)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    armadillo
    GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
    GIT_TAG aa8ee5e1c0f6288419445aed21e9aff8fba4ea1c
    GIT_SHALLOW TRUE
)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not use shared lib.")
FetchContent_MakeAvailable(armadillo)

add_library(armadillo::armadillo ALIAS armadillo)
