include_guard()

if (TARGET armadillo::armadillo)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    armadillo
    GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
    GIT_TAG 6a312a98f6ca5450b6ab04361abd6d36ec0759f7
    GIT_SHALLOW TRUE
)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not use shared lib.")
FetchContent_MakeAvailable(armadillo)

add_library(armadillo::armadillo ALIAS armadillo)
