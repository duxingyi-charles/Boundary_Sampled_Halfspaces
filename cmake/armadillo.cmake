include_guard()

if (TARGET armadillo::armadillo)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    armadillo
    GIT_REPOSITORY https://gitlab.com/conradsnicta/armadillo-code.git
    GIT_TAG 7b9215a35870e2ff3a7522eacf815f8f4a8489f7
    GIT_SHALLOW TRUE
)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not use shared lib.")
FetchContent_MakeAvailable(armadillo)

add_library(armadillo::armadillo ALIAS armadillo)
