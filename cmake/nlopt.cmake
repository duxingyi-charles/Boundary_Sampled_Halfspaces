include_guard()

if (TARGET nlopt::nlopt)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    nlopt
    GIT_REPOSITORY git@github.com:stevengj/nlopt.git
    GIT_TAG        v2.7.0
    GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(nlopt)

add_library(nlopt::nlopt ALIAS nlopt)
