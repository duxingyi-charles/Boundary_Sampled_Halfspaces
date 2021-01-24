include_guard()

include(FetchContent)
FetchContent_Declare(
    dirent
    GIT_REPOSITORY  https://github.com/tronkko/dirent.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)
option(DIRENT_BUILD_TESTS "Build bundled tests" OFF)
FetchContent_MakeAvailable(dirent)
if (WIN32)
	include_directories (${dirent_SOURCE_DIR}/include)
endif (WIN32)
