if(TARGET nlohmann_json::nlohmann_json)
    return()
endif()

set(NLOHMANNJSON_VERSION "v3.7.3")

include(FetchContent)
FetchContent_Declare(
    nlohmann_json
    URL "https://github.com/nlohmann/json/releases/download/${NLOHMANNJSON_VERSION}/include.zip"
    URL_HASH SHA256=87b5884741427220d3a33df1363ae0e8b898099fbc59f1c451113f6732891014
)
FetchContent_MakeAvailable(nlohmann_json)

add_library(nlohmann_json INTERFACE)
target_include_directories(nlohmann_json INTERFACE
    ${nlohmann_json_SOURCE_DIR}/include)
add_library(nlohmann::json ALIAS nlohmann_json)
