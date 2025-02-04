enable_testing()
include(GoogleTest)

include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.14.0
)
FetchContent_MakeAvailable(googletest)

add_executable(test_method_progonky test.cpp)
target_link_libraries(test_method_progonky gtest gtest_main ${PROJECT_NAME})
set_target_properties(test_method_progonky PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

gtest_discover_tests(test_method_progonky)

