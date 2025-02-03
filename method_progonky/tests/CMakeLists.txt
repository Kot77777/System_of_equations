enable_testing()
include(GoogleTest)

add_executable(test_method_progonky test.cpp)
target_link_libraries(test_method_progonky gtest gtest_main ${PROJECT_NAME})
set_target_properties(test_method_progonky PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

gtest_discover_tests(test_method_progonky)

