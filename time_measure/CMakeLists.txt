file(GLOB_RECURSE COMPARE_SOURCES compare_*.cpp)

foreach (COMPARE_FILE IN LISTS COMPARE_SOURCES)
    get_filename_component(COMPARE_FILE_NAME ${COMPARE_FILE} NAME_WE)
    add_executable(${COMPARE_FILE_NAME} ${COMPARE_FILE})
    target_link_libraries(${COMPARE_FILE_NAME} ${PROJECT_NAME})
    set_target_properties(${COMPARE_FILE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endforeach ()