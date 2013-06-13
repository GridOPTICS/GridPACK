add_executable(parser1 test/PTI23_test.cpp)
target_link_libraries(parser1 ${Boost_LIBRARIES})

add_test(parser1 PTI23_test)

