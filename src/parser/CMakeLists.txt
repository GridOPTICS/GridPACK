add_executable(greetings test/PTI23_test.cpp)
target_link_libraries(greetings ${Boost_LIBRARIES})

add_test(greetings_serial PTI23_test)

