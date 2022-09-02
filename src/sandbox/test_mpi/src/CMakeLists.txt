
pybind11_add_module(Test test_mpi_wrap.cpp)

target_link_libraries(Test PRIVATE ${MPI_CXX_LIBRARIES})
