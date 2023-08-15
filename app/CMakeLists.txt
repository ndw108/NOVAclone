add_executable (TG ./TG/TG.cpp)
target_link_libraries (TG  NOVA ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LIBS}) 
install(TARGETS TG DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
