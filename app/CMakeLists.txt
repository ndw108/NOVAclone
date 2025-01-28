add_executable (TG ./TG/TG.cpp)
target_link_libraries (TG  NOVA ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${HEFFTE_LIBRARIES} ${LIBS}) 
install(TARGETS TG DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_executable (TG_Poisson ./TG_Piosson/TG.cpp)
target_link_libraries (TG_Poisson  NOVA ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LIBS}) 
install(TARGETS TG_Poisson DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_executable (OT ./OT/OT.cpp)
target_link_libraries (OT  NOVA ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${HEFFTE_LIBRARIES} ${LIBS})
install(TARGETS OT DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_executable (OT_Poisson ./OT_Poisson/OT.cpp)
target_link_libraries (OT_Poisson  NOVA ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${HEFFTE_LIBRARIES} ${LIBS})
install(TARGETS OT DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_executable (shearLayer ./shearLayer/shearLayer.cpp)
target_link_libraries (shearLayer  NOVA ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LIBS}) 
install(TARGETS shearLayer DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_executable (mhdjet ./mhdjet/mhdjet.cpp)
target_link_libraries (mhdjet  NOVA ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${LIBS}) 
install(TARGETS mhdjet DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

