file(GLOB SOURCES "./*/*.cpp")
file(GLOB SOURCES2 "./*/*/*.cpp")
add_library (NOVA ${SOURCES} ${SOURCES2})
install( TARGETS NOVA DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
