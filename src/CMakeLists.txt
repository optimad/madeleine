#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
set(MADELEINE_EXECUTABLE_NAME madeleine CACHE INTERNAL "Executable name of the solver" FORCE)

include_directories(${PETSC_INCLUDES})

file(GLOB sources "*.cpp")
add_executable(${MADELEINE_EXECUTABLE_NAME} ${sources})

target_link_libraries(${MADELEINE_EXECUTABLE_NAME} ${BITPIT_LIBRARIES})
target_link_libraries(${MADELEINE_EXECUTABLE_NAME} ${LIBXML2_LIBRARIES})
target_link_libraries(${MADELEINE_EXECUTABLE_NAME} ${PETSC_LIBRARIES})
target_link_libraries(${MADELEINE_EXECUTABLE_NAME} ${LAPACKE_LIBRARIES})
target_link_libraries(${MADELEINE_EXECUTABLE_NAME} ${LAPACK_LIBRARIES})

INSTALL (TARGETS ${MADELEINE_EXECUTABLE_NAME} DESTINATION bin)

