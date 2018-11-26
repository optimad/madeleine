#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

#Name your project here
#set(MADELEINE_EXECUTABLE_NAME madeleine CACHE INTERNAL "Executable name of the solver" FORCE)
set(MADELEINE_LIBRARY madeleine CACHE INTERNAL "Library name" FORCE)

include_directories(${PETSC_INCLUDES})

file(GLOB sources "*.cpp")
add_library(${MADELEINE_LIBRARY} ${sources})

target_link_libraries(${MADELEINE_LIBRARY} ${MADELEINE_EXTERNAL_LIBRARIES})

INSTALL (TARGETS ${MADELEINE_LIB_NAME} DESTINATION lib)

# Data files
#add_custom_command(
#    TARGET ${MADELEINE_EXECUTABLE_NAME} PRE_BUILD
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/sphere1.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/sphere1.stl"
#    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/data/sphere2.stl" "${CMAKE_CURRENT_BINARY_DIR}/data/sphere2.stl"
#)

