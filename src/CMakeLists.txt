cmake_minimum_required(VERSION 2.8)

# Set verbose output while testing CMake
#set(CMAKE_VERBOSE_MAKEFILE 1)
# Set CMake behavior
cmake_policy(SET CMP0004 OLD)

# Get DOLFIN configuration data (DOLFINConfig.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
find_package(DOLFIN)

# Default build type (can be overridden by user)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING
    "Choose the type of build, options are: Debug MinSizeRel Release RelWithDebInfo." FORCE)
endif()

# Compiler definitions
add_definitions(${DOLFIN_CXX_DEFINITIONS})

# Compiler flags
set(CMAKE_CXX_FLAGS "${DOLFIN_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")

# Include directories
include_directories(${DOLFIN_INCLUDE_DIRS})
include_directories(SYSTEM ${DOLFIN_3RD_PARTY_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})


set(SRC SMSDSubDomains.cpp SMSDetector.cpp
    Carrier.cpp CarrierMobility.cpp CarrierTransport.cpp
    CarrierCollection.cpp utilities.cpp
    qcustomplot.cpp qcustomplot.h)

set(HEADERS qcustomplot.h)

# Run MOC from Qt over headers
# OPTIONS are a  workaround to moc+boost problem(QT bug #22829).
qt4_wrap_cpp(NONGUI_MOC ${HEADERS} OPTIONS   -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_AND_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_AND_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_OR_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_OR_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_XOR_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_XOR_HPP_INCLUDED
                                              -DBOOST_TT_HAS_DIVIDES_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_DIVIDES_HPP_INCLUDED
                                              -DBOOST_TT_HAS_EQUAL_TO_HPP_INCLUDED
                                              -DBOOST_TT_HAS_GREATER_EQUAL_HPP_INCLUDED
                                              -DBOOST_TT_HAS_GREATER_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LEFT_SHIFT_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LEFT_SHIFT_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LESS_EQUAL_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LESS_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LOGICAL_AND_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LOGICAL_OR_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MINUS_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MINUS_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MODULUS_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MODULUS_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MULTIPLIES_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MULTIPLIES_HPP_INCLUDED
                                              -DBOOST_TT_HAS_NOT_EQUAL_TO_HPP_INCLUDED
                                              -DBOOST_TT_HAS_PLUS_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_PLUS_HPP_INCLUDED
                                              -DBOOST_TT_HAS_RIGHT_SHIFT_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_RIGHT_SHIFT_HPP_INCLUDED
                                              )


# Executable
add_executable(WeightFEM main.cpp ${SRC} ${HEADERS} ${NONGUI_MOC})

# Target libraries
target_link_libraries(WeightFEM ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${LIBRARIES} ${QT_LIBRARIES})


set(GUI_HEADERS mainWindow.h qcustomplot.h)
set(GUI_SRC mainWindow.cpp SMSDSubDomains.cpp SMSDetector.cpp Carrier.cpp CarrierMobility.cpp CarrierTransport.cpp CarrierCollection.cpp utilities.cpp qcustomplot.cpp)
set(GUI_UIS mainWindow.ui)


# Run MOC from Qt over headers
# OPTIONS are a  workaround to moc+boost problem(QT bug #22829).
qt4_wrap_cpp(GUI_MOC ${GUI_HEADERS} OPTIONS   -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_AND_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_AND_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_OR_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_OR_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_XOR_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_BIT_XOR_HPP_INCLUDED
                                              -DBOOST_TT_HAS_DIVIDES_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_DIVIDES_HPP_INCLUDED
                                              -DBOOST_TT_HAS_EQUAL_TO_HPP_INCLUDED
                                              -DBOOST_TT_HAS_GREATER_EQUAL_HPP_INCLUDED
                                              -DBOOST_TT_HAS_GREATER_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LEFT_SHIFT_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LEFT_SHIFT_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LESS_EQUAL_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LESS_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LOGICAL_AND_HPP_INCLUDED
                                              -DBOOST_TT_HAS_LOGICAL_OR_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MINUS_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MINUS_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MODULUS_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MODULUS_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MULTIPLIES_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_MULTIPLIES_HPP_INCLUDED
                                              -DBOOST_TT_HAS_NOT_EQUAL_TO_HPP_INCLUDED
                                              -DBOOST_TT_HAS_PLUS_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_PLUS_HPP_INCLUDED
                                              -DBOOST_TT_HAS_RIGHT_SHIFT_ASSIGN_HPP_INCLUDED
                                              -DBOOST_TT_HAS_RIGHT_SHIFT_HPP_INCLUDED
                                              )

qt4_wrap_ui(GUI_UIS_HEADERS ${GUI_UIS})


# Executable
add_executable(WeightFEM-GUI main_gui.cpp ${GUI_SRC} ${GUI_MOC} ${GUI_UIS_HEADERS} )
target_link_libraries(WeightFEM-GUI ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${LIBRARIES} ${QT_LIBRARIES} QVTK )



# vim:ft=cmake


