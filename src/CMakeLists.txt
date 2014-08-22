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

set(SRC SMSDSubDomains.cpp SMSDetector.cpp Carrier.cpp CarrierMobility.cpp CarrierTransport.cpp)

# Executable
add_executable(WeightFEM main.cpp ${SRC})

# Target libraries
target_link_libraries(WeightFEM ${DOLFIN_LIBRARIES} ${DOLFIN_3RD_PARTY_LIBRARIES} ${LIBRARIES})


set(GUI_HEADERS mainWindow.h)
set(GUI_SRC mainWindow.cpp)
set(GUI_UIS mainWindow.ui)

qt4_wrap_cpp(GUI_MOC ${GUI_HEADERS})
qt4_wrap_ui(GUI_UIS_HEADERS ${GUI_UIS})


# Executable
add_executable(WeightFEM-GUI main_gui.cpp ${GUI_SRC} ${GUI_MOC} ${GUI_UIS_HEADERS} )
target_link_libraries(WeightFEM-GUI ${QT_LIBRARIES})



# vim:ft=cmake


