#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(GRATE)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" OFF)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitio
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)
find_package(ROOT CONFIG REQUIRED)
find_package(ROOT COMPONENTS MathMore)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

include_directories(ROOT_BUG)

include_directories(${PROJECT_SOURCE_DIR}/include  
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS})
add_subdirectory(TGlauber)
link_directories(${ROOT_LIBRARY_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/SMM/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

#add_library(TGlauNucleon SHARED ${PROJECT_SOURCE_DIR}/src/TGlauNucleon.cc ${PROJECT_SOURCE_DIR}/include/TGlauNucleon.hh)
#add_library(TGlauNucleus SHARED ${PROJECT_SOURCE_DIR}/src/TGlauNucleus.cc ${PROJECT_SOURCE_DIR}/include/TGlauNucleus.hh)
#add_library(TGlauberMC SHARED ${PROJECT_SOURCE_DIR}/src/TGlauberMC.cc ${PROJECT_SOURCE_DIR}/include/TGlauberMC.hh)
add_executable(GRATE GRATE.cc ${sources} ${headers})
#target_link_libraries(GRATE TGlauNucleon TGlauNucleus TGlauberMC ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} TGlauberLib)
target_link_libraries(GRATE ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} TGlauberLib)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS GRATE DESTINATION bin)

