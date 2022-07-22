# CMakeLists.txt for the "event" package. It creates a library and a main program. If ROOT is not installed in a default system location you need to tell CMake where to find it. Sourcing `thisroot.sh` already sets the required environment variables. Otherwise, you must tell the build system where to look 
# for ROOT, for example by passing `-DROOT_DIR="/path/to/root/installation` at CMake configuration time.



   cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
   project(trigger)

# Locate the ROOT package and define a number of useful targets and variables.
   find_package(ROOT 6.20 CONFIG REQUIRED)
   #set(CMAKE_BUILD_TYPE Debug)
   #find_package(FFTW REQUIRED)
#   include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
   include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include )
include("${ROOT_USE_FILE}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -g")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg -g")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg -g")

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -fno-omit-frame-pointer ")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ggdb -fno-omit-frame-pointer")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ggdb -fno-omit-frame-pointer")



# Create a shared library.
# Passing cmake targets such as `ROOT::RIO` as dependencies (rather than plain
# library names for example via ${ROOT_LIBRARIES}) ensures that properties such as required
# include directories and C++ standard are propagated to our libraries or executables.
# Note: To ensure compatibility with Cling, targets *must* be compiled using the
# same C++ standard as ROOT was compiled with.
#   add_library(Event SHARED Event.cxx )
#   target_link_libraries(Event PUBLIC ROOT::RIO ROOT::Net)

# Create the main program using the library.
#   set(SOURCES FTPair.cc trigger_shared.cc ${CMAKE_SOURCE_DIR}/include/trigger.h)
#   target_sources(triggerStandalone PRIVATE FTPair.cc trigger_standalone.cc)
   add_executable(triggerStandalone trigger_standalone.cc)
   target_include_directories(triggerStandalone PRIVATE ${CMAKE_SOURCE_DIR}/include)
   target_sources(triggerStandalone PRIVATE FTPair.cc trigger_shared.cc)
   target_link_libraries(triggerStandalone PUBLIC fftw3 fftw ROOT::Core ROOT::Tree ROOT::Graf3d ROOT::Graf ROOT::Hist) 
   #target_link_libraries(Main PUBLIC fftw3 fftw ROOT::Core ROOT::MathCore ROOT::Gpad ROOT::Graf3d ROOT::Graf ROOT::Hist ROOT::Imt ROOT::Matrix ROOT::Net ROOT::Physics ROOT::RIO ROOT::ROOTDataFrame ROOT::ROOTVecOps ROOT::Rint ROOT::Thread ROOT::TreePlayer ROOT::Tree) 



   target_link_libraries(triggerStandalone PUBLIC ${ROOT_LIBRARIES} ${ROOT_EXE_LINKER_FLAGS})
