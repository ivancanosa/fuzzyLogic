cmake_minimum_required(VERSION 3.10)
project(FuzzyLogic VERSION 1.0 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest) 
enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)

add_library(fuzzyLogic SHARED FuzzyLogic.cpp)
target_include_directories(fuzzyLogic PUBLIC 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
add_library(fuzzyLogic::fuzzyLogic ALIAS fuzzyLogic)

add_subdirectory(test)
add_subdirectory(example)

# Install the target
install(TARGETS fuzzyLogic
        EXPORT fuzzyLogicTarget
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

#install the includes
install(FILES FuzzyLogic.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#install the cmake file for importing
install(EXPORT fuzzyLogicTarget
        FILE fuzzyLogicTarget.cmake
        NAMESPACE fuzzyLogic::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fuzzyLogic
)

# Install the config file
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/fuzzyLogicConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fuzzyLogic
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/fuzzyLogicConfig.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fuzzyLogic
)
