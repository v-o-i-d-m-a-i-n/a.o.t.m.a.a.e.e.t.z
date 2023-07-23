# cmake super build mode, root dir

cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
project ("AOTMAAEETZ" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
# set(CMAKE_BUILD_TYPE "Debug") # for vscode extension cmake tools, this can be set using GUI
## set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -std=gnu++11")
set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "${PROJECT_NAME} staged install: ${STAGED_INSTALL_PREFIX}")

# if the EP_BASE directory property is set then components of an external project are stored under the specified base
# If no PREFIX, EP_PREFIX, or EP_BASE is specified, then the default is to set PREFIX to <name>-prefix.
# https://cmake.org/cmake/help/latest/module/ExternalProject.html#command:externalproject_add
set_property(DIRECTORY PROPERTY EP_BASE ${CMAKE_BINARY_DIR}/out)

list(APPEND BOOST_COMPONENTS_REQUIRED )
list(APPEND BOOST_COMPONENTS_NOREQUIRED atomic chrono container context contract coroutine date_time exception fiber filesystem graph graph_parallel headers iostreams json locale log math mpi nowide program_options python random regex serialization stacktrace system test thread timer type_erasure wave)
set(Boost_MINIMUM_REQUIRED 1.60)
set(CMAKE_DISABLE_FIND_PACKAGE_Boost OFF) # set this to on will force superbuild
add_subdirectory(external/upstream)

set(EX_NAME "algs")

# refer https://github.com/martijnberger/blender_dependencies/blob/master/CMakeLists.txt
# ExternalProject is a cmake standard module. The module provide function ExternalProject_Add
# detail check https://cmake.org/cmake/help/latest/module/ExternalProject.html
# set(nlohmann_json_DIR ${JSON_PKGCONFIGDIR})
include(ExternalProject)

ExternalProject_Add(${EX_NAME}
	DEPENDS
		boost_external
	SOURCE_DIR
		${CMAKE_CURRENT_LIST_DIR}/src
	CMAKE_ARGS
		-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
		-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
		-DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
		-DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
		-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
	CMAKE_CACHE_ARGS
		-DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
		-DBoost_MINIMUM_REQUIRED:STRING=${Boost_MINIMUM_REQUIRED}
		-DEX_NAME:STRING=${EX_NAME}
		-DSTAGED_INSTALL_PREFIX:PATH=${STAGED_INSTALL_PREFIX}
		-DCMAKE_INCLUDE_PATH:PATH=${BOOST_INCLUDEDIR} #${STAGED_INSTALL_PREFIX}/boost/include 
		-DCMAKE_LIBRARY_PATH:PATH=${BOOST_LIBRARYDIR} 
		# -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH} ${BOOST_INCLUDEDIR}
	BUILD_ALWAYS
		1
	# INSTALL_COMMAND
	# 	""
)

add_executable(${PROJECT_NAME} IMPORTED)

set_property(TARGET ${PROJECT_NAME} PROPERTY IMPORTED_LOCATION "${STAGED_INSTALL_PREFIX}/bin/${EX_NAME}")

add_dependencies(${PROJECT_NAME} ${EX_NAME})
