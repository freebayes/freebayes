# Specify the minimum version for CMake

cmake_minimum_required(VERSION 3.1)

# Project's name
project(freebayes)
# We build using c++14
set(CMAKE_CXX_STANDARD 14)

# Use all standard-compliant optimizations
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")

# Add external projects
include(ExternalProject)

# Find our dependencies
include(FindZLIB)
include(FindBZip2)
include(FindLibLZMA)
include(FindThreads)

ExternalProject_Add(vcflib
  SOURCE_DIR "${CMAKE_SOURCE_DIR}/vcflib"
  CMAKE_ARGS "${CMAKE_ARGS};-DBUILD_ONLY_LIB=1;-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>")
ExternalProject_Get_property(vcflib INSTALL_DIR)
set(vcflib_INCLUDE "${INSTALL_DIR}/include")
set(vcflib_LIB "${INSTALL_DIR}/lib")

ExternalProject_Add(seqlib
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/SeqLib
  BUILD_IN_SOURCE 1 # why do we have to build in source?
  CONFIGURE_COMMAND ./configure
  BUILD_COMMAND $(MAKE) CXXFLAGS='-std=c++14'
  INSTALL_COMMAND $(MAKE) install)
ExternalProject_Get_property(seqlib SOURCE_DIR)
set(seqlib_INCLUDE "${SOURCE_DIR}")
set(seqlib_LIB "${SOURCE_DIR}/lib")

#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

# set up our target executable and specify its dependencies and includes
add_library(freebayes_OBJS OBJECT
  ${CMAKE_SOURCE_DIR}/src/Contamination.cpp
  ${CMAKE_SOURCE_DIR}/src/Marginals.cpp
  ${CMAKE_SOURCE_DIR}/src/DataLikelihood.cpp
  ${CMAKE_SOURCE_DIR}/src/NonCall.cpp
  ${CMAKE_SOURCE_DIR}/src/Genotype.cpp
  ${CMAKE_SOURCE_DIR}/src/Utility.cpp
  ${CMAKE_SOURCE_DIR}/src/AlleleParser.cpp
  ${CMAKE_SOURCE_DIR}/src/Dirichlet.cpp
  ${CMAKE_SOURCE_DIR}/src/split.cpp
  ${CMAKE_SOURCE_DIR}/src/Sample.cpp
  ${CMAKE_SOURCE_DIR}/src/SegfaultHandler.cpp
  ${CMAKE_SOURCE_DIR}/src/BedReader.cpp
  ${CMAKE_SOURCE_DIR}/src/Multinomial.cpp
  ${CMAKE_SOURCE_DIR}/src/LeftAlign.cpp
  ${CMAKE_SOURCE_DIR}/src/CNV.cpp
  ${CMAKE_SOURCE_DIR}/src/ResultData.cpp
  ${CMAKE_SOURCE_DIR}/src/Bias.cpp
  ${CMAKE_SOURCE_DIR}/src/Result.cpp
  ${CMAKE_SOURCE_DIR}/src/Allele.cpp
  ${CMAKE_SOURCE_DIR}/src/IndelAllele.cpp
  ${CMAKE_SOURCE_DIR}/src/Fasta.cpp
  ${CMAKE_SOURCE_DIR}/src/Ewens.cpp
  ${CMAKE_SOURCE_DIR}/src/Parameters.cpp
  ${CMAKE_SOURCE_DIR}/src/BGZF.cpp)

set(freebayes_DEPS
  vcflib
  seqlib)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/version_git.txt
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/scripts/update_version.sh ${CMAKE_CURRENT_SOURCE_DIR})
add_custom_target(SetVersion ALL
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/version_git.txt)

add_dependencies(freebayes_OBJS ${freebayes_DEPS} SetVersion)

set(freebayes_INCLUDES
  "${CMAKE_SOURCE_DIR}/src"
  "${CMAKE_SOURCE_DIR}/ttmath"
  "${vcflib_INCLUDE}"
  "${seqlib_INCLUDE}"
  "${seqlib_INCLUDE}/htslib")

set(freebayes_LIBS
  "${vcflib_LIB}/libvcflib.a"
  "${seqlib_LIB}/libseqlib.a"
  "${seqlib_LIB}/libhts.a"
  ${CMAKE_THREAD_LIBS_INIT}
  ${ZLIB_LIBRARIES}
  ${LIBLZMA_LIBRARIES}
  ${BZIP2_LIBRARIES})

target_include_directories(freebayes_OBJS PUBLIC ${freebayes_INCLUDES})
set_target_properties(freebayes_OBJS PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

add_executable(freebayes
  $<TARGET_OBJECTS:freebayes_OBJS>
  ${CMAKE_SOURCE_DIR}/src/freebayes.cpp)
target_include_directories(freebayes PUBLIC ${freebayes_INCLUDES})
target_link_libraries(freebayes ${freebayes_LIBS})
set_target_properties(freebayes PROPERTIES OUTPUT_NAME "freebayes")

add_executable(bamleftalign
  $<TARGET_OBJECTS:freebayes_OBJS>
  ${CMAKE_SOURCE_DIR}/src/bamleftalign.cpp)
target_include_directories(bamleftalign PUBLIC ${freebayes_INCLUDES})
target_link_libraries(bamleftalign ${freebayes_LIBS})
set_target_properties(bamleftalign PROPERTIES OUTPUT_NAME "bamleftalign")

install(TARGETS freebayes DESTINATION bin)
install(TARGETS bamleftalign DESTINATION bin)

if (APPLE)
elseif (TRUE)
  if (BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
  endif()
endif()
