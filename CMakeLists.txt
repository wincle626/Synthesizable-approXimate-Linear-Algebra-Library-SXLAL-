## cmake options
cmake_minimum_required( VERSION 3.4 )
if( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8 )
    message( FATAL_ERROR "GCC 4.8 or later is required." )
endif()
# set c++ language standard
set( CMAKE_CXX_STANDARD 14 )
message( "Using c++14 standard" )
# set cmake configuration directory
set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" )


## project options
set( PROJECT_NAME "Synthesizable approXimate Linear Algebra Library ( SXLAL )" )
set( PROJECT_LANGUAGE CXX )
set( PROJECT_VERSION 0.0.1)
set( PROJECT_DESCRIPTION "This is a initial test project for proximal grediant descent algorithm" )
project( XPRECISION_LINEAR_ALGEBRA
	 VERSION ${PROJECT_VERSION}
	 LANGUAGES ${PROJECT_LANGUAGE}
	 DESCRIPTION ${PROJECT_DESCRIPTION}
)
message( "${PROJECT_DESCRIPTION} of version ${PROJECT_VERSION} using ${PROJECT_LANGUAGE}" )



## project folder
if( NOT HEAD_DIR )
	set( HEAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/header )
	set( XILINX_HEADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/Xilinx)
	set( FLOATX_HEADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/FloatX)
	set( SOFTPOSIT_HEADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/SoftPosit/include)
	set( SOFTPOSIT_SSE_HEADERS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/SoftPosit/8086-SSE)
endif()
if( NOT SRC_DIR )
	set( SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/source )
endif()
aux_source_directory( ${CMAKE_CURRENT_SOURCE_DIR}/external/SoftPosit/ SOFTPOSIT_SRC)
add_subdirectory( ${HEAD_DIR} )
add_subdirectory( ${SRC_DIR} )
include_directories( ${HEAD_DIR} )
include_directories( ${SOFTPOSIT_HEADERS_DIR} )
include_directories( ${SOFTPOSIT_SSE_HEADERS_DIR} )
set( PROJECT_SRC
     ${SOFTPOSIT_SRC}
     source/common.cpp
     source/fxpt_bitwise.cpp
     source/fft.cpp
)
set( PROJECT_INC
     header/common.hpp
     header/data.hpp
     header/eigen_algebra.hpp
     header/fpt_algebra.hpp
     header/softposit_algebra.hpp
     header/xfxpt_algebra.hpp
     header/xfxpt.hpp
     header/matplotlibcpp.hpp
     header/fpt_bitwise.hpp
     header/fxpt_bitwise.hpp
     header/softapprox.hpp
     header/fft.hpp
     header/qrd.hpp
)



## project dependence 
find_package( Eigen REQUIRED )
if( EIGEN_FOUND )
	message( "Found Eigen header directory: ${EIGEN_INCLUDE_DIRS}" )
	message( "Found Eigen library directory: ${EIGEN_LIBRARIES}" )
endif()
find_package( FFTW REQUIRED )
if( FFTW_FOUND )
	message( "Found Eigen header directory: ${FFTW_INCLUDE_DIRS}" )
	message( "Found Eigen library directory: ${FFTW_LIBRARIES}" )
endif()
find_package(PythonInterp)
if ( PYTHONINTERP_FOUND )
	message( "Found Python interpretor directory: ${PYTHON_EXECUTABLE}" )
	message( "Found Python interpretor version: ${PYTHON_VERSION_STRING}" )
endif()
find_package(PythonLibs)
if( PYTHONLIBS_FOUND )
	message( "Found Python header directory: ${PYTHON_INCLUDE_DIRS}" )
	message( "Found Python library directory: ${PYTHON_LIBRARIES}" )
endif()



## project build options
# Build type ?
include( "BuildsystemFunctions" )
mrg_setup_cmake_build_type( "Release" )
mrg_setup_executable_rpath()
# Support RPATH ?
option( ENABLE_RPATH "Enable RPATH for this library" ON )
mark_as_advanced( ENABLE_RPATH )
include( AddInstallRPATHSupport )
add_install_rpath_support( BIN_DIRS "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}"
                           LIB_DIRS "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
                           INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}"
                           DEPENDS ENABLE_RPATH
                           USE_LINK_PATH )
# compiler option
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# Add the uninstall target
include(AddUninstallTarget)
# install locally
option(INSTALL_LOCAL "Make the installation at local folder instead of system directory" ON)

## project library
add_library( sxlal ${PROJECT_SRC} 
		   ${PROJECT_INC} 
		   ${PYTHON_INCLUDE_DIRS} 
		   ${SOFTPOSIT_HEADERS_DIR}
		   ${SOFTPOSIT_SSE_HEADERS_DIR})
target_include_directories( sxlal SYSTEM
	PUBLIC ${EIGEN_INCLUDE_DIRS}
	PUBLIC ${FFTW_INCLUDE_DIRS}
	PUBLIC ${PYTHON_INCLUDE_DIRS}
	PUBLIC ${HEAD_DIR}
	PUBLIC ${XILINX_HEADERS_DIR}
	PUBLIC ${FLOATX_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_SSE_HEADERS_DIR}
)
target_link_libraries( sxlal
	PUBLIC m
	PUBLIC fftw3
	PUBLIC pthread
	PUBLIC python2.7
)

## project executable
# fixed point sqrt
add_executable( xfxpt_sqrt examples/c/xsqrt/xfxpt_sqrt.cpp examples/c/xsqrt/xfxpt_sqrt.hpp )
target_include_directories( xfxpt_sqrt SYSTEM
	PUBLIC ${EIGEN_INCLUDE_DIRS}
	PUBLIC ${PYTHON_INCLUDE_DIRS}
	PUBLIC ${HEAD_DIR}
	PUBLIC ${XILINX_HEADERS_DIR}
	PUBLIC ${FLOATX_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_SSE_HEADERS_DIR}
)
target_link_libraries( xfxpt_sqrt
	PUBLIC m
)
# matrix multiplication
add_executable( mm examples/c/mm/mm.cpp examples/c/mm/mm.hpp )
target_include_directories( mm SYSTEM
	PUBLIC ${EIGEN_INCLUDE_DIRS}
	PUBLIC ${PYTHON_INCLUDE_DIRS}
	PUBLIC ${HEAD_DIR}
	PUBLIC ${XILINX_HEADERS_DIR}
	PUBLIC ${FLOATX_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_SSE_HEADERS_DIR}
)
target_link_libraries( mm
	PUBLIC m
)
# proximal gradient descend (PGD)
add_executable( pgd examples/c/pgd/pgd_eigend.cpp examples/c/pgd/pgd_eigend.hpp
		    examples/c/pgd/pgd_eigenf.cpp examples/c/pgd/pgd_eigenf.hpp
		    examples/c/pgd/pgd_gdouble.cpp examples/c/pgd/pgd_gdouble.hpp
		    examples/c/pgd/pgd_gfloat.cpp examples/c/pgd/pgd_gfloat.hpp
		    examples/c/pgd/pgd_xfpt.cpp examples/c/pgd/pgd_xfpt.hpp
		    examples/c/pgd/pgd_xfxpt.cpp examples/c/pgd/pgd_xfxpt.hpp
		    examples/c/pgd/pgd_softposit8.cpp examples/c/pgd/pgd_softposit8.hpp
		    examples/c/pgd/pgd_softposit16.cpp examples/c/pgd/pgd_softposit16.hpp
		    examples/c/pgd/pgd_softposit32.cpp examples/c/pgd/pgd_softposit32.hpp
		    examples/c/pgd/pgd.cpp examples/c/pgd/pgd.hpp
		    source/common.cpp header/common.hpp header/data.hpp
     		    ${SOFTPOSIT_SRC}
)
target_include_directories( pgd SYSTEM
	PUBLIC ${EIGEN_INCLUDE_DIRS}
	PUBLIC ${PYTHON_INCLUDE_DIRS}
	PUBLIC ${HEAD_DIR}
	PUBLIC ${XILINX_HEADERS_DIR}
	PUBLIC ${FLOATX_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_SSE_HEADERS_DIR}
)
target_link_libraries( pgd
	PUBLIC m
	PUBLIC pthread
	PUBLIC python2.7
)
# alternating direction method of multipliers (ADMM) - ali's case
add_executable( admm examples/c/admm/admm.cpp examples/c/admm/admm.hpp
		    source/common.cpp header/common.hpp header/data.hpp
     		    ${SOFTPOSIT_SRC}
)
target_include_directories( admm SYSTEM
	PUBLIC ${EIGEN_INCLUDE_DIRS}
	PUBLIC ${PYTHON_INCLUDE_DIRS}
	PUBLIC ${HEAD_DIR}
	PUBLIC ${XILINX_HEADERS_DIR}
	PUBLIC ${FLOATX_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_SSE_HEADERS_DIR}
)
target_link_libraries( admm
	PUBLIC m
	PUBLIC fftw3
	PUBLIC pthread
	PUBLIC python2.7
)
# ADMM LASSO
add_executable( admm_lasso examples/c/admm_lasso/admm_lasso.cpp examples/c/admm_lasso/admm_lasso.hpp
		    source/common.cpp header/common.hpp header/data.hpp
     		    ${SOFTPOSIT_SRC}
)
target_include_directories( admm_lasso SYSTEM
	PUBLIC ${EIGEN_INCLUDE_DIRS}
	PUBLIC ${PYTHON_INCLUDE_DIRS}
	PUBLIC ${HEAD_DIR}
	PUBLIC ${XILINX_HEADERS_DIR}
	PUBLIC ${FLOATX_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_HEADERS_DIR}
	PUBLIC ${SOFTPOSIT_SSE_HEADERS_DIR}
)
target_link_libraries( admm_lasso
	PUBLIC m
	PUBLIC fftw3
	PUBLIC pthread
	PUBLIC python2.7
)
