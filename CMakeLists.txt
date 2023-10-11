  
OPTION(USE_MKL "MKL" OFF)#TODO: some fixes is needed
OPTION(USE_OPENBLAS "OpenBLAS" OFF)
OPTION(COMPILE_R "R" ON)

# BUILD_COMPLETE : If a package is not found, it throws error. Otherwise, it builds the project but without the required package

if (BUILD_COMPLETE)
  message(STATUS "BUILD_COMPLETE is ON. If a required package is not found, it throws error.")
else()
  message(WARNING "BUILD_COMPLETE is OFF. If a required package is not found, it does not stop the build.")
endif()

cmake_minimum_required(VERSION 3.12)
set(CXX_STD_VERSION cxx_std_17)

message(STATUS "MINGW_ARCH: ${MINGW_ARCH}")
if (CMAKEPRESETS) 
   set(PACKAGES_ROOT C:/msys64/${MINGW_ARCH})
   message(STATUS "Package roots: ${PACKAGES_ROOT}")
else()
   message(STATUS "Package roots: Unavailable")
endif()
  
set(LDT_VERSION_MAJOR 0)
set(LDT_VERSION_MINOR 1)
set(LDT_DESCRIPTION "LDT v. ${LDT_VERSION_MAJOR}.${LDT_VERSION_MINOR} by Ramin Mojab (rmojab63@gmail.com)")
 
project(ldt VERSION ${LDT_VERSION_MAJOR}.${LDT_VERSION_MINOR} DESCRIPTION "${LDT_DESCRIPTION}")
 
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules")
 
message(STATUS "Project source directory: ${PROJECT_SOURCE_DIR}")

file(GLOB SOURCES
 src/api.cpp  
 
 #src/ldt_statistics/modeling/time_series/expsmoothing.cpp
 #src/ldt_statistics/modeling/time_series/tsfilter.cpp 
 
 src/ldt_statistics/modeling/linear/sur.cpp
 src/ldt_statistics/modeling/linear/sur_projection.cpp
 src/ldt_statistics/modeling/linear/sur_extended.cpp
 src/ldt_statistics/modeling/linear/sur_simulation.cpp
 src/ldt_statistics/modeling/linear/sur_modelset.cpp

 src/ldt_statistics/modeling/time_series/varma/varma_sizes.cpp 
 src/ldt_statistics/modeling/time_series/varma/varma_arma.cpp
 src/ldt_statistics/modeling/time_series/varma/varma.cpp
 src/ldt_statistics/modeling/time_series/varma/varma_forecast.cpp
 src/ldt_statistics/modeling/time_series/varma/varma_restriction.cpp
 src/ldt_statistics/modeling/time_series/varma/varma_extended.cpp
 src/ldt_statistics/modeling/time_series/varma/varma_simulation.cpp
 src/ldt_statistics/modeling/time_series/varma/varma_modelset.cpp

 src/ldt_statistics/modeling/discrete/discrete_choice.cpp
 src/ldt_statistics/modeling/discrete/discrete_choice_sim.cpp
 src/ldt_statistics/modeling/discrete/discrete_choice_modelset.cpp
 src/ldt_statistics/modeling/discrete/discrete_choice_extended.cpp 

 src/ldt_statistics/modeling/search/modelset.cpp
 src/ldt_statistics/modeling/search/searcher.cpp
 src/ldt_statistics/modeling/search/searcher_reg.cpp
 src/ldt_statistics/modeling/search/summary.cpp

 src/ldt_statistics/clustering/HCluster.cpp
 src/ldt_statistics/clustering/distance.cpp
 src/ldt_statistics/clustering/group_data.cpp

  src/ldt_statistics/scoring/scoring.cpp
  src/ldt_statistics/scoring/frequency_cost.cpp
  src/ldt_statistics/scoring/roc.cpp

 src/ldt_statistics/correlation.cpp
 src/ldt_statistics/pca.cpp
 src/ldt_statistics/data_split.cpp

 src/ldt_statistics/distributions/distributions_uni.cpp
 src/ldt_statistics/distributions/distribution_gld.cpp
 src/ldt_statistics/distributions/histogram.cpp
 src/ldt_statistics/distributions/mixture.cpp
 src/ldt_statistics/distributions/mnormal.cpp 
 src/ldt_statistics/distributions/empirical103.cpp 
 
 src/ldt_statistics/statistics/statistics.cpp 
 src/ldt_statistics/statistics/descriptive.cpp
 

 src/ldt_math/optimization/lmbfgsb.cpp
 src/ldt_math/optimization/nelder_mead.cpp
 src/ldt_math/optimization/newton.cpp

 src/ldt_math/optimization/lbfgsb/lbfgsb.c
 src/ldt_math/optimization/lbfgsb/linesearch.c
 src/ldt_math/optimization/lbfgsb/linpack.c
 src/ldt_math/optimization/lbfgsb/miniCBLAS.c
 src/ldt_math/optimization/lbfgsb/print.c
 src/ldt_math/optimization/lbfgsb/subalgorithms.c
 src/ldt_math/optimization/lbfgsb/timer.c
 
 src/ldt_math/polynomials/mpolynomial.cpp
 src/ldt_math/polynomials/polynomial.cpp
 src/ldt_math/derivative.cpp
 src/ldt_math/auc.cpp

 src/matrix/matrix_la.cpp
 src/matrix/matrix.cpp
 src/matrix/matrix_sym.cpp
 src/matrix/dataset.cpp 
 src/matrix/decompositions.cpp
 src/matrix/vmatrix.cpp

 src/variable/variables.cpp
 src/variable/variable.cpp
 
 src/frequency/f_base.cpp
 src/frequency/f_cross_section.cpp
 src/frequency/f_year_based.cpp
 src/frequency/f_week_based.cpp
 src/frequency/f_list.cpp
 src/frequency/f_day_based.cpp

 src/helpers.cpp
 src/array.cpp
  
)

SET(LDT_HEADER_DIR ${PROJECT_SOURCE_DIR}/include 
                   ${PROJECT_SOURCE_DIR}/include/math 
				   ${PROJECT_SOURCE_DIR}/include/blas 
				   ${PROJECT_SOURCE_DIR}/include/math/statistics
				   ${PROJECT_SOURCE_DIR}/include/math/statistics/inference
				   ${PROJECT_SOURCE_DIR}/include/math/statistics/inference/models
				   )

if(COMPILE_R)
    list(APPEND SOURCES
         Rpackage/src/r_ldt.h
         Rpackage/src/r_statistics.cpp
		 Rpackage/src/r_clustering.cpp
		 Rpackage/src/r_ldt.cpp
		 Rpackage/src/r_sur.cpp
		 Rpackage/src/r_varma.cpp
		 Rpackage/src/r_dc.cpp
		 Rpackage/src/r_gen.cpp
    )
  
	# currently, I hard code the paths:  
	set(R_PATH "C:/Program Files/R/R-4.3.1")
	set(R_CPP_PATH "C:/Users/rmoja/AppData/Local/R/win-library/4.3/Rcpp")
 
  include_directories(${R_PATH}/include)
  include_directories(${R_CPP_PATH}/include)
 
  message(STATUS "R Include Dir: ${R_INCLUDE_DIR}/include")
  message(STATUS "Rcpp Include Dir: ${R_CPP_PATH}/include")
  
endif()

# ------ Rcpp: Frequency package

file(GLOB SOURCES_F 
 include/ldt_base.h
 include/array.h
 include/helpers.h
 include/frequency.h
 include/variable.h
 src/helpers.cpp
 src/array.cpp
 src/variable/variables.cpp
 src/variable/variable.cpp
 src/frequency/f_base.cpp
 src/frequency/f_cross_section.cpp
 src/frequency/f_year_based.cpp
 src/frequency/f_week_based.cpp
 src/frequency/f_list.cpp
 src/frequency/f_day_based.cpp
)
file(COPY ${SOURCES_F} DESTINATION ${PROJECT_SOURCE_DIR}/Rpackage_F/src)
file(COPY ${SOURCES_F} DESTINATION ${PROJECT_SOURCE_DIR}/Rpackage_F/src)
message(STATUS "Source files and headers are copied to Rpackage_F/src directory.")

# ------

# ------ Rcpp: Copy files to R/src directory (They are ignored by git)
file(COPY ${SOURCES} DESTINATION ${PROJECT_SOURCE_DIR}/Rpackage/src)
foreach (hd ${LDT_HEADER_DIR})
file(GLOB HD "${hd}/*.h" )
file(COPY ${HD} DESTINATION ${PROJECT_SOURCE_DIR}/Rpackage/src)
endforeach()
file(GLOB REMOVE_FILES  ${PROJECT_SOURCE_DIR}/Rpackage/src/api.h
                        ${PROJECT_SOURCE_DIR}/Rpackage/src/api.cpp 
						${PROJECT_SOURCE_DIR}/Rpackage/src/frequency.h
						${PROJECT_SOURCE_DIR}/Rpackage/src/variable.h 
						${PROJECT_SOURCE_DIR}/Rpackage/src/variables.cpp
						${PROJECT_SOURCE_DIR}/Rpackage/src/variable.cpp
						${PROJECT_SOURCE_DIR}/Rpackage/src/f_base.cpp
						${PROJECT_SOURCE_DIR}/Rpackage/src/f_cross_section.cpp
						${PROJECT_SOURCE_DIR}/Rpackage/src/f_year_based.cpp
						${PROJECT_SOURCE_DIR}/Rpackage/src/f_week_based.cpp
						${PROJECT_SOURCE_DIR}/Rpackage/src/f_list.cpp
						${PROJECT_SOURCE_DIR}/Rpackage/src/f_day_based.cpp
						
						${PROJECT_SOURCE_DIR}/Rpackage/src/empirical103.cpp
						${PROJECT_SOURCE_DIR}/Rpackage/src/histogram.cpp)
file(REMOVE ${REMOVE_FILES})
message(STATUS "Source files and headers are copied to Rpackage/src directory.")
# Copy copyright file to R package
file(GLOB COPYRIGHT  ${PROJECT_SOURCE_DIR}/readme.md)
file(COPY ${COPYRIGHT} DESTINATION ${PROJECT_SOURCE_DIR}/Rpackage/inst/COPYRIGHTS)
# ------

add_library(ldt SHARED ${SOURCES})

message(STATUS "Header Directories: ${LDT_HEADER_DIR}") 
include_directories(${LDT_HEADER_DIR})
	 
target_compile_features(ldt PUBLIC ${CXX_STD_VERSION})
 
if(MSVC)
	target_compile_options(ldt PRIVATE "/wd4251;/wd4661")
endif()
if ( CMAKE_COMPILER_IS_GNUCC )
    #target_compile_options(ldt PRIVATE -DUSE_TYPE_CHECKING_STRICT -D_FORTIFY_SOURCE=2 -fpic  -g -O3 -Wall -pedantic)
endif()
 
  if (CMAKEPRESETS)
    set(BOOST_ROOT ${PACKAGES_ROOT})
    set(Boost_NO_SYSTEM_PATHS ON)
  endif()

  set(Boost_DETAILED_FAILURE_MSG ON)
  find_package(Boost)
  if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    target_link_libraries(ldt ${Boost_LIBRARIES})
    message(STATUS "BOOST Include Dir: ${Boost_INCLUDE_DIRS}")
	ADD_DEFINITIONS(-DUSE_BOOST)
  elseif(BUILD_COMPLETE)
    message(FATAL_ERROR "Boost is not found.")
  else()
    message(WARNING "Boost is not found.")
  endif()


if(USE_MKL) #MKL
	set(MKL_ARCH intel64) #ia32 intel64
	set(MKL_INTERFACE_FULL intel_ilp64) #intel_ilp64 intel_lp64
	set(MKL_LINK static) #static dynamic sdl
	set(MKL_THREADING intel_thread) #sequential intel_thread tbb_thread
	set(MKL_MPI intelmpi) #intelmpi mshpc msmpi
	find_package(MKL CONFIG REQUIRED) 

	include_directories(${MKL_ROOT}/include) 
	message(STATUS "Imported =  ${MKL_ROOT}")
	target_link_libraries(ldt ${MKL_IMPORTED_TARGETS}) 

	 ADD_DEFINITIONS(-DUSE_MKL)
elseif(USE_OPENBLAS)
   
    set(BLA_STATIC ON)
	set(BLA_VENDOR OpenBLAS)

	if (CMAKEPRESETS)
	  set(BLAS_ROOT ${PACKAGES_ROOT}) 
	  set(BLAS_NO_SYSTEM_PATHS ON)
	endif()

	if (USE_INTEGER64)
       #TODO
	endif(USE_INTEGER64)
    find_package(BLAS)
	if (BLAS_FOUND)
	  target_link_libraries(ldt ${BLAS_LIBRARIES}) 
	  message(STATUS  "BLAS STATUS: ${BLAS_LIBRARIES}")
	  find_path(BLAS_INCLUDE_DIRS cblas.h PATHS CMAKE_SYSTEM_PREFIX_PATH PATH_SUFFIXES include/openblas)
	  message(STATUS  "BLAS Include Directory: ${BLAS_INCLUDE_DIRS}") 
	  include_directories(${BLAS_INCLUDE_DIRS})
    elseif(BUILD_COMPLETE)
	  message(FATAL_ERROR "BLAS is not found.")
	else()
	  message(WARNING "BLAS is not found.")
	endif()
	ADD_DEFINITIONS(-DUSE_OPENBLAS)

else()

set(BLA_STATIC ON) 
if (CMAKEPRESETS)
	  set(BLAS_ROOT ${PACKAGES_ROOT}) 
	  set(BLAS_NO_SYSTEM_PATHS ON)
	endif()
    find_package(BLAS)
	if (BLAS_FOUND)
	  target_link_libraries(ldt ${BLAS_LIBRARIES}) 
	  message(STATUS  "BLAS STATUS: ${BLAS_LIBRARIES}")
	  
    elseif(BUILD_COMPLETE)
	  message(FATAL_ERROR "BLAS is not found.")
	else()
	  message(WARNING "BLAS is not found.")
	endif()

	find_package(LAPACK)
	if (LAPACK_FOUND)
	  target_link_libraries(ldt ${LAPACK_LIBRARIES}) 
	  message(STATUS  "LAPACK STATUS: ${LAPACK_LIBRARIES}")
	  
    elseif(BUILD_COMPLETE)
	  message(FATAL_ERROR "LAPACK is not found.")
	else()
	  message(WARNING "LAPACK is not found.")
	endif()

endif () 


if(COMPILE_R) 
  target_link_libraries(ldt "${R_PATH}/bin/x64/R.dll")
  target_link_libraries(ldt "${R_CPP_PATH}/libs/x64/Rcpp.dll")
endif()

if(BUILD_TEST)
     message(STATUS "Building tests...")
	  
     file(GLOB TEST_SOURCES 
	 test/distributions_t.cpp 
	 test/functionx_t.cpp  
	 test/matrix_t.cpp  
	 test/optimization_t.cpp  
	 test/polynomial_t.cpp 
	 test/running_t.cpp 
	 test/scoring_t.cpp
	 test/statistics_t.cpp 
	 test/correlation_t.cpp 
	 test/distance_t.cpp 
	 test/varma_t.cpp  
	 #test/sur_t.cpp  
	 test/cluster_t.cpp  
	 #test/discrete_choice_t.cpp
	 #test/searcher_t.cpp
	 test/pca_t.cpp  
	 test/data_split_t.cpp 
	 test/frequency_t.cpp 
	 test/variable_t.cpp 
	 test/array_t.cpp 
   ) 
  
	 enable_testing()
	 set(GTest_ROOT ${PACKAGES_ROOT})
	 find_package(GTest CONFIG REQUIRED)

     add_executable(testsldt ${TEST_SOURCES})
     target_compile_features(testsldt PUBLIC ${CXX_STD_VERSION})
 
	 target_link_libraries(testsldt ldt GTest::gtest GTest::gtest_main)
	 
	 #include_directories(${GTEST_INCLUDE_DIRS})
	 #message(STATUS "GTest Include Directory: ${GTEST_INCLUDE_DIRS}")

     if (USE_MKL)
        target_link_libraries(testsldt ${MKL_IMPORTED_TARGETS})
        #set_tests_properties(AppTest PROPERTIES ENVIRONMENT "${MKL_ENV}")
     else()
        target_link_libraries(testsldt ${BLAS_LIBRARIES}) 
     endif(USE_MKL) 
 	 
	 add_test(TESTS testsldt)

endif(BUILD_TEST)


install(TARGETS ldt
		RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
		LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
		ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY ${LDT_HEADER_DIR}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/ldt)
 
 
include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md")
set(CPACK_PACKAGE_VERSION_MAJOR "${LDT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${LDT_VERSION_MINOR}")
set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)
 