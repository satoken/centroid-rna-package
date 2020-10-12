add_subdirectory(contrafold)
add_subdirectory(contralign)
add_subdirectory(probconsRNA)

add_library(centroid STATIC
	folding_engine.cpp centroid.cpp 
	mea.cpp bp.cpp 
	centroid_fold_interface.cpp 
	fa.cpp aln.cpp 
	engine/aux.cpp 
	engine/alifold.cpp 
	engine/contrafold.cpp 
	engine/mccaskill.cpp 
	engine/averaged.cpp 
	engine/contrafoldm.cpp 
	engine/mixture.cpp 
	engine/pfold.cpp 
	engine/contrafoldhom.cpp 
	engine/mccaskillhom.cpp 
    diana.cpp ps_plot.cpp naview.c
    engine/boltzmann_param.c )

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${VIENNARNA_INCLUDE_DIRS} 
    ${Boost_INCLUDE_DIRS})

add_executable(centroid_fold centroid_fold.cpp)
target_link_libraries(centroid_fold
    centroid contrafold
    ${VIENNARNA_LDFLAGS}
    Boost::program_options Boost::random Boost::system)

add_executable(centroid_alifold centroid_alifold.cpp)
target_link_libraries(centroid_alifold
    centroid contrafold
    ${VIENNARNA_LDFLAGS}
    Boost::program_options Boost::random Boost::system)

add_executable(centroid_homfold centroid_homfold.cpp)
target_link_libraries(centroid_homfold
    centroid contrafold contralign probcons
    ${VIENNARNA_LDFLAGS}
    Boost::program_options Boost::random Boost::system)

install(TARGETS centroid_fold RUNTIME DESTINATION bin)
install(TARGETS centroid_alifold RUNTIME DESTINATION bin)
install(TARGETS centroid_homfold RUNTIME DESTINATION bin)