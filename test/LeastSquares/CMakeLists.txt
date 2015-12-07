SET(test_leastsquares_HEADERS include/test_leastsquares_dialog.h include/test_leastsquares_dialog_control.h)
SET(test_leastsquares_FORMS include/test_leastsquares_dialog.ui)
QT4_WRAP_CPP(test_leastsquares_HEADERS_MOC ${test_leastsquares_HEADERS})
QT4_WRAP_UI(test_leastsquares_FORMS_HEADERS ${test_leastsquares_FORMS})
SET(COMMON_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/common/include)
SET(LEASTSQUARES_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/LeastSquares/include)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR} 
	${CMAKE_CURRENT_BINARY_DIR} 
	${COMMON_INCLUDE_DIR} 
	${LEASTSQUARES_INCLUDE_DIR})

IF(GMNR_SHARED_LIBS)
	ADD_EXECUTABLE(test_leastsquares
					src/test_leastsquares.cc
					include/test_leastsquares_dialog.h
					src/test_leastsquares_dialog.cc
					include/test_leastsquares_dialog_control.h
					src/test_leastsquares_dialog_control.cc
					${test_leastsquares_HEADERS_MOC}
					${test_leastsquares_FORMS_HEADERS})
	IF(CMAKE_HOST_WIN32)
		SET_TARGET_PROPERTIES(test_leastsquares PROPERTIES COMPILE_FLAGS "-DLIBLEASTSQUARES_DYNAMIC")
	ENDIF(CMAKE_HOST_WIN32)
	ADD_DEPENDENCIES(test_leastsquares leastsquares)
	LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/lib)
	TARGET_LINK_LIBRARIES(test_leastsquares leastsquares Qt4::QtCore Qt4::QtGui)
	INSTALL(TARGETS test_leastsquares RUNTIME DESTINATION bin)
ELSE(GMNR_SHARED_LIBS)
	ADD_EXECUTABLE(test_leastsquares_s 
					src/test_leastsquares.cc
					include/test_leastsquares_dialog.h
					src/test_leastsquares_dialog.cc
					include/test_leastsquares_dialog_control.h
					src/test_leastsquares_dialog_control.cc
					${test_leastsquares_HEADERS_MOC}
					${test_leastsquares_FORMS_HEADERS})
	ADD_DEPENDENCIES(test_leastsquares_s leastsquares_static)
	IF(CMAKE_HOST_WIN32)
		LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/libstatic)
	ELSEIF(CMAKE_HOST_WIN32)
		LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/lib)
	ENDIF(CMAKE_HOST_WIN32)
	TARGET_LINK_LIBRARIES(test_leastsquares_s leastsquares_static Qt4::QtCore Qt4::QtGui)
	INSTALL(TARGETS test_leastsquares_s RUNTIME DESTINATION bin)
ENDIF(GMNR_SHARED_LIBS)