SET(COMMON_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/common/include)
SET(MATH_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/math/include)
SET(TPS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/TPS/include)
SET(DUALTPS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/DualTPS/include)
SET(MULTITPS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/MultiTPS/include)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR}  
	${OpenCV_INCLUDE_DIRS}
	${COMMON_INCLUDE_DIR} 
	${MATH_INCLUDE_DIR} 
	${TPS_INCLUDE_DIR}
	${DUALTPS_INCLUDE_DIR}
	${MULTITPS_INCLUDE_DIR})

IF(GMNR_SHARED_LIBS)
	ADD_EXECUTABLE(test_multitps src/test_multitps.cc)
	IF(CMAKE_HOST_WIN32)
		SET_TARGET_PROPERTIES(test_multitps PROPERTIES COMPILE_FLAGS "-DLIBMATH_DYNAMIC -DLIBTPS_DYNAMIC")
	ENDIF(CMAKE_HOST_WIN32)
	ADD_DEPENDENCIES(test_multitps math tps multitps)
	LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/lib ${OpenCV_LIB_DIR})
	TARGET_LINK_LIBRARIES(test_multitps multitps dualtps tps math ${OpenCV_LIBS})
	INSTALL(TARGETS test_multitps RUNTIME DESTINATION bin)
ELSE(GMNR_SHARED_LIBS)
	ADD_EXECUTABLE(test_multitps_s src/test_multitps.cc)
	ADD_DEPENDENCIES(test_multitps_s multitps_static dualtps_static tps_static math_static)
	IF(CMAKE_HOST_WIN32)
		LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/libstatic ${OpenCV_LIB_DIR})
	ELSEIF(CMAKE_HOST_WIN32)
		LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/lib ${OpenCV_LIB_DIR})
	ENDIF(CMAKE_HOST_WIN32)
	TARGET_LINK_LIBRARIES(test_multitps_s multitps_static dualtps_static tps_static math_static ${OpenCV_LIBS})
	INSTALL(TARGETS test_multitps_s RUNTIME DESTINATION bin)
ENDIF(GMNR_SHARED_LIBS)

SET(test_multitps_qt_HEADERS include/test_multitps_qt_dialog.h include/test_multitps_qt_dialog_control.h)
SET(test_multitps_qt_FORMS include/test_multitps_qt_dialog.ui)
QT4_WRAP_CPP(test_multitps_qt_HEADERS_MOC ${test_multitps_qt_HEADERS})
QT4_WRAP_UI(test_multitps_qt_FORMS_HEADERS ${test_multitps_qt_FORMS})
SET(COMMON_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/common/include)
SET(MULTITPS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/MultiTPS/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF(GMNR_SHARED_LIBS)
	ADD_EXECUTABLE(test_multitps_qt 
					src/test_multitps_qt.cc
					include/test_multitps_qt_dialog.h
					src/test_multitps_qt_dialog.cc
					include/test_multitps_qt_dialog_control.h
					src/test_multitps_qt_dialog_control.cc
					${test_multitps_qt_HEADERS_MOC}
					${test_multitps_qt_FORMS_HEADERS})
	IF(CMAKE_HOST_WIN32)
		SET_TARGET_PROPERTIES(test_multitps_qt 
								PROPERTIES 
									COMPILE_FLAGS 
										"-DLIBMATH_DYNAMIC -DLIBTPS_DYNAMIC -DLIBMULTITPS_DYNAMIC")
	ENDIF(CMAKE_HOST_WIN32)
	ADD_DEPENDENCIES(test_multitps_qt math tps multitps)
	LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/lib)
	TARGET_LINK_LIBRARIES(test_multitps_qt multitps dualtps tps math Qt4::QtCore Qt4::QtGui)
	INSTALL(TARGETS test_multitps_qt RUNTIME DESTINATION bin)
ELSE(GMNR_SHARED_LIBS)
	ADD_EXECUTABLE(test_multitps_qt_s 
					src/test_multitps_qt.cc
					include/test_multitps_qt_dialog.h
					src/test_multitps_qt_dialog.cc
					include/test_multitps_qt_dialog_control.h
					src/test_multitps_qt_dialog_control.cc
					${test_multitps_qt_HEADERS_MOC}
					${test_multitps_qt_FORMS_HEADERS})
	ADD_DEPENDENCIES(test_multitps_qt_s multitps_static dualtps_static tps_static math_static)
	IF(CMAKE_HOST_WIN32)
		LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/libstatic)
	ELSEIF(CMAKE_HOST_WIN32)
		LINK_DIRECTORIES(${PROJECT_BINARY_DIR}/lib)
	ENDIF(CMAKE_HOST_WIN32)
	TARGET_LINK_LIBRARIES(test_multitps_qt_s multitps_static dualtps_static tps_static math_static Qt4::QtCore Qt4::QtGui)
	INSTALL(TARGETS test_multitps_qt_s RUNTIME DESTINATION bin)
ENDIF(GMNR_SHARED_LIBS)