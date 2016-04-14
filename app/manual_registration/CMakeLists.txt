INCLUDE_DIRECTORIES(${PCL_INCLUDE_DIRS})
LINK_DIRECTORIES(${PCL_LIBRARY_DIRS} ${Boost_LIBRARY_DIR_DEBUG} ${Boost_LIBRARY_DIR_RELEASE} ${VTK_DIR})

SET(manual_registration_qt_HEADERS manual_registration.h)
SET(manual_registration_qt_FORMS manual_registration.ui)

QT4_WRAP_CPP(manual_registration_qt_HEADERS_MOC ${manual_registration_qt_HEADERS})
QT4_WRAP_UI(manual_registration_qt_FORMS_HEADERS ${manual_registration_qt_FORMS})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_EXECUTABLE(manual_registration 
					main.cpp
					manual_registration.cpp
					manual_registration.h
					${manual_registration_qt_HEADERS_MOC}
					${manual_registration_qt_FORMS_HEADERS})
TARGET_LINK_LIBRARIES(manual_registration Qt4::QtCore Qt4::QtGui ${VTK_LIBRARIES} ${PCL_LIBRARIES} debug QVTK_debug optimized QVTK_release)