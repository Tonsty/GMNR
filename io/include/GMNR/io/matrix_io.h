#ifndef GMNR_IO_MATRIX_IO_H
#define GMNR_IO_MATRIX_IO_H

#include <GMNR/common.h>

#ifndef LIBIO_API
	#ifdef _WIN32
		#ifdef LIBIO_DYNAMIC
			#if LIBIO_BUILD
				#define LIBIO_API __declspec(dllexport)
			#else 
				#define LIBIO_API __declspec(dllimport)
			#endif
		#else
			#define LIBIO_API
		#endif
	#else
		#define LIBIO_API
	#endif
#endif

namespace gmnr {
	LIBIO_API void write_matrix(const std::string &file_name, const Matrix& mat);
}

#endif
