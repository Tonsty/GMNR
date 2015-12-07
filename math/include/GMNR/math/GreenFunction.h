#ifndef GMNR_MATH_GREENFUCTION_H
#define GMNR_MATH_GREENFUCTION_H

#ifndef LIBMATH_API
	#ifdef _WIN32
		#ifdef LIBMATH_DYNAMIC
			#if LIBMATH_BUILD
				#define LIBMATH_API __declspec(dllexport)
			#else 
				#define LIBMATH_API __declspec(dllimport)
			#endif
		#else
			#define LIBMATH_API
		#endif
	#else
		#define LIBMATH_API
	#endif
#endif

#include <GMNR/common.h>

namespace gmnr{
	LIBMATH_API Scalar greenFunc( Scalar x, int d);

	LIBMATH_API Matrix greenFunc(const Matrix &_Xi, const Matrix &_Xj);
};

#endif