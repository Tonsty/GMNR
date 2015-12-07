#ifndef GMNR_DUALTPS_DUALTPS_H
#define GMNR_DUALTPS_DUALTPS_H

#ifdef _WIN32
	#ifndef LIBDUALTPS_API
		#ifdef LIBDUALTPS_DYNAMIC
			#if LIBDUALTPS_BUILD
				#define LIBDUALTPS_API __declspec(dllexport)
			#else 
				#define LIBDUALTPS_API __declspec(dllimport)
			#endif
		#else
			#define LIBDUALTPS_API
		#endif
	#endif
#else
	#define LIBDUALTPS_API
#endif

#ifdef _DEBUG
#include <iostream>
#endif

#include <GMNR/ThinPlateSplines.h>

namespace gmnr{

	class DualTPS{

#ifdef _DEBUG
		// used for debug
		friend LIBDUALTPS_API std::ostream& operator<<(std::ostream& os, const DualTPS& tps);
#endif

	public:
		LIBDUALTPS_API DualTPS();
		
		LIBDUALTPS_API DualTPS(const Matrix &_X, 
			const Matrix &_Y, 
			const Scalar &_kappa1, 
			const Scalar &_kappa2, 
			const Scalar &_lambda1, 
			const Scalar &_lambda2);

		inline const TPSFunction& getf1() const { return f1_;}

		inline const TPSFunction& getf2() const { return f2_;}

	protected:
		TPSFunction f1_, f2_;

	private:
		LIBTPS_API void formula_solve_(const Matrix &_X, 
			const Matrix &_Y, 
			const Scalar &_kappa1, 
			const Scalar &_kappa2, 
			const Scalar &_lambda1, 
			const Scalar &_lambda2);
		LIBTPS_API void direct_inverse_solve_(const Matrix &_X, 
			const Matrix &_Y, 
			const Scalar &_kappa1, 
			const Scalar &_kappa2, 
			const Scalar &_lambda1, 
			const Scalar &_lambda2);
		LIBTPS_API void svd_solve_(const Matrix &_X, 
			const Matrix &_Y, 
			const Scalar &_kappa1, 
			const Scalar &_kappa2, 
			const Scalar &_lambda1, 
			const Scalar &_lambda2);
		LIBTPS_API void qr_solve_(const Matrix &_X, 
			const Matrix &_Y, 
			const Scalar &_kappa1, 
			const Scalar &_kappa2, 
			const Scalar &_lambda1, 
			const Scalar &_lambda2);
	};
}

#endif