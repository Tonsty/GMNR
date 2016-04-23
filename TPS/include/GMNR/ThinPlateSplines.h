#ifndef GMNR_TPS_THINPLATESPLINES_H
#define GMNR_TPS_THINPLATESPLINES_H

#ifndef LIBTPS_API
	#ifdef _WIN32
		#ifdef LIBTPS_DYNAMIC
			#if LIBTPS_BUILD
				#define LIBTPS_API __declspec(dllexport)
			#else 
				#define LIBTPS_API __declspec(dllimport)
			#endif
		#else 
			#define LIBTPS_API
		#endif
	#else
		#define LIBTPS_API
	#endif
#endif

#ifdef _DEBUG
#include <iostream>
#endif

#include <GMNR/common.h>
#include <GMNR/math/GreenFunction.h>

namespace gmnr{

	//-----------------------------------------------------------------
	//     The Implementation Code Referred To The Following Paper
	//
	//                  Thin-Plate Splines
	//
	//                    David Eberly
	//
	//                 Geometric Tools, LLC
	//              http://www.geometrictools.com/
	//                Created: March 1, 1996
	//             Last Modified: January 5, 2015
	//-----------------------------------------------------------------

	class TPSFunction{

#ifdef _DEBUG
		// used for debug
		friend LIBTPS_API std::ostream& operator<<(std::ostream& os, const TPSFunction& tps);
#endif

	public:

		LIBTPS_API TPSFunction();

		LIBTPS_API TPSFunction(const Matrix&_X, const Matrix&_A, const Matrix&_B);

		LIBTPS_API TPSFunction(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda = 0.0, const Scalar &_kappa = 0.0);

		LIBTPS_API Matrix evaluate(const Matrix &_X) const;

		inline int variable_Dim() const { return X_.cols();}

		inline int value_Dim() const { return A_.cols(); }

		inline void setX(const Matrix&_X){ X_ = _X;}

		inline void setA(const Matrix&_A){ A_ = _A;}

		inline void setB(const Matrix&_B){ B_ = _B;}

		inline const Matrix& getX() const { return X_;}

		inline const Matrix& getA() const { return A_;}

		inline const Matrix& getB() const { return B_;}

	protected:
		Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> X_;
		Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> A_;
		Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> B_;

	private:
		LIBTPS_API void formula_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa);
		LIBTPS_API void direct_inverse_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa);
		LIBTPS_API void svd_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa);
		LIBTPS_API void qr_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa);

	};
}



#endif