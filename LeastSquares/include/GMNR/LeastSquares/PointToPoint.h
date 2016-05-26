#ifndef GMNR_LEASTSQUARES_POINTTOPOINT_H
#define GMNR_LEASTSQUARES_POINTTOPOINT_H

#ifdef _WIN32
	#ifndef LIBLEASTSQUARES_API
		#ifdef LIBLEASTSQUARES_DYNAMIC
			#if LIBLEASTSQUARES_BUILD
				#define LIBLEASTSQUARES_API __declspec(dllexport)
			#else 
				#define LIBLEASTSQUARES_API __declspec(dllimport)
			#endif
		#else
			#define LIBLEASTSQUARES_API
		#endif
	#endif
#else 
	#define LIBLEASTSQUARES_API
#endif

#include <GMNR/common.h>

namespace gmnr{

	//-----------------------------------------------------------------------------------------------------------------------------------
	//             The Implementation Code Referred To The Following Papers
	//
	//             1. "Least-Squares Fitting of Two 3-D Point Sets" (PAMI 1987)
	//
	//                 K. S. ARUN, T. S. HUANG, AND S. D. BLOSTEIN
	//
	//             2. "Simultaneously Registration of Multiple Corresponding Point Sets" (Computer Vision and Image Understanding 2001)
	//
	//                 John Williams and Mohammed Bennamoun
	//
	//                 Space Center for Satellite Navigation, Queensland University of Technology, GPO BOX 2434, Brisbane 4001, Australia
	//                 E-mail: j2.williams@qut.edu.au; m.bennamoun@qut.edu.au 
	//------------------------------------------------------------------------------------------------------------------------------------

	class PointToPointArun{
	public:
		LIBLEASTSQUARES_API PointToPointArun(const PointSet3D &_X, 
			const PointSet3D &_Y, 
			const bool _calculate_rms_error = false, 
			const Vector &_W = Vector());

		inline Matrix3D rotation(){
			return R_;
		}

		inline Vector3D translation(){
			return t_;
		}

		inline Matrix4D transfomation(){
			Matrix4D tf = Matrix4D::Zero(4, 4);
			tf.block(0,0,3,3) = R_;
			tf.block(0,3,3,1) = t_;
			tf(3, 3) = 1.0;
			return tf;
		}

		//rms_error_ == -1 means the algorithm is faild!
		//rms_error_ == -2 means the RMS error is not calculated!
		inline Scalar rmsError(){
			return rms_error_;
		}

	protected:
		Matrix3D R_;
		Vector3D t_;
		Scalar rms_error_;
	};

	//-----------------------------------------------------------------------------------------------------------------------------------
	//             The Implementation Code Referred To The Following Papers
	//
	//             "Least-Squares Estimation of Transformation Parameters Between Two Point Patterns" (PAMI 1991)
	//
	//                                            Shinji Umeyama
	//
	//------------------------------------------------------------------------------------------------------------------------------------

	class PointToPointUmeyama{
	public:
		LIBLEASTSQUARES_API PointToPointUmeyama(const Matrix &_X, 
			const Matrix &_Y, 
			const bool _with_scaling = false, 
			const bool _calculate_rms_error = false, 
			const Vector &_W = Vector());

		inline Matrix rotation(){
			return R_;
		}

		inline Vector translation(){
			return t_;
		}

		inline Scalar scale(){
			return c_;
		}

		inline Matrix transfomation(){
			int dim = t_.size();
			Matrix tf = Matrix::Zero(dim+1, dim+1);
			tf.block(0,0,dim,dim) = c_*R_;
			tf.block(0,dim,dim,1) = t_;
			tf(dim, dim) = 1.0;
			return tf;
		}

		inline Scalar rmsError(){
			//rms_error_ == -1 means the algorithm is faild!
			//rms_error_ == -2 means the RMS error is not calculated!
			return rms_error_;
		}

	protected:
		Matrix R_;
		Vector t_;
		Scalar c_;
		Scalar rms_error_;
	};
};

#endif