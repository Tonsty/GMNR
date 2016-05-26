#ifndef GMNR_LEASTSQUARES_POINTTOPLANE_H
#define GMNR_LEASTSQUARES_POINTTOPLANE_H

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

	//------------------------------------------------------------------------------------------------------------
	//             The Implementation Code Referred To The Following Paper
	//
	//"Linear Least-Squares Optimization for Point-to-Plane ICP Surface Registration" (Technical Report 2004)
	//
	//                               Kok-Lim Low
	//
	//                      Department of Computer Science
	//                University of North Carnolina at Chapel Hill
	//                          Email: lowk@cs.unc.edu
	//------------------------------------------------------------------------------------------------------------

	class PointToPlaneLinear{
	public:
		LIBLEASTSQUARES_API PointToPlaneLinear(const PointSet3D &_s, 
			const PointSet3D &_d , 
			const NormalSet3D &_n, 
			const bool _calculate_rms_error = false, 
			const Vector &_w = Vector());

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

		inline Scalar rmsError(){
			return rms_error_;
		}

		inline Matrix3D rotationApproximate(){
			return R_approximate_;
		}

		inline Scalar rmsErrorApproximate(){
			return rms_error_approximate_;
		}

	protected:
		Matrix3D R_;
		Vector3D t_;
		Scalar rms_error_;
		
		Matrix3D R_approximate_;
		Scalar rms_error_approximate_;
	};
};

#endif