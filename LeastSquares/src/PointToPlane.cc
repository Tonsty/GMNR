#include <iostream>

#include <GMNR/common.h>
#include <GMNR/LeastSquares/PointToPlane.h>

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


	PointToPlaneLinear::PointToPlaneLinear(const PointSet3D &_s, 
		const PointSet3D &_d , 
		const NormalSet3D &_n, 
		const bool _calculate_rms_error, 
		const Vector &_w){

		//std::cout << "_s = \n" << _s << std::endl; 
		//std::cout << "_d = \n" << _d << std::endl; 
		//std::cout << "_n = \n" << _n << std::endl; 
			
		int N = _s.rows();
		Matrix A(N, 6);
		Vector b(N), x(6);

		for (int i = 0; i < N; i++){
			A(i, 0) = _n(i, 2) * _s(i, 1) - _n(i, 1) * _s(i, 2);
			A(i, 1) = _n(i, 0) * _s(i, 2) - _n(i, 2) * _s(i, 0);
			A(i, 2) = _n(i, 1) * _s(i, 0) - _n(i, 0) * _s(i, 1);
		}
		A.block(0, 3, N, 3) = _n;

		for (int i = 0; i < N; i++){
			b(i, 0) = _n.row(i).dot(_d.row(i)) - _n.row(i).dot(_s.row(i));
		}
		A.block(0, 3, N, 3) = _n; 

		Eigen::JacobiSVD<Matrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		x = svd.solve(b);

		Scalar alpha = x(0), beta = x(1), gamma = x(2);
		Quaternion R = AngleAxis(gamma, Vector3D::UnitZ()) * AngleAxis(beta, Vector3D::UnitY()) * AngleAxis(alpha, Vector3D::UnitX());
		R_ = R.matrix();

		//std::cout << "R = \n" << R_ << std::endl;

		//R_(0,0) = cos(gamma) * cos(beta);
		//R_(0,1) = -sin(gamma) * cos(alpha) + cos(gamma) * sin(beta) * sin(alpha);
		//R_(0,2) = sin(gamma) * sin(alpha) + cos(gamma) * sin(beta) * cos(alpha);

		//R_(1,0) = sin(gamma) * cos(beta);
		//R_(1,1) = cos(gamma) * cos(alpha) + sin(gamma) * sin(beta) * sin(alpha);
		//R_(1,2) = -cos(gamma) * sin(alpha) + sin(gamma) * sin(beta) * cos(alpha);

		//R_(2,0) = -sin(beta);
		//R_(2,1) = cos(beta) * sin(alpha);
		//R_(2,2) = cos(beta) * cos(alpha);

		//std::cout << "R = \n" << R_ << std::endl;

		t_ = x.block<3,1>(3,0);

		R_approximate_ << 1, -gamma, beta, 
			              gamma, 1, -alpha,
						  -beta, alpha, 1;

		if (_calculate_rms_error){
			rms_error_ = 0;
			for(int i = 0; i < N; i++){
				rms_error_ += (R_ * _s.row(i).transpose() + t_ - _d.row(i).transpose()).dot(_n.row(i));
			}
			rms_error_ = sqrt((double)rms_error_/N);

			rms_error_approximate_ = sqrt((A * x - b).squaredNorm()/N);
		}else {
			rms_error_ = -2; // uncalculated
			rms_error_approximate_ = -2; // uncalculated
		}

		std::cout << "rms_error = \n" << rms_error_ << std::endl;
		std::cout << "rms_error_approximate = \n" << rms_error_approximate_ << std::endl; 
	}
};