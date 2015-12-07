#include <iostream>

#include <GMNR/common.h>
#include <GMNR/LeastSquares/PointToPoint.h>

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

	PointToPointArun::PointToPointArun(
		const PointSet3D &_X, 
		const PointSet3D &_Y, 
		const bool _calculate_rms_error, 
		const Vector &_W){

		int N = _X.rows(); //Only be able to deal with 3D data

		const Scalar one_over_N = 1.0 / N;

		Point3D X_mean = _X.colwise().sum() * one_over_N;
		Point3D Y_mean = _Y.colwise().sum() * one_over_N;

		//std::cout << "X_mean = \n" << X_mean << std::endl; 
		//std::cout << "Y_mean = \n" << Y_mean << std::endl; 

		PointSet3D P = _X, Q = _Y;
		for (int i = 0; i < N; i++) P.row(i) -= X_mean.transpose();
		for (int i = 0; i < N; i++) Q.row(i) -= Y_mean.transpose();
		
		Matrix3D H = Matrix3D::Zero();
		for (int i = 0; i < N; i++){
			H += Q.row(i).transpose() * P.row(i);
			//H += ( _Y.row(i).transpose() - Y_mean ) * (_X.row(i) - X_mean.transpose());
		}
		//std::cout << "H = \n" << H << std::endl; 

		Eigen::JacobiSVD<Matrix> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
		R_ = svd.matrixV() * svd.matrixU().transpose();

		//std::cout << "R.determinant() = \n" << R_.determinant() << std::endl; 
		if (R_.determinant() < 0 ){ //det(R) == -1
			if ( svd.nonzeroSingularValues() < 3) //One of the singular values (lambda3, say) of H is zero
			{
				Matrix3D V_new = svd.matrixV();
				V_new.col(2) = V_new.col(2) * -1;
				R_ = V_new * svd.matrixU().transpose();
			}
			else{//None of the singular values of H is zero
				std::cerr << "class \"PointToPointArun\" constructor failed!!!" << std::endl;
				rms_error_ = -1;
				return;
			}
		}
		t_ = Y_mean - R_ * X_mean;

		if (_calculate_rms_error){
			rms_error_ = (Q - P * R_.transpose()).squaredNorm();
			rms_error_ = sqrt((double)rms_error_ * one_over_N);
		}else{
			rms_error_ = -2; // uncalculated
		}
		std::cout << "rms_error = \n" << rms_error_ << std::endl; 
	}

	//-----------------------------------------------------------------------------------------------------------------------------------
	//             The Implementation Code Referred To The Following Papers
	//
	//             "Least-Squares Estimation of Transformation Parameters Between Two Point Patterns" (PAMI 1991)
	//
	//                                            Shinji Umeyama
	//
	//------------------------------------------------------------------------------------------------------------------------------------

	PointToPointUmeyama::PointToPointUmeyama(
		const Matrix &_X, 
		const Matrix &_Y, 
		const bool _with_scaling, 
		const bool _calculate_rms_error, 
		const Vector &_W){

		const int n = _X.rows();
		const int m = _X.cols(); //Be able to deal with any m-dimensional data.

		const Scalar one_over_n = 1.0 / n;

		const Vector X_mean = _X.colwise().sum() * one_over_n;
		const Vector Y_mean = _Y.colwise().sum() * one_over_n;

		//std::cout << "X_mean = \n" << X_mean << std::endl; 
		//std::cout << "Y_mean = \n" << Y_mean << std::endl; 

		Matrix X_demean = _X;
		Matrix Y_demean = _Y;
		for (int i = 0; i < n; i++) X_demean.row(i) -= X_mean.transpose();
		for (int i = 0; i < n; i++) Y_demean.row(i) -= Y_mean.transpose();

		Scalar X_delta2 = 0;
		Scalar Y_delta2 = 0;
		for (int i = 0; i < n; i++) X_delta2 += X_demean.row(i).squaredNorm();
		for (int i = 0; i < n; i++) Y_delta2 += Y_demean.row(i).squaredNorm();
		X_delta2 *= one_over_n;
		Y_delta2 *= one_over_n;

		Matrix XY_covariance = Matrix::Zero(m, m);
		for (int i = 0; i < n; i++){
			XY_covariance += (_Y.row(i).transpose() - Y_mean) * (_X.row(i) - X_mean.transpose());
		}
		XY_covariance *= one_over_n;
		//std::cout << "XY_covariance = \n" << XY_covariance << std::endl;

		Eigen::JacobiSVD<Matrix> svd(XY_covariance, Eigen::ComputeThinU | Eigen::ComputeThinV);
		//std::cout << "UDV^t = \n" << svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().transpose() << std::endl;

		Vector D = svd.singularValues();
		Vector S = Vector::Ones(m);
		//std::cout << "S =\n" << S << std::endl;
		//std::cout << "D =\n" << D << std::endl;

		if (XY_covariance.determinant() <0) S(m-1) = -1;

		//nsigned int svd_rank = svd.rank();
		unsigned int svd_rank = 0; 
		for (unsigned int i=0; i<m; ++i) if (!Eigen::internal::isMuchSmallerThan(D.coeff(i),D.coeff(0))) ++svd_rank;
		//std::cout << "svd_rank = \n" << svd_rank << std::endl;
		if (svd_rank == m-1){
			if ( svd.matrixU().determinant() * svd.matrixV().determinant() > 0){ //det(U) * det(V) = 1 
				R_.noalias() = svd.matrixU() * svd.matrixV().transpose();
			}else{ //det(U) * det(V) = -1 
				const Scalar s = S(m-1); S(m-1) = -1;
				R_.noalias() = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
				S(m-1) = s;
			}
		}else if(svd_rank == m){ // svd_rank > m-1
			R_.noalias() = svd.matrixU() * S.asDiagonal() * svd.matrixV().transpose();
		}else{ // svd_rank < m-1
			std::cerr << "class \"PointToPointUmeyama\" constructor failed!!!" << std::endl;
			rms_error_ = -1;
			return;
		}

		if (_with_scaling){
			c_ = 1.0/X_delta2 * D.dot(S);
			t_ = Y_mean;
			t_.noalias() -= c_ * R_ * X_mean; 
				
			if(_calculate_rms_error) rms_error_ = sqrt(Y_delta2 - D.dot(S) * D.dot(S) / X_delta2);
			else rms_error_ = -2; // uncalculated
		}else{
			c_ = 1;
			t_ = Y_mean;
			t_.noalias() -= R_ * X_mean;

			if(_calculate_rms_error) rms_error_ = sqrt((X_demean - Y_demean * R_.transpose()).squaredNorm() * one_over_n);
			else rms_error_ = -2; // uncalculated
		}

		std::cout << "rms_error = \n" << rms_error_ << std::endl; 
	}

};
