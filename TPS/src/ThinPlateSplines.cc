#include <iostream>

#include <GMNR/common.h>
#include <GMNR/math/GreenFunction.h>
#include <GMNR/ThinPlateSplines.h>

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


	TPSFunction::TPSFunction(){}

	TPSFunction::TPSFunction(const Matrix&_X, const Matrix&_A, const Matrix&_B): X_(_X), A_(_A), B_(_B){}

	TPSFunction::TPSFunction(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa) {
		qr_solve_(_X, _Y, _lambda, _kappa);
		//formula_solve_(_X, _Y, _lambda, _kappa);
		//svd_solve_(_X, _Y, _lambda, _kappa);
		//direct_inverse_solve_(_X, _Y, _lambda, _kappa);
	}

	void TPSFunction::formula_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa){
		Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
		Matrix M = (1 + _kappa) * greenFunc(_X, _X) + _lambda * Matrix::Identity(_X.rows(), _X.rows());

		Eigen::JacobiSVD<Matrix> svd_1(M);
		std::cout << "formula solver Condition Number 1: " << svd_1.singularValues()[0] / svd_1.singularValues().reverse()[0] << std::endl;

		Matrix inverseM = M.inverse();
		Matrix transposehomoX_inverseM = X.transpose() * inverseM;

		Eigen::JacobiSVD<Matrix> svd_2(transposehomoX_inverseM * (1+_kappa) * X);
		std::cout << "formula solver Condition Number 2: " << svd_2.singularValues()[0] / svd_2.singularValues().reverse()[0] << std::endl;

		B_ = (transposehomoX_inverseM * (1+_kappa) * X).inverse() * transposehomoX_inverseM * (_Y + _kappa * _X);
		A_ = inverseM * (_Y + _kappa * _X - (1+_kappa) * X * B_);
		X_ = _X;

		int m = _X.rows(), d = _X.cols();

		std::cout << "formula solver TPS energy : " << _lambda * (A_.transpose() * (M + _lambda * Matrix::Identity(m, m)) * A_).trace() << std::endl;	
	}

	void TPSFunction::direct_inverse_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa){
		int m = _X.rows(), d = _X.cols();
		
		Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
		Matrix M = (1 + _kappa) * greenFunc(_X, _X) + _lambda * Matrix::Identity(m, m);

		Matrix equationL(m + (d+1), m + (d+1)), equationR(m + (d+1), d), solution(m + (d+1), d);
		equationL << M, (1 + _kappa) * X,
			X.transpose(), Matrix::Zero(d+1, d+1);
		equationR << _Y + _kappa * _X, Matrix::Zero(d+1, d);

		Eigen::JacobiSVD<Matrix> svd(equationL);
		std::cout << "direct inverse solver Condition Number: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;

		solution = equationL.inverse() * equationR;

		A_ = solution.block(0, 0, m, d);
		B_ = solution.block(m, 0, d+1, d);
		X_ = _X;

		std::cout << "direct inverse solver TPS energy : " << _lambda * (A_.transpose() * (M + _lambda * Matrix::Identity(m, m)) * A_).trace() << std::endl;
	}

	void TPSFunction::svd_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa){
		int m = _X.rows(), d = _X.cols();

		Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
		Matrix M = (1 + _kappa) * greenFunc(_X, _X) + _lambda * Matrix::Identity(m, m);

		Matrix equationL(m + (d+1), m + (d+1)), equationR(m + (d+1), d), solution(m + (d+1), d);
		equationL << M, (1 + _kappa) * X,
			X.transpose(), Matrix::Zero(d+1, d+1);
		equationR << _Y + _kappa * _X, Matrix::Zero(d+1, d);

		Eigen::JacobiSVD<Matrix> svd(equationL, Eigen::ComputeThinU | Eigen::ComputeThinV);
		std::cout << "svd solver Condition Number: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
		//std::cout << "Rows = " << equationL.rows() << " Columns = " << equationL.cols() << std::endl;
		//std::cout << "before reset threshold, svd_rank = " << svd.rank() << std::endl;
		//svd.setThreshold(Eigen::NumTraits<Scalar>::epsilon());
		//std::cout << "after reset threshold, svd_rank = " << svd.rank() << std::endl;
		solution = svd.solve(equationR);

		A_ = solution.block(0, 0, m, d);
		B_ = solution.block(m, 0, d+1, d);
		X_ = _X;

		std::cout << "svd solver TPS energy : " << _lambda * (A_.transpose() * (M + _lambda * Matrix::Identity(m, m)) * A_).trace() << std::endl;
	}

	void TPSFunction::qr_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa){
		int m = _X.rows(), d = _X.cols();

		Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
		Matrix M = greenFunc(_X, _X);

		Eigen::ColPivHouseholderQR<Matrix> qr(X);
		//std::cout << "Rows = " << X.rows() << " Columns = " << X.cols() << std::endl;
		//std::cout << "qr_rank = " << qr.rank() << std::endl;
		Matrix P = qr.colsPermutation();
		//std::cout << "P = \n" << P << std::endl;
		Matrix Q = qr.matrixQ();
		//std::cout << "Q = \n" << Q << std::endl;
		Matrix Q1 = Q.block(0, 0, m, d+1), Q2 = Q.block(0, d+1, m, m-(d+1));
		Matrix R = qr.matrixR().topLeftCorner(d+1, d+1).template triangularView<Eigen::Upper>();
		//Matrix R = qr.matrixR().topLeftCorner(qr.rank(), qr.rank()).template triangularView<Eigen::Upper>();
		//std::cout << "R = \n" << R << std::endl;

		Eigen::JacobiSVD<Matrix> svd_formula( Q2.transpose() * (1 + _kappa) * M * Q2 + _lambda * Matrix::Identity(m - (d+1), m - (d+1)) );
		std::cout << "qr solver formula subsolver Condition Number: " << svd_formula.singularValues()[0] / svd_formula.singularValues().reverse()[0] << std::endl;

		A_ = Q2 * (Q2.transpose() * (1 + _kappa) * M * Q2 + _lambda * Matrix::Identity(m - (d+1), m - (d+1))).inverse() * Q2.transpose() * (_Y + _kappa * _X);
		B_ = P * 1.0 / (1 + _kappa) * R.inverse() * Q1.transpose() * ( _Y + _kappa * _X - (1 + _kappa) * M * A_); //Remember to permute B_ back with P!
		X_ = _X;

		std::cout << "qr solver formula subsolver TPS energy : " << _lambda * (A_.transpose() * (M + _lambda * Matrix::Identity(m, m)) * A_).trace() << std::endl;

		//Matrix equationL(m, m), solution(m, d), equationR(m, d);
		//equationL << Q1.transpose() * M * Q2, R,
		//	Q2.transpose() * M * Q2 + _lambda * Matrix::Identity(m-d-1, m-d-1), Matrix::Zero(m-d-1, d+1);
		//equationR << Q1.transpose() * _Y, Q2.transpose() * _Y;
		//
		//Eigen::JacobiSVD<Matrix> svd(equationL, Eigen::ComputeThinU | Eigen::ComputeThinV);
		//std::cout << "qr solver svd subsolver Condition Number: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
		//solution = svd.solve(equationR);

		////solution = equationL.inverse() * equationR;

		//A_ = Q2 * solution.block(0, 0, m-d-1, d);
		//B_ = P * solution.block(m-d-1, 0, d+1, d); //Remember to permute B_ back with P!
		//X_ = _X;

		//std::cout << "qr solver svd subsolver TPS energy : " << _lambda * A_.transpose() * (M + _lambda * Matrix::Identity(m, m)) * A_ << std::endl;
	}

	Matrix TPSFunction::evaluate(const Matrix &_X) const {
		Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);

		//Very space-consuming!
		//Matrix M = greenFunc(_X, X_);
		//return M * A_ + X * B_;

		Matrix rt = Matrix::Zero(_X.rows(), _X.cols());
		for (int i = 0; i < _X.rows(); i++) {
			for (int j = 0; j < X_.rows(); j++) {
				rt.row(i) += greenFunc( (_X.row(i) - X_.row(j)).norm(), _X.cols() ) * A_.row(j);
			}
		}
		rt += X * B_;

		return rt;
	}

#ifdef _DEBUG
	std::ostream& operator<<(std::ostream& os, const TPSFunction& tps){
		os << "X = \n" << tps.getX() << std::endl;
		os << "A = \n" << tps.getA() << std::endl;
		os << "B = \n" << tps.getB() << std::endl;
		return os;
	}
#endif
	
}