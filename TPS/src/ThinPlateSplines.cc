#include <iostream>

#include <GMNR/common.h>
#include <GMNR/math/GreenFunction.h>
#include <GMNR/time.h>
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
		//qr_solve_(_X, _Y, _lambda, _kappa);
		//formula_solve_(_X, _Y, _lambda, _kappa);
		direct_inverse_solve_(_X, _Y, _lambda, _kappa);
	}

	void TPSFunction::formula_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa){
		Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
		Matrix M = greenFunc(_X, _X);
		Matrix S = (1 + _kappa) * M + _lambda * Matrix::Identity(_X.rows(), _X.rows());

		Matrix inverseS = S.inverse();
		Matrix transposehomoX_inverseS = X.transpose() * inverseS;

		B_ = (transposehomoX_inverseS * (1+_kappa) * X).inverse() * transposehomoX_inverseS * (_Y + _kappa * _X);
		A_ = inverseS * (_Y + _kappa * _X - (1+_kappa) * X * B_);
		X_ = _X;

		int m = _X.rows(), d = _X.cols();

		std::cout << "formula solver TPS bending energy : " << _lambda * (A_.transpose() * M * A_).trace() << std::endl;	
	}

	void TPSFunction::direct_inverse_solve_(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa){
		int m = _X.rows(), d = _X.cols();
		
		Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
		Matrix M = greenFunc(_X, _X);
		Matrix S = (1 + _kappa) * M + _lambda * Matrix::Identity(m, m);

		Matrix equationL(m + (d+1), m + (d+1)), equationR(m + (d+1), d), solution(m + (d+1), d);
		equationL << S, (1 + _kappa) * X,
			X.transpose(), Matrix::Zero(d+1, d+1);
		equationR << _Y + _kappa * _X, Matrix::Zero(d+1, d);

		timestamp opt_start = now();

		//Eigen::FullPivLU<Matrix> fullLU(equationL);
		//solution = fullLU.solve(equationR);

		//Eigen::PartialPivLU<Matrix> partLU(equationL);
		//solution = partLU.solve(equationR);

		Eigen::HouseholderQR<Matrix> qr(equationL);
		solution = qr.solve(equationR);

		//Eigen::LDLT<Matrix> ldlt(equationL);
		//solution = ldlt.solve(equationR);

		float equation_solution_time = now() - opt_start;
		std::cout << "TPS equation " << equationL.rows() << " * " << equationL.cols() << " solution time: " << equation_solution_time << std::endl;

		A_ = solution.block(0, 0, m, d);
		B_ = solution.block(m, 0, d+1, d);
		X_ = _X;

		//std::cout << "direct inverse solver TPS bending energy : " << _lambda * (A_.transpose() * M * A_).trace() << std::endl;
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

		Matrix equationL = Q2.transpose() * (1 + _kappa) * M * Q2 + _lambda * Matrix::Identity(m - (d+1), m - (d+1));
		Eigen::FullPivLU<Matrix> fullLU(equationL);
		Matrix solution = fullLU.solve(Q2.transpose() * (_Y + _kappa * _X));
		A_ = Q2 * solution;
		B_ = P * 1.0 / (1 + _kappa) * R.inverse() * Q1.transpose() * ( _Y + _kappa * _X - (1 + _kappa) * M * A_); //Remember to permute B_ back with P!
		X_ = _X;

		std::cout << "qr solver TPS bending energy : " << _lambda * (A_.transpose() * M * A_).trace() << std::endl;
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

	ApproxiTPSFunction::ApproxiTPSFunction() : TPSFunction() {}

	ApproxiTPSFunction::ApproxiTPSFunction(const Matrix &_X, const Matrix &_Y, const Scalar &_lambda, const Scalar &_kappa, const int &n) {
		if (n < 0 || n > _X.rows()) {
			TPSFunction tps = TPSFunction(_X, _Y, _lambda, _kappa);
			X_ = tps.getX();
			A_ = tps.getA();
			B_ = tps.getB();
		} else {
			int m = _X.rows(), d = _X.cols();

			Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X);
			Matrix M = greenFunc(_X, _X.topRows(n));
			Matrix S = (1 + _kappa) * M + _lambda * Matrix::Identity(m, n);

			Matrix equationL(m + (d+1), n + (d+1)), equationR(m + (d+1), d), solution(n + (d+1), d);
			equationL << S, (1 + _kappa) * X,
				X.topRows(n).transpose(), Matrix::Zero(d+1, d+1);
			equationR << _Y + _kappa * _X, Matrix::Zero(d+1, d);

			timestamp opt_start = now();

			//Eigen::FullPivLU<Matrix> fullLU(equationL);
			//solution = fullLU.solve(equationR);

			Eigen::HouseholderQR<Matrix> qr(equationL);
			solution = qr.solve(equationR);
			
			float equation_solution_time = now() - opt_start;
			std::cout << "TPS equation " << equationL.rows() << " * " << equationL.cols() << " solution time: " << equation_solution_time << std::endl;

			A_ = solution.block(0, 0, n, d);
			B_ = solution.block(n, 0, d+1, d);
			X_ = _X.topRows(n);

			//std::cout << "direct inverse solver TPS bending energy : " << _lambda * (A_.transpose() * M * A_).trace() << std::endl;
		}
	}

	std::ostream& operator<<(std::ostream& os, const TPSFunction& tps){
		os << "X = \n" << tps.getX() << std::endl;
		os << "A = \n" << tps.getA() << std::endl;
		os << "B = \n" << tps.getB() << std::endl;
		return os;
	}
	
}