#include <iostream>

#include <GMNR/DualTPS.h>

namespace gmnr{

	DualTPS::DualTPS(){}

	DualTPS::DualTPS(
		const Matrix &_X, 
		const Matrix &_Y, 
		const Scalar &_kappa1, 
		const Scalar &_kappa2, 
		const Scalar &_lambda1, 
		const Scalar &_lambda2){
			qr_solve_(_X, _Y, _kappa1, _kappa2, _lambda1, _lambda2);			
			//formula_solve_(_X, _Y, _kappa1, _kappa2, _lambda1, _lambda2);
			//svd_solve_(_X, _Y, _kappa1, _kappa2, _lambda1, _lambda2);
			//direct_inverse_solve_(_X, _Y, _kappa1, _kappa2, _lambda1, _lambda2);
	}

	void DualTPS::formula_solve_(const Matrix &_X, 
		const Matrix &_Y, 
		const Scalar &_kappa1, 
		const Scalar &_kappa2, 
		const Scalar &_lambda1, 
		const Scalar &_lambda2){

			Matrix M1, M2;
			Matrix S1, S2, T1, T2, Z1, Z2;
			Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X), Y = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_Y);

			M1 = greenFunc(_X, _X);
			M2 = greenFunc(_Y, _Y);
			S1 = (1 + _kappa1) * M1 + _lambda1 * Matrix::Identity(_X.rows(), _X.rows());
			S2 = (1 + _kappa2) * M2 + _lambda2 * Matrix::Identity(_Y.rows(), _Y.rows());
			T1 = (1 + _kappa1) * X;
			T2 = (1 + _kappa2) * Y;
			Z1 = _kappa1 * _X;
			Z2 = _kappa2 * _Y;

			Matrix P1, P2, P3, P4, P5, P6, P7, P8;
			Eigen::JacobiSVD<Matrix> svd;
			svd.compute(S1);
			std::cout << "formula Condition Number 1: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			Matrix S1in = S1.inverse();
			Matrix M1S1in = M1 * S1in;
			P1 = S2 - M1S1in * M2;
			P2 = T2 - M1S1in * Y;
			P3 = X - M1S1in * T1;
			P4 = M1S1in * Z1 + Z2;
			Matrix Xtr = X.transpose();
			Matrix XtrS1in = Xtr * S1in;
			svd.compute(XtrS1in * T1);
			std::cout << "formula Condition Number 2: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			P5 = (XtrS1in * T1).inverse() * XtrS1in;
			Matrix P3P5 = P3 * P5;
			P6 = P1 - P3P5 * M2;
			P7 = P3P5 * Y - P2;
			P8 = P3P5 * Z1 + P4;
			Matrix Ytr = Y.transpose();
			svd.compute(P6);
			std::cout << "formula Condition Number 3: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			Matrix P6in = P6.inverse();
			Matrix YtrP6in = Ytr * P6in;

			Matrix W1, W2, B1, B2, P9;
			svd.compute(YtrP6in * P7);
			std::cout << "formula Condition Number 4: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			B2 = (YtrP6in * P7).inverse() * YtrP6in * P8 * (-1);
			W2 = P6in * (P7 * B2 + P8);
			P9 = Y * B2 + M2 * W2 + Z1;
			B1 = P5 * P9;
			W1 = S1in * (P9 - T1 * B1);

			f1_ = TPSFunction(_X, W1, B1);
			f2_ = TPSFunction(_Y, W2, B2);

			std::cout << "formula solver RTPS energy : " << 
			(X * B1 + M1 * W1 - Y * B2 - M2 * W2).squaredNorm() 
				+ _kappa1 * (X * B1 + M1 * W1 - _X).squaredNorm() 
				+ _kappa2 * (Y * B1 + M2 * W2 - _Y).squaredNorm()
				+ _lambda1 * (W1.transpose() * M1 * W1).trace()
				+ _lambda2 * (W2.transpose() * M2 * W2).trace() << std::endl;
	}

	void DualTPS::direct_inverse_solve_(const Matrix &_X, 
		const Matrix &_Y, 
		const Scalar &_kappa1, 
		const Scalar &_kappa2, 
		const Scalar &_lambda1, 
		const Scalar &_lambda2){

			Matrix M1, M2;
			Matrix S1, S2, T1, T2, Z1, Z2;
			Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X), Y = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_Y);

			M1 = greenFunc(_X, _X);
			M2 = greenFunc(_Y, _Y);
			S1 = (1 + _kappa1) * M1 + _lambda1 * Matrix::Identity(_X.rows(), _X.rows());
			S2 = (1 + _kappa2) * M2 + _lambda2 * Matrix::Identity(_Y.rows(), _Y.rows());
			T1 = (1 + _kappa1) * X;
			T2 = (1 + _kappa2) * Y;
			Z1 = _kappa1 * _X;
			Z2 = _kappa2 * _Y;

			int m = _X.rows(), d = _X.cols();

			Matrix equationL(2*m + 2*(d+1), 2*m + 2*(d+1)), equationR(2*m + 2*(d+1), d), solution(2*m + 2*(d+1), d);
			equationL << S1, -M2, T1, -Y,
				-M1, S2, -X, T2,
				X.transpose(), Matrix::Zero(d+1, m), Matrix::Zero(d+1, d+1), Matrix::Zero(d+1, d+1),
				Matrix::Zero(d+1, m), Y.transpose(), Matrix::Zero(d+1, d+1), Matrix::Zero(d+1, d+1);
			equationR << Z1, Z2, Matrix::Zero(d+1, d), Matrix::Zero(d+1, d);

			Eigen::JacobiSVD<Matrix> svd(equationL);
			std::cout << "direct inverse Condition Number: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;

			solution = equationL.inverse() * equationR;

			Matrix W1, W2, B1, B2;
			W1 = solution.block(0, 0, m, d);
			W2= solution.block(m, 0, m, d);
			B1 = solution.block(2*m, 0, d+1, d);
			B2= solution.block(2*m + (d+1), 0, d+1, d);

			f1_ = TPSFunction(_X, W1, B1);
			f2_ = TPSFunction(_Y, W2, B2);

			std::cout << "direct inverse solver RTPS energy : " << 
				(X * B1 + M1 * W1 - Y * B2 - M2 * W2).squaredNorm() 
				+ _kappa1 * (X * B1 + M1 * W1 - _X).squaredNorm() 
				+ _kappa2 * (Y * B1 + M2 * W2 - _Y).squaredNorm()
				+ _lambda1 * (W1.transpose() * M1 * W1).trace()
				+ _lambda2 * (W2.transpose() * M2 * W2).trace() << std::endl;
	}

	void DualTPS::svd_solve_(const Matrix &_X, 
		const Matrix &_Y, 
		const Scalar &_kappa1, 
		const Scalar &_kappa2, 
		const Scalar &_lambda1, 
		const Scalar &_lambda2){

			Matrix M1, M2;
			Matrix S1, S2, T1, T2, Z1, Z2;
			Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X), Y = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_Y);

			M1 = greenFunc(_X, _X);
			M2 = greenFunc(_Y, _Y);
			S1 = (1 + _kappa1) * M1 + _lambda1 * Matrix::Identity(_X.rows(), _X.rows());
			S2 = (1 + _kappa2) * M2 + _lambda2 * Matrix::Identity(_Y.rows(), _Y.rows());
			T1 = (1 + _kappa1) * X;
			T2 = (1 + _kappa2) * Y;
			Z1 = _kappa1 * _X;
			Z2 = _kappa2 * _Y;

			int m = _X.rows(), d = _X.cols();

			Matrix equationL(2*m + 2*(d+1), 2*m + 2*(d+1)), equationR(2*m + 2*(d+1), d), solution(2*m + 2*(d+1), d);
			equationL << S1, -M2, T1, -Y,
				-M1, S2, -X, T2,
				X.transpose(), Matrix::Zero(d+1, m), Matrix::Zero(d+1, d+1), Matrix::Zero(d+1, d+1),
				Matrix::Zero(d+1, m), Y.transpose(), Matrix::Zero(d+1, d+1), Matrix::Zero(d+1, d+1);
			equationR << Z1, Z2, Matrix::Zero(d+1, d), Matrix::Zero(d+1, d);

			solution = equationL.inverse() * equationR;

			Matrix W1, W2, B1, B2;
			W1 = solution.block(0, 0, m, d);
			W2= solution.block(m, 0, m, d);
			B1 = solution.block(2*m, 0, d+1, d);
			B2= solution.block(2*m + (d+1), 0, d+1, d);
			f1_ = TPSFunction(_X, W1, B1);
			f2_ = TPSFunction(_Y, W2, B2);

			Eigen::JacobiSVD<Matrix> svd(equationL, Eigen::ComputeThinU | Eigen::ComputeThinV);
			std::cout << "svd Condition Number: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			//std::cout << "Rows = " << equationL.rows() << " Columns = " << equationL.cols() << std::endl;
			//std::cout << "before reset threshold, svd_rank = " << svd.rank() << std::endl;
			//svd.setThreshold(Eigen::NumTraits<Scalar>::epsilon());
			//std::cout << "after reset threshold, svd_rank = " << svd.rank() << std::endl;
			solution = svd.solve(equationR);

			W1 = solution.block(0, 0, m, d);
			W2= solution.block(m, 0, m, d);
			B1 = solution.block(2*m, 0, d+1, d);
			B2= solution.block(2*m + (d+1), 0, d+1, d);

			f1_ = TPSFunction(_X, W1, B1);
			f2_ = TPSFunction(_Y, W2, B2);

			std::cout << "svd solver RTPS energy : " << 
				(X * B1 + M1 * W1 - Y * B2 - M2 * W2).squaredNorm() 
				+ _kappa1 * (X * B1 + M1 * W1 - _X).squaredNorm() 
				+ _kappa2 * (Y * B1 + M2 * W2 - _Y).squaredNorm()
				+ _lambda1 * (W1.transpose() * M1 * W1).trace()
				+ _lambda2 * (W2.transpose() * M2 * W2).trace() << std::endl;
	}

	void DualTPS::qr_solve_(const Matrix &_X, 
		const Matrix &_Y, 
		const Scalar &_kappa1, 
		const Scalar &_kappa2, 
		const Scalar &_lambda1, 
		const Scalar &_lambda2){
			Matrix M1, M2;
			Matrix X = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_X), Y = Eigen::Homogeneous<Matrix, Eigen::Horizontal>(_Y);

			M1 = greenFunc(_X, _X);
			M2 = greenFunc(_Y, _Y);

			int m = _X.rows(), d = _X.cols();

			Eigen::ColPivHouseholderQR<Matrix> qr_x(X), qr_y(Y);
			Matrix Px = qr_x.colsPermutation();
			Matrix Py = qr_y.colsPermutation();
			Matrix Qx = qr_x.matrixQ(), Rx = qr_x.matrixR().topLeftCorner(d+1, d+1).template triangularView<Eigen::Upper>();
			Matrix Qy = qr_y.matrixQ(), Ry = qr_y.matrixR().topLeftCorner(d+1, d+1).template triangularView<Eigen::Upper>();
			Matrix Qx1 = Qx.block(0, 0, m, d+1), Qx2 = Qx.block(0, d+1, m, m-d-1);
			Matrix Qy1 = Qy.block(0, 0, m, d+1), Qy2 = Qy.block(0, d+1, m, m-d-1);
			Matrix Rx_hat = (Rx * Px.transpose()).block(0, 0, d+1, d), Ry_hat = (Ry * Py.transpose()).block(0, 0, d+1, d);

			Matrix equationL(2*m, 2*m), solution(2*m, d), equationR(2*m, d);
			equationL <<
				(1 + _kappa1) * Qx1.transpose() * M1 * Qx2, 
				-Qx1.transpose() * M2 * Qy2, 
				(1+_kappa1) * Rx, 
				-Qx1.transpose() * Y * Py, //Permute Y

				-Qy1.transpose() * M1 * Qx2, 
				(1 + _kappa2) * Qy1.transpose() * M2 * Qy2, 
				-Qy1.transpose() * X * Px,
				(1+_kappa2) * Ry,//Permute X

				(1 + _kappa1) * Qx2.transpose() * M1 * Qx2 + _lambda1 * Matrix::Identity(m-d-1, m-d-1), 
				-Qx2.transpose() * M2 * Qy2, 
				Matrix::Zero(m-d-1, d+1),  
				-Qx2.transpose() * Y * Py, //Permute Y

				-Qy2.transpose() * M1 * Qx2, 
				(1 + _kappa2) * Qy2.transpose() * M2 * Qy2 + _lambda2 * Matrix::Identity(m-d-1, m-d-1), 
				-Qy2.transpose() * X * Px, 
				Matrix::Zero(m-d-1, d+1); //Permute X

			equationR << _kappa1 * Rx_hat, _kappa2 * Ry_hat, Matrix::Zero(m-d-1, d), Matrix::Zero(m-d-1, d);

			Eigen::JacobiSVD<Matrix> svd(equationL, Eigen::ComputeThinU | Eigen::ComputeThinV);
			std::cout << "qr Condition Number: " << svd.singularValues()[0] / svd.singularValues().reverse()[0] << std::endl;
			solution = svd.solve(equationR);

			//solution = equationL.inverse() * equationR;

			Matrix W1, W2, Gammax, Gammay, B1, B2;
			Gammax = solution.block(0, 0, m-d-1, d);
			Gammay = solution.block(m-d-1, 0, m-d-1, d);
			W1 = Qx2 * Gammax;
			W2 = Qy2 * Gammay;
			B1 = Px * solution.block(2*(m-d-1), 0, d+1, d);
			B2 = Py * solution.block(2*m-d-1, 0, d+1, d);

			f1_ = TPSFunction(_X, W1, B1);
			f2_ = TPSFunction(_Y, W2, B2);

			std::cout << "qr solver RTPS energy : " << 
				(X * B1 + M1 * W1 - Y * B2 - M2 * W2).squaredNorm() 
				+ _kappa1 * (X * B1 + M1 * W1 - _X).squaredNorm() 
				+ _kappa2 * (Y * B1 + M2 * W2 - _Y).squaredNorm()
				+ _lambda1 * (W1.transpose() * M1 * W1).trace()
				+ _lambda2 * (W2.transpose() * M2 * W2).trace() << std::endl;
	}

#ifdef _DEBUG
	std::ostream& operator<<(std::ostream& os, const DualTPS& tps){
		os << "f1 = \n" << tps.getf1() << std::endl;
		os << "f2 = \n" << tps.getf2() << std::endl;

		return os;
	}
#endif

}