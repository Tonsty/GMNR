#include <iostream>

#ifdef _DEBUG
#include <QDebug>
#endif
#include <QList>
#include <QPoint>
#include <QMessageBox>

#include "../include/test_multitps_qt_dialog_control.h"

#include <GMNR/MultiTPS.h>

namespace gmnr{

	void TestMultiTPSDialogControl::on_dialog_updatePoints(QList<QPoint> &_points1, QList<QPoint> &_points2, QList<QPoint> &_points3){
		if (_points1.size() != 8){
			QMessageBox::warning(this, "Error", "wrong number of points1!!!");
			return;
		}

		if (_points2.size() != 8){
			QMessageBox::warning(this, "Error", "Error", "wrong number of points2!!!");
			return;
		}

		if (_points3.size() != 8){
			QMessageBox::warning(this, "Error", "Error", "wrong number of points3!!!");
			return;
		}

		gmnr::PointSet2D X, Y, Z	;
		QListPoint_to_PointSet2D(_points1, X);
		QListPoint_to_PointSet2D(_points2, Y);
		QListPoint_to_PointSet2D(_points3, Z);

		Vector2D t;
		t << -300, -300;
		for (int i = 0; i < 8; i++){
			X.row(i) += t;
			Y.row(i) += t;
			Z.row(i) += t;
		}
		X /= 300.0;
		Y /= 300.0;
		Z /= 300.0;

		gmnr::PointSet2D input_X(12, 2), input_Y(12, 2);
		input_X << X.block(4, 0, 4, 2), 
					Y.block(4, 0, 4, 2), 
					Z.block(4, 0, 4, 2); 
		input_Y << Y.block(0, 0, 4, 2),
					Z.block(0, 0, 4, 2),
					X.block(0, 0, 4, 2);

		std::vector<int> m(3);
		m[0] = 4;
		m[1] = 4;
		m[2] = 4;
		std::vector<int> alpha(3);
		alpha[0] = 0;
		alpha[1] = 1;
		alpha[2] = 2;
		std::vector<int> beta(3);
		beta[0] = 1;
		beta[1] = 2;
		beta[2] = 0;
		Vector kappa(3);
		kappa << 0.01, 0.01, 0.01;
		Vector lambda(3);
		lambda << 0.001, 0.001, 0.001;

		MultiTPS mtps(input_X, input_Y, m, alpha, beta, kappa, lambda);
		//std::cout << "mtps = \n" << mtps << std::endl;
		TPSFunction f1 = mtps.getfs()[0], f2 = mtps.getfs()[1], f3 = mtps.getfs()[2];

		X = f1.evaluate(X);
		Y = f2.evaluate(Y);
		Z = f3.evaluate(Z);

		Scalar XtA_norm_1 = (f1.getX().transpose() * f1.getA()).norm();
		std::cout << "XtA_norm_1 is: " << XtA_norm_1 << std::endl;

		Scalar XtA_norm_2 = (f2.getX().transpose() * f2.getA()).norm();
		std::cout << "XtA_norm_2 is: " << XtA_norm_2 << std::endl;

		Scalar XtA_norm_3 = (f1.getX().transpose() * f3.getA()).norm();
		std::cout << "XtA_norm_3 is: " << XtA_norm_1 << std::endl;

		X *= 300.0;
		Y *= 300.0;
		Z *= 300.0;
		for (int i = 0; i < 8; i++){
			X.row(i) -= t;
			Y.row(i) -= t;
			Z.row(i) -= t;
		}

		PointSet2D_to_QListPoint(X, _points1);
		PointSet2D_to_QListPoint(Y, _points2);
		PointSet2D_to_QListPoint(Z, _points3);

		repaint();
	}
};