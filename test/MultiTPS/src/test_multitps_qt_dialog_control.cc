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
			QMessageBox::warning(this, "Error", "the number of points of view1 must be equal to 8!!!");
			return;
		}

		if (_points2.size() != 8){
			QMessageBox::warning(this, "Error", "the number of points of view2 must be equal to 8!!!");
			return;
		}

		if (_points3.size() != 8){
			QMessageBox::warning(this, "Error", "the number of points of view3 must be equal to 8!!!");
			return;
		}

		gmnr::PointSet2D X, Y, Z	;
		QListPoint_to_PointSet2D(_points1, X);
		QListPoint_to_PointSet2D(_points2, Y);
		QListPoint_to_PointSet2D(_points3, Z);

		float kappa_1 = kappa_1_dsb->value();
		float kappa_2 = kappa_2_dsb->value();
		float kappa_3 = kappa_3_dsb->value();
		float lambda_1 = lambda_1_dsb->value();
		float lambda_2 = lambda_2_dsb->value();
		float lambda_3 = lambda_3_dsb->value();
		bool iterative = checkBox_iterative->isChecked();
		int max_iter_num = iterative ? max_iter_num_sb->value() : 0;
		float iter_rate = iter_rate_dsb->value();

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

		//gmnr::PointSet2D input_X(12, 2), input_Y(12, 2);
		//input_X << X.block(4, 0, 4, 2), 
		//	Y.block(4, 0, 4, 2), 
		//	Z.block(4, 0, 4, 2); 
		//input_Y << Y.block(0, 0, 4, 2),
		//	Z.block(0, 0, 4, 2),
		//	X.block(0, 0, 4, 2);

		//std::vector<int> m(3);
		//m[0] = 4;
		//m[1] = 4;
		//m[2] = 4;
		//std::vector<int> alpha(3);
		//alpha[0] = 0;
		//alpha[1] = 1;
		//alpha[2] = 2;
		//std::vector<int> beta(3);
		//beta[0] = 1;
		//beta[1] = 2;
		//beta[2] = 0;
		//Vector kappa(3);
		//kappa << kappa_1, kappa_2, kappa_3;
		//Vector lambda(3);
		//lambda << lambda_1, lambda_2, lambda_3;

		gmnr::PointSet2D input_X(24, 2), input_Y(24, 2);
		input_X << X.block(4, 0, 4, 2), 
			Y.block(4, 0, 4, 2), 
			Z.block(4, 0, 4, 2),
			X.block(4, 0, 4, 2), 
			Y.block(4, 0, 4, 2), 
			Z.block(4, 0, 4, 2); 
		input_Y << Y.block(0, 0, 4, 2),
			Z.block(0, 0, 4, 2),
			X.block(0, 0, 4, 2),
			Y.block(0, 0, 4, 2),
			Z.block(0, 0, 4, 2),
			X.block(0, 0, 4, 2);

		std::vector<int> m(6);
		m[0] = 4;
		m[1] = 4;
		m[2] = 4;
		m[3] = 4;
		m[4] = 4;
		m[5] = 4;
		std::vector<int> alpha(6);
		alpha[0] = 0;
		alpha[1] = 1;
		alpha[2] = 2;
		alpha[3] = 0;
		alpha[4] = 1;
		alpha[5] = 2;
		std::vector<int> beta(6);
		beta[0] = 1;
		beta[1] = 2;
		beta[2] = 0;
		beta[3] = 1;
		beta[4] = 2;
		beta[5] = 0;
		Vector kappa(3);
		kappa << kappa_1, kappa_2, kappa_3;
		Vector lambda(3);
		lambda << lambda_1, lambda_2, lambda_3;

		MultiTPS mtps(input_X, input_Y, m, alpha, beta, kappa, lambda, max_iter_num, iter_rate);
		//std::cout << "mtps = \n" << mtps << std::endl;
		TPSFunction f1 = mtps.getfs()[0], f2 = mtps.getfs()[1], f3 = mtps.getfs()[2];

		X = f1.evaluate(X);
		Y = f2.evaluate(Y);
		Z = f3.evaluate(Z);

		Scalar XtA_norm_1 = (f1.getX().transpose() * f1.getA()).norm();
		std::cout << "XtA_norm_1 is: " << XtA_norm_1 << std::endl;

		Scalar XtA_norm_2 = (f2.getX().transpose() * f2.getA()).norm();
		std::cout << "XtA_norm_2 is: " << XtA_norm_2 << std::endl;

		Scalar XtA_norm_3 = (f3.getX().transpose() * f3.getA()).norm();
		std::cout << "XtA_norm_3 is: " << XtA_norm_3 << std::endl;

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