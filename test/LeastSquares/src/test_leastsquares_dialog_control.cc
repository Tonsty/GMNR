#ifdef _DEBUG
#include <QDebug>
#endif
#include <QList>
#include <QPoint>
#include <QMessageBox>

#include "../include/test_leastsquares_dialog_control.h"

#include <GMNR/LeastSquares/PointToPoint.h>
#include <GMNR/LeastSquares/PointToPlane.h>

namespace gmnr{

	void TestLeastSquaresDialogControl::QListPoint_to_PointSet3D(QList<QPoint> &_points, gmnr::PointSet3D &_pointset3d){
		_pointset3d.resize(_points.size(), Eigen::NoChange);
		for (int i = 0; i < _points.size(); i++){
			_pointset3d(i, 0) = _points[i].x();
			_pointset3d(i, 1) = _points[i].y();
			_pointset3d(i, 2) = 1.0f;
		}
	}

	void TestLeastSquaresDialogControl::QListPoint_to_NormalSet3D(QList<QPoint> &_points, gmnr::NormalSet3D &_normalset3d){
		int n = _points.size();
		_normalset3d.resize(n, Eigen::NoChange);
		for (int i = 0; i < n; i++){
			_normalset3d(i, 0) = 0.0f;
			_normalset3d(i, 1) = -1.0;
			_normalset3d(i, 2) = 0.0f;
		}
	}

	void TestLeastSquaresDialogControl::PointSet3D_to_QListPoint(gmnr::PointSet3D &_pointset3d, QList<QPoint> &_points){
		_points.clear();
		for (int i = 0; i < _pointset3d.rows(); i++){
			QPoint point;
			point.setX(_pointset3d(i, 0));
			point.setY(_pointset3d(i, 1));
			_points.push_back(point);
		}
	}

	void TestLeastSquaresDialogControl::on_dialog_updatePoints(QList<QPoint> &_points1, QList<QPoint> &_points2, PushButton _pushButton){
		if (_points1.size() != _points2.size()){
			QMessageBox::warning(this, "Error", "numbers of source and target point must be the same!!!");
			return;
		}
		if (_points1.size() == 0){
			QMessageBox::warning(this, "Error", "no point added so far!!!");
			return;
		}

		switch(_pushButton){
		case PointToPointArun:{
				gmnr::PointSet3D X, Y;
				QListPoint_to_PointSet3D(_points1, X);
				QListPoint_to_PointSet3D(_points2, Y);
				gmnr::PointToPointArun arun(X, Y, true);
				X = (X * arun.rotation().transpose()).rowwise() + arun.translation().transpose();
				PointSet3D_to_QListPoint(X, _points1);
			}
			break;
		case PointToPointUmeyama:{
			if (_points1.size() < 3){
				QMessageBox::warning(this, "Error", "umeyama need 3 points at least!!!");
				return;
			}
			gmnr::PointSet3D X, Y;
			QListPoint_to_PointSet3D(_points1, X);
			QListPoint_to_PointSet3D(_points2, Y);
			gmnr::PointToPointUmeyama umeyama(X, Y, true, true);
			//gmnr::PointToPointUmeyama umeyama(X, Y, false, true);
			X = (X  * umeyama.scale() * umeyama.rotation().transpose()).rowwise() + umeyama.translation().transpose();
			PointSet3D_to_QListPoint(X, _points1);
			}
			break;
		case PointToPlaneLinear:{
			gmnr::PointSet3D X, Y;
			QListPoint_to_PointSet3D(_points1, X);
			QListPoint_to_PointSet3D(_points2, Y);
			gmnr::NormalSet3D N;
			QListPoint_to_NormalSet3D(_points2, N);
			gmnr::PointToPlaneLinear linear(X, Y, N, true);
			//X = (X * linear.rotationApproximate().transpose()).rowwise() + linear.translation().transpose();
			X = (X * linear.rotation().transpose()).rowwise() + linear.translation().transpose();
			PointSet3D_to_QListPoint(X, _points1);
			}
			break;
		}
		repaint();
	}
};