#include <QPainter>
#include <QMouseEvent>
#include <QKeyEvent>

#include <unsupported/Eigen/Splines>

#include <GMNR/common.h>

#include "../include/test_multitps_qt_dialog.h"

namespace gmnr{


	typedef Eigen::Spline<Scalar, 2> Spline2D;
	typedef Eigen::SplineFitting<Spline2D> SplineFitting2D;

	TestMultiTPSDialog::TestMultiTPSDialog(QWidget *parent): QDialog(parent){
		setupUi(this);
		points1_active_ = false;
		points2_active_ = false;
		points3_active_ = false;

		setSample();
	}

	void TestMultiTPSDialog::setSample() {
		points1_.clear();
		for (int i = 0; i < 8; i++) {
			int x, y;
			x = 300 + 100 * cos(M_PI / 180.0 * i * 240.0 / 8);
			y = 300 + 100 * sin(M_PI / 180.0 * i * 240.0 / 8);
			QPoint point;
			point.setX(x);
			point.setY(y);
			points1_.push_back(point);
		}
		points2_.clear();
		for (int i = 0; i < 8; i++) {
			int x, y;
			x = 300 + 120 * cos(M_PI / 180.0 * (i+4) * 240.0 / 8);
			y = 300 + 120 * sin(M_PI / 180.0 * (i+4) * 240.0 / 8);
			QPoint point;
			point.setX(x);
			point.setY(y);
			points2_.push_back(point);
		}
		points3_.clear();
		for (int i = 0; i < 8; i++) {
			int x, y;
			x = 300 + 140 * cos(M_PI / 180.0 * (i+8) * 240.0 / 8);
			y = 300 + 140 * sin(M_PI / 180.0 * (i+8) * 240.0 / 8);
			QPoint point;
			point.setX(x);
			point.setY(y);
			points3_.push_back(point);
		}
	}

	void TestMultiTPSDialog::on_pushButton_reset_clicked() {
		setSample();
		repaint();
	}

	void TestMultiTPSDialog::QListPoint_to_PointSet2D(QList<QPoint> &_points, gmnr::PointSet2D &_pointset2d){
		_pointset2d.resize(_points.size(), Eigen::NoChange);
		for (int i = 0; i < _points.size(); i++){
			_pointset2d(i, 0) = _points[i].x();
			_pointset2d(i, 1) = _points[i].y();
		}
	}

	void TestMultiTPSDialog::PointSet2D_to_QListPoint(gmnr::PointSet2D &_pointset2d, QList<QPoint> &_points){
		_points.clear();
		for (int i = 0; i < _pointset2d.rows(); i++){
			QPoint point;
			point.setX(_pointset2d(i, 0));
			point.setY(_pointset2d(i, 1));
			_points.push_back(point);
		}
	}

	void TestMultiTPSDialog::paintEvent(QPaintEvent *event){
		QPainter painter(this);
		painter.setRenderHint(QPainter::Antialiasing, true);
		//painter.setRenderHint(QPainter::HighQualityAntialiasing, true);

		painter.setPen(Qt::NoPen);
		painter.setBrush(QBrush(Qt::red, Qt::SolidPattern));
		for (int i = 0; i < points1_.size(); i++){
			QPoint point = points1_[i];
			int radius = 8;
			painter.drawEllipse(point, radius, radius);
		}
		if(points1_.size() == 2){
			painter.setPen(Qt::red);
			painter.drawPolyline((QPoint*)(&points1_[0]), points1_.size());
		}

		painter.setPen(Qt::NoPen);
		painter.setBrush(QBrush(Qt::blue, Qt::SolidPattern));
		for (int i = 0; i < points2_.size(); i++){
			QPoint point = points2_[i];
			int radius = 6;
			painter.drawEllipse(point, radius, radius);
		}
		if(points2_.size() == 2) {
			painter.setPen(Qt::blue);
			painter.drawPolyline((QPoint*)(&points2_[0]), points2_.size());
		}

		painter.setPen(Qt::NoPen);
		painter.setBrush(QBrush(Qt::darkGreen, Qt::SolidPattern));
		for (int i = 0; i < points3_.size(); i++){
			QPoint point = points3_[i];
			int radius = 4;
			painter.drawEllipse(point, radius, radius);
		}
		if(points3_.size() == 2) {
			painter.setPen(Qt::darkGreen);
			painter.drawPolyline((QPoint*)(&points3_[0]), points3_.size());
		}

		int num_of_subdivision = 500;

		if(points1_.size() > 2){
			PointSet2D pointset2d_1;
			QListPoint_to_PointSet2D(points1_, pointset2d_1);
			Spline2D spline_2d_1 = SplineFitting2D::Interpolate(pointset2d_1.transpose(), std::min(points1_.size() -1, 3));
			PointSet2D pointset2d_1_more;
			pointset2d_1_more.resize(num_of_subdivision, Eigen::NoChange);
			for (int i = 0; i < num_of_subdivision; i++){
				pointset2d_1_more.row(i) = spline_2d_1(1.0/(num_of_subdivision-1) * i);
			}
			QList<QPoint> points1_more;
			PointSet2D_to_QListPoint(pointset2d_1_more, points1_more);
			painter.setPen(Qt::red);
			painter.drawPolyline((QPoint*)(&points1_more[0]), points1_more.size());
			//painter.drawPoints((QPoint*)(&points1_more[0]), points1_more.size());
		}

		if (points2_.size() > 2){
			PointSet2D pointset2d_2;
			QListPoint_to_PointSet2D(points2_, pointset2d_2);
			Spline2D cubic_spline_2d_2 = SplineFitting2D::Interpolate(pointset2d_2.transpose(), std::min(points2_.size() -1, 3));
			PointSet2D pointset2d_2_more;
			pointset2d_2_more.resize(num_of_subdivision, Eigen::NoChange);
			for (int i = 0; i < num_of_subdivision; i++){
				pointset2d_2_more.row(i) = cubic_spline_2d_2(1.0/(num_of_subdivision-1) * i);
			}
			QList<QPoint> points2_more;
			PointSet2D_to_QListPoint(pointset2d_2_more, points2_more);
			painter.setPen(Qt::blue);
			painter.drawPolyline((QPoint*)(&points2_more[0]), points2_more.size());
			//painter.drawPoints((QPoint*)(&points2_more[0]), points2_more.size());
		}

		if (points3_.size() > 2){
			PointSet2D pointset2d_3;
			QListPoint_to_PointSet2D(points3_, pointset2d_3);
			Spline2D cubic_spline_2d_3 = SplineFitting2D::Interpolate(pointset2d_3.transpose(), std::min(points3_.size() -1, 3));
			PointSet2D pointset2d_3_more;
			pointset2d_3_more.resize(num_of_subdivision, Eigen::NoChange);
			for (int i = 0; i < num_of_subdivision; i++){
				pointset2d_3_more.row(i) = cubic_spline_2d_3(1.0/(num_of_subdivision-1) * i);
			}
			QList<QPoint> points3_more;
			PointSet2D_to_QListPoint(pointset2d_3_more, points3_more);
			painter.setPen(Qt::darkGreen);
			painter.drawPolyline((QPoint*)(&points3_more[0]), points3_more.size());
			//painter.drawPoints((QPoint*)(&points3_more[0]), points3_more.size());
		}
	}

	void TestMultiTPSDialog::mousePressEvent(QMouseEvent *event){
		//#ifdef _DEBUG
		//		qDebug() << "x = " << event->pos().x() << " y = " << event->pos().y();
		//#endif
		switch(event->button()){
		case Qt::LeftButton:
			{
				if(points1_active_) points1_.push_back(event->pos());
				if(points2_active_) points2_.push_back(event->pos());
				if(points3_active_) points3_.push_back(event->pos());
			}
			break;
		case  Qt::RightButton:
			{
				if(points1_active_) points1_.removeLast();
				if(points2_active_) points2_.removeLast();
				if(points3_active_) points3_.removeLast();
			}
			break;
		}

		repaint();
		//QDialog::mousePressEvent(event);
	}

	void TestMultiTPSDialog::keyPressEvent(QKeyEvent *event){
		switch (event->key()){
		case Qt::Key_1:
			points1_active_ = true;
			break;
		case Qt::Key_2:
			points2_active_ = true;
			break;
		case Qt::Key_3:
			points3_active_ = true;
			break;
		default:
			QDialog::keyReleaseEvent(event);
		}
	}

	void TestMultiTPSDialog::keyReleaseEvent(QKeyEvent *event){
		switch (event->key()){
		case Qt::Key_1:
			points1_active_ = false;
			break;
		case Qt::Key_2:
			points2_active_ = false;
			break;
		case Qt::Key_3:
			points3_active_ = false;
			break;
		default:
			QDialog::keyReleaseEvent(event);
		}
	}

	void TestMultiTPSDialog::on_pushButton_multitps_clicked(){
		emit updatePoints(points1_, points2_, points3_);
	}

	void TestMultiTPSDialog::on_pushButton_clear_1_clicked(){
		points1_.clear();
		repaint();
	}

	void TestMultiTPSDialog::on_pushButton_clear_2_clicked(){
		points2_.clear();
		repaint();
	}

	void TestMultiTPSDialog::on_pushButton_clear_3_clicked(){
		points3_.clear();
		repaint();
	}

};

