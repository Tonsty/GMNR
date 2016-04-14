#include <QPainter>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QMessageBox>

#include <unsupported/Eigen/Splines>

#include <GMNR/common.h>

#include "../include/test_multitps_qt_dialog.h"

namespace gmnr{


	typedef Eigen::Spline<Scalar, 2> Spline2D;
	typedef Eigen::SplineFitting<Spline2D> SplineFitting2D;

	TestMultiTPSDialog::TestMultiTPSDialog(QWidget *parent): QDialog(parent){
		setupUi(this);

		setDefaultSample();

		points1_last_ = points1_;
		points2_last_ = points2_;
		points3_last_ = points3_;

		after_multitps = false;
	}

	void TestMultiTPSDialog::setDefaultSample() {
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

	void TestMultiTPSDialog::setLastSample() {
		points1_ = points1_last_;
		points2_ = points2_last_;
		points3_ = points3_last_;
	}

	void TestMultiTPSDialog::on_pushButton_reset_d_clicked() {
		setDefaultSample();
		points1_last_ = points1_;
		points2_last_ = points2_;
		points3_last_ = points3_;
		after_multitps = false;
		repaint();
	}

	void TestMultiTPSDialog::on_pushButton_reset_l_clicked() {
		setLastSample();
		after_multitps = false;
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

	void TestMultiTPSDialog::paint(QList<QPoint> &_points1, QList<QPoint> &_points2, QList<QPoint> &_points3, bool _draw_correspondences, float opacity) {
		QPainter painter(this);
		painter.setRenderHint(QPainter::Antialiasing, true);
		//painter.setRenderHint(QPainter::HighQualityAntialiasing, true);

		painter.setOpacity(opacity);

		//Qt::GlobalColor colors[] = {Qt::red, Qt::blue, Qt::magenta};
		QColor colors[3];
		colors[0] = QColor(255, 0, 0);
		colors[1] = QColor(0, 0, 255);
		colors[2] = QColor(255, 128, 0);

		painter.setPen(Qt::NoPen);
		painter.setBrush(QBrush(colors[0], Qt::SolidPattern));
		for (int i = 0; i < _points1.size(); i++){
			QPoint point = _points1[i];
			int radius = 8;
			painter.drawEllipse(point, radius, radius);
		}
		if(_points1.size() == 2){
			painter.setPen(colors[0]);
			painter.drawPolyline((QPoint*)(&_points1[0]), _points1.size());
		}

		painter.setPen(Qt::NoPen);
		painter.setBrush(QBrush(colors[1], Qt::SolidPattern));
		for (int i = 0; i < _points2.size(); i++){
			QPoint point = _points2[i];
			int radius = 6;
			painter.drawEllipse(point, radius, radius);
		}
		if(_points2.size() == 2) {
			painter.setPen(colors[1]);
			painter.drawPolyline((QPoint*)(&_points2[0]), _points2.size());
		}

		painter.setPen(Qt::NoPen);
		painter.setBrush(QBrush(colors[2], Qt::SolidPattern));
		for (int i = 0; i < _points3.size(); i++){
			QPoint point = _points3[i];
			int radius = 4;
			painter.drawEllipse(point, radius, radius);
		}
		if(points3_.size() == 2) {
			painter.setPen(colors[2]);
			painter.drawPolyline((QPoint*)(&_points3[0]), _points3.size());
		}

		if (_draw_correspondences) {
			if (_points1.size() == 8 && _points2.size() == 8 && _points3.size() == 8) {
				painter.setPen(QPen(QBrush(QColor(0, 176, 80)), 2, Qt::DotLine));

				for (int i = _points1.size()/2, j = 0; i < _points1.size(); i++, j++) {
					painter.drawLine(_points1[i], _points2[j]);
				}

				for (int i = _points2.size()/2, j = 0; i < _points2.size(); i++, j++) {
					painter.drawLine(_points2[i], _points3[j]);
				}

				for (int i = _points3.size()/2, j = 0; i < _points3.size(); i++, j++) {
					painter.drawLine(_points3[i], _points1[j]);
				}
			}
		}

		int num_of_subdivision = 500;

		if(_points1.size() > 2){
			PointSet2D pointset2d_1;
			QListPoint_to_PointSet2D(_points1, pointset2d_1);
			Spline2D spline_2d_1 = SplineFitting2D::Interpolate(pointset2d_1.transpose(), std::min(_points1.size() -1, 3));
			PointSet2D pointset2d_1_more;
			pointset2d_1_more.resize(num_of_subdivision, Eigen::NoChange);
			for (int i = 0; i < num_of_subdivision; i++){
				pointset2d_1_more.row(i) = spline_2d_1(1.0/(num_of_subdivision-1) * i);
			}
			QList<QPoint> points1_more;
			PointSet2D_to_QListPoint(pointset2d_1_more, points1_more);
			painter.setPen(QPen(QBrush(colors[0]), 3));
			painter.drawPolyline((QPoint*)(&points1_more[0]), points1_more.size());
			//painter.drawPoints((QPoint*)(&points1_more[0]), points1_more.size());
		}

		if (_points2.size() > 2){
			PointSet2D pointset2d_2;
			QListPoint_to_PointSet2D(_points2, pointset2d_2);
			Spline2D cubic_spline_2d_2 = SplineFitting2D::Interpolate(pointset2d_2.transpose(), std::min(_points2.size() -1, 3));
			PointSet2D pointset2d_2_more;
			pointset2d_2_more.resize(num_of_subdivision, Eigen::NoChange);
			for (int i = 0; i < num_of_subdivision; i++){
				pointset2d_2_more.row(i) = cubic_spline_2d_2(1.0/(num_of_subdivision-1) * i);
			}
			QList<QPoint> points2_more;
			PointSet2D_to_QListPoint(pointset2d_2_more, points2_more);
			painter.setPen(QPen(QBrush(colors[1]), 3));
			painter.drawPolyline((QPoint*)(&points2_more[0]), points2_more.size());
			//painter.drawPoints((QPoint*)(&points2_more[0]), points2_more.size());
		}

		if (_points3.size() > 2){
			PointSet2D pointset2d_3;
			QListPoint_to_PointSet2D(_points3, pointset2d_3);
			Spline2D cubic_spline_2d_3 = SplineFitting2D::Interpolate(pointset2d_3.transpose(), std::min(_points3.size() -1, 3));
			PointSet2D pointset2d_3_more;
			pointset2d_3_more.resize(num_of_subdivision, Eigen::NoChange);
			for (int i = 0; i < num_of_subdivision; i++){
				pointset2d_3_more.row(i) = cubic_spline_2d_3(1.0/(num_of_subdivision-1) * i);
			}
			QList<QPoint> points3_more;
			PointSet2D_to_QListPoint(pointset2d_3_more, points3_more);
			painter.setPen(QPen(QBrush(colors[2]), 3));
			painter.drawPolyline((QPoint*)(&points3_more[0]), points3_more.size());
			//painter.drawPoints((QPoint*)(&points3_more[0]), points3_more.size());
		}
	}

	void TestMultiTPSDialog::on_checkBox_showLast_stateChanged(int state) {
		repaint();
	}

	void TestMultiTPSDialog::paintEvent(QPaintEvent *event){
		if(checkBox_showLast->checkState() == Qt::Checked && after_multitps) paint(points1_last_, points2_last_, points3_last_, true, 0.2f);
		paint(points1_, points2_, points3_, !after_multitps);
	}

	void TestMultiTPSDialog::mousePressEvent(QMouseEvent *event){
		//#ifdef _DEBUG
		//		qDebug() << "x = " << event->pos().x() << " y = " << event->pos().y();
		//#endif
		
		if(QApplication::keyboardModifiers () != Qt::ControlModifier) return;
		
		switch(event->button()){
		case Qt::LeftButton:
			{
				if(radioButton_1->isChecked()) {
					if(points1_.size() < 8) points1_.push_back(event->pos());
					else {
						QMessageBox::warning(this, "Error", "the number of points of view1 cannot be more than 8!!!");
					}
				}
				if(radioButton_2->isChecked()) {
					if(points2_.size() < 8) points2_.push_back(event->pos());
					else {
						QMessageBox::warning(this, "Error", "the number of points of view2 cannot be more than 8!!!");
					}
				}
				if(radioButton_3->isChecked()) {
					if(points3_.size() < 8) points3_.push_back(event->pos());
					else {
						QMessageBox::warning(this, "Error", "the number of points of view3 cannot be more than 8!!!");
					}
				}
				//if(points1_active_ || points2_active_ || points3_active_) after_multitps = false;
			}
			break;
		case  Qt::RightButton:
			{
				if(radioButton_1->isChecked() && !points1_.empty() ) points1_.removeLast();
				if(radioButton_2->isChecked() && !points2_.empty() ) points2_.removeLast();
				if(radioButton_3->isChecked() && !points3_.empty() ) points3_.removeLast();
				after_multitps = false;
			}
			break;
		}

		repaint();
		//QDialog::mousePressEvent(event);
	}

	void TestMultiTPSDialog::on_pushButton_multitps_clicked(){
		points1_last_ = points1_;
		points2_last_ = points2_;
		points3_last_ = points3_;
		after_multitps = true;
		emit updatePoints(points1_, points2_, points3_);
	}

	void TestMultiTPSDialog::on_pushButton_clear_1_clicked(){
		points1_.clear();
		after_multitps = false;
		repaint();
	}

	void TestMultiTPSDialog::on_pushButton_clear_2_clicked(){
		points2_.clear();
		after_multitps = false;
		repaint();
	}

	void TestMultiTPSDialog::on_pushButton_clear_3_clicked(){
		points3_.clear();
		after_multitps = false;
		repaint();
	}

};

