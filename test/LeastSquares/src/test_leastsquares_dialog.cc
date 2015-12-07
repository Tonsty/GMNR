#include <QPainter>
#include <QMouseEvent>
#include <QKeyEvent>

#include "../include/test_leastsquares_dialog.h"

namespace gmnr{

	void TestLeastSquaresDialog::paintEvent(QPaintEvent *event){
		QPainter painter(this);
		painter.setRenderHint(QPainter::Antialiasing, true);
		painter.setPen(Qt::NoPen);

		painter.setBrush(QBrush(Qt::red, Qt::SolidPattern));
		for (int i = 0; i < points1_.size(); i++){
			QPoint point = points1_[i];
			int radius = 8;
			painter.drawEllipse(point, radius, radius);
		}
		if(points1_.size() > 1){
			painter.setPen(Qt::red);
			painter.drawPolyline((QPoint*)(&points1_[0]), points1_.size());
		}

		painter.setPen(Qt::NoPen);
		painter.setBrush(QBrush(Qt::blue, Qt::SolidPattern));
		for (int i = 0; i < points2_.size(); i++){
			QPoint point = points2_[i];
			int radius = 8;
			painter.drawEllipse(point, radius, radius);
		}
		if(points2_.size() > 1) {
			painter.setPen(Qt::blue);
			painter.drawPolyline((QPoint*)(&points2_[0]), points2_.size());
		}
	}

	void TestLeastSquaresDialog::mousePressEvent(QMouseEvent *event){
		//#ifdef _DEBUG
		//		qDebug() << "x = " << event->pos().x() << " y = " << event->pos().y();
		//#endif
		switch(event->button()){
		case Qt::LeftButton:
			{
				if(points1_active_) points1_.push_back(event->pos());
				if(points2_active_) points2_.push_back(event->pos());
			}
			break;
		case  Qt::RightButton:
			{
				if(points1_active_) points1_.removeLast();
				if(points2_active_) points2_.removeLast();
			}
			break;
		}

		repaint();
		//QDialog::mousePressEvent(event);
	}

	void TestLeastSquaresDialog::keyPressEvent(QKeyEvent *event){
		switch (event->key()){
		case Qt::Key_1:
			points1_active_ = true;
			break;
		case Qt::Key_2:
			points2_active_ = true;
		default:
			QDialog::keyReleaseEvent(event);
		}
	}

	void TestLeastSquaresDialog::keyReleaseEvent(QKeyEvent *event){
		switch (event->key()){
		case Qt::Key_1:
			points1_active_ = false;
			break;
		case Qt::Key_2:
			points2_active_ = false;
		default:
			QDialog::keyReleaseEvent(event);
		}
	}

	void TestLeastSquaresDialog::on_pushButton_PointToPointArun_clicked(){
		emit updatePoints(points1_, points2_, PointToPointArun);
	}

	void TestLeastSquaresDialog::on_pushButton_PointToPointUmeyama_clicked(){
		emit updatePoints(points1_, points2_, PointToPointUmeyama);
	}

	void TestLeastSquaresDialog::on_pushButton_PointToPlaneLinear_clicked(){
		emit updatePoints(points1_, points2_, PointToPlaneLinear);
	}

	void TestLeastSquaresDialog::on_pushButton_clear_source_clicked(){
		points1_.clear();
		repaint();
	}

	void TestLeastSquaresDialog::on_pushButton_clear_target_clicked(){
		points2_.clear();
		repaint();
	}

};

