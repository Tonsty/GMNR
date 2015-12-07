#ifndef GMNR_TEST_TEST_LEASTSQUARES_DIALOG_H
#define GMNR_TEST_TEST_LEASTSQUARES_DIALOG_H

#include <QDialog>
#include <QPoint>

#include "ui_test_leastsquares_dialog.h"

class QWidget;
class QPaintEvent;
class QKeyEvent;

namespace gmnr{

	class TestLeastSquaresDialog: public QDialog, public Ui::test_leastsquares{	

		Q_OBJECT

	public:
		TestLeastSquaresDialog(QWidget *parent = 0): QDialog(parent){
			setupUi(this);
			points1_active_ = false;
			points2_active_ = false;
		}

		void paintEvent(QPaintEvent *event);

		void mousePressEvent(QMouseEvent *event);

		void keyPressEvent(QKeyEvent *event);

		void keyReleaseEvent(QKeyEvent *event);

	protected:
		enum PushButton {PointToPointArun, PointToPointUmeyama, PointToPlaneLinear};
		QList<QPoint> points1_;
		QList<QPoint> points2_;
		bool points1_active_;
		bool points2_active_;

	signals:
		void updatePoints(QList<QPoint> &_points1, QList<QPoint> &_points2, PushButton _pushButton);

	private slots:
		void on_pushButton_PointToPointArun_clicked();

		void on_pushButton_PointToPointUmeyama_clicked();

		void on_pushButton_PointToPlaneLinear_clicked();

		void on_pushButton_clear_source_clicked();

		void on_pushButton_clear_target_clicked();

	};
};

#endif