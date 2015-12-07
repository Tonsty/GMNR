#ifndef GMNR_TEST_TEST_MULTITPS_DIALOG_H
#define GMNR_TEST_TEST_MULTITPS_DIALOG_H

#include <QDialog>
#include <QPoint>

#include <GMNR/common.h>

#include "ui_test_multitps_qt_dialog.h"

class QWidget;
class QPaintEvent;
class QKeyEvent;

namespace gmnr{

	class TestMultiTPSDialog: public QDialog, public Ui::test_multitps{	

		Q_OBJECT

	public:
		TestMultiTPSDialog(QWidget *parent = 0);

		void paintEvent(QPaintEvent *event);

		void mousePressEvent(QMouseEvent *event);

		void keyPressEvent(QKeyEvent *event);

		void keyReleaseEvent(QKeyEvent *event);

	protected:
		void QListPoint_to_PointSet2D(QList<QPoint> &_points, gmnr::PointSet2D &_pointset2d);

		void PointSet2D_to_QListPoint(gmnr::PointSet2D &_pointset2d, QList<QPoint> &_points);

		void setSample();

		enum PushButton {PointToPointArun, PointToPointUmeyama, PointToPlaneLinear};
		
		QList<QPoint> points1_;
		QList<QPoint> points2_;
		QList<QPoint> points3_;

		bool points1_active_;
		bool points2_active_;
		bool points3_active_;

	signals:
		void updatePoints(QList<QPoint> &_points1, QList<QPoint> &_points2, QList<QPoint> &_points3);

	private slots:
		void on_pushButton_multitps_clicked();

		void on_pushButton_clear_1_clicked();

		void on_pushButton_clear_2_clicked();

		void on_pushButton_clear_3_clicked();

		void on_pushButton_reset_clicked();
	};
};

#endif