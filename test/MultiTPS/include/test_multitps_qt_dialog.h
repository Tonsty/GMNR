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

	protected:
		void QListPoint_to_PointSet2D(QList<QPoint> &_points, gmnr::PointSet2D &_pointset2d);

		void PointSet2D_to_QListPoint(gmnr::PointSet2D &_pointset2d, QList<QPoint> &_points);

		void setDefaultSample();

		void setLastSample();

		void paint(QList<QPoint> &_points1, QList<QPoint> &_points2, QList<QPoint> &_points3, bool _draw_correspondences = false, float opacity = 1.0f);

		enum PushButton {PointToPointArun, PointToPointUmeyama, PointToPlaneLinear};
		
		QList<QPoint> points1_;
		QList<QPoint> points2_;
		QList<QPoint> points3_;

		QList<QPoint> points1_last_;
		QList<QPoint> points2_last_;
		QList<QPoint> points3_last_;

		bool after_multitps;

	signals:
		void updatePoints(QList<QPoint> &_points1, QList<QPoint> &_points2, QList<QPoint> &_points3);

	private slots:
		void on_pushButton_multitps_clicked();

		void on_pushButton_clear_1_clicked();

		void on_pushButton_clear_2_clicked();

		void on_pushButton_clear_3_clicked();

		void on_pushButton_reset_d_clicked();

		void on_pushButton_reset_l_clicked();

		void on_checkBox_showLast_stateChanged(int state);
	};
};

#endif