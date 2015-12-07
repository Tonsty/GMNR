#ifndef GMNR_TEST_TEST_LEASTSQUARES_CONTROL_H
#define GMNR_TEST_TEST_LEASTSQUARES_CONTROL_H

#include <QList>
#include <QPoint>

#include "test_leastsquares_dialog.h"
#include <GMNR/common.h>

class QWidget;

namespace gmnr{

	class TestLeastSquaresDialogControl: public TestLeastSquaresDialog{	

		Q_OBJECT

	public:
		TestLeastSquaresDialogControl(QWidget *parent = 0): TestLeastSquaresDialog(parent){
			connect(this, SIGNAL(updatePoints(QList<QPoint>&, QList<QPoint>&, PushButton)), 
				this, SLOT(on_dialog_updatePoints(QList<QPoint>&, QList<QPoint>&, PushButton)));
		}

	private:
		void QListPoint_to_PointSet3D(QList<QPoint> &_points, gmnr::PointSet3D &_pointset3d);

		void PointSet3D_to_QListPoint(gmnr::PointSet3D &_pointset3d, QList<QPoint> &_points);

		void QListPoint_to_NormalSet3D(QList<QPoint> &_points, gmnr::NormalSet3D &_normalset3d);

		private slots:
			void on_dialog_updatePoints(QList<QPoint> &_points1, QList<QPoint> &_points2, PushButton _pushButton);
	};
};

#endif