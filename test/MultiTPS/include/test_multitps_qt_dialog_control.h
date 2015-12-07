#ifndef GMNR_TEST_TEST_MULTITPS_CONTROL_H
#define GMNR_TEST_TEST_MULTITPS_CONTROL_H

#include <QList>
#include <QPoint>

#include "test_multitps_qt_dialog.h"
#include <GMNR/common.h>

class QWidget;

namespace gmnr{

	class TestMultiTPSDialogControl: public TestMultiTPSDialog{	

		Q_OBJECT

	public:
		TestMultiTPSDialogControl(QWidget *parent = 0): TestMultiTPSDialog(parent){
			connect(this, SIGNAL(updatePoints(QList<QPoint>&, QList<QPoint>&, QList<QPoint>&)), 
				this, SLOT(on_dialog_updatePoints(QList<QPoint>&, QList<QPoint>&, QList<QPoint>&)));
		}

	protected:

		private slots:
			void on_dialog_updatePoints(QList<QPoint> &_points1, QList<QPoint> &_points2, QList<QPoint>& _points3);
	};
};

#endif