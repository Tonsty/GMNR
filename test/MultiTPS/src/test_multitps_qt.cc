#include <iostream>

#include "../include/test_multitps_qt_dialog_control.h"

int main(int argc, char** argv){
	QApplication app(argc, argv);

	gmnr::TestMultiTPSDialogControl *dialog = new gmnr::TestMultiTPSDialogControl;
	dialog->setStyleSheet("background-color:white");
	dialog->show();
	return app.exec();
}