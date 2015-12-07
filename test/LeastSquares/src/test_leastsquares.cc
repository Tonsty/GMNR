#include <iostream>

#include "../include/test_leastsquares_dialog_control.h"

int main(int argc, char** argv){
	QApplication app(argc, argv);

	gmnr::TestLeastSquaresDialogControl *dialog = new gmnr::TestLeastSquaresDialogControl;
	dialog->show();
	return app.exec();
}