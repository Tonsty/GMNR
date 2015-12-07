#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

#include <GMNR/ThinPlateSplines.h>

using namespace gmnr;

void my_mouse_callback(int event, int x, int y, int flags, void* param){

}

void drawCircles(const Matrix &X, const Matrix &Y, const CvScalar &color, IplImage* img, bool wait_and_show = false){

	int num = X.rows();

	CvPoint **pts = new CvPoint*[1];
	pts[0] = new CvPoint[num];
	int *npts = new int[1];
	npts[0] = num;

	for (int i = 0; i < num; i++)
	{
		(*pts)[i].x = (int)((X(i) + 1.0f) * 400);
		(*pts)[i].y = (int)((Y(i) + 1.0f) * 400);
	}
	//cvPolyLine(img, pts, npts, 1, 0, color, 1);
	for (int i = 0; i < num; i++)
	{
		cvCircle(img, (*pts)[i], 5, color);
	}
	if (wait_and_show)
	{
		while(1){
			cvShowImage("TPS_test", img);
			unsigned char c = cvWaitKey(15);
			if (c == 32) break; 
			else if (c == 27) exit(0); 
		}
	}
}

void drawLine(const Vector &X_more, const Vector &Y_more, const CvScalar &color, IplImage* img, bool wait_and_show = false){
	int num = X_more.rows();

	CvPoint **pts = new CvPoint*[1];
	pts[0] = new CvPoint[num];
	int *npts = new int[1];
	npts[0] = num;

	for (int i = 0; i < num; i++)
	{
		(*pts)[i].x = (int)((X_more(i) + 1.0f) * 400);
		(*pts)[i].y = (int)((Y_more(i) + 1.0f) * 400);
	}
	cvPolyLine(img, pts, npts, 1, 0, color, 1);
	if (wait_and_show)
	{
		while(1){
			cvShowImage("TPS_test", img);
			unsigned char c = cvWaitKey(15);
			if (c == 32) break; 
			else if (c == 27) exit(0); 
		}
	}
}

int main(int argc, char** argv){

	int kappa1 = 100, kappa2 = 100, lambda1 = 10, lambda2 = 10;
	int kappa1_last = 100, kappa2_last = 100, lambda1_last = 10, lambda2_last = 10;

	cvNamedWindow("TPS_test", CV_WINDOW_AUTOSIZE);

	cvSetMouseCallback("TPS_test", my_mouse_callback);
	
	CvSize size;
	size.height = 800;
	size.width = 800;
	IplImage* img = cvCreateImage(size, IPL_DEPTH_8U, 3);
	cvZero(img);
	CvScalar color;
	color.val[0] = 255.0;
	color.val[1] = 255.0;
	color.val[2] = 255.0;
	color.val[3] = 0.0;
	cvSet(img, color);

//------------------------------
	Matrix X1(11, 1);
	for (int i = 0; i < X1.size(); i++) X1(i) = 0.2f * i - 1.0f;
	Matrix Y1 = 0.1f * Vector::Random(11, 1);
	//std::cout << "X1 = \n" << X1 << std::endl;
	//std::cout << "Y1 = \n" << Y1 << std::endl;

	color.val[0] = 0.0;
	color.val[1] = 0.0;
	color.val[2] = 255.0;
	color.val[3] = 0.0;
	drawCircles(X1, Y1, color, img, true);

	TPSFunction tps(X1, Y1, 1.e-6f);
	//std::cout << "tps = \n" << tps << std::endl;
	//std::cout << "tps(X) = \n" << tps.evaluate(X1) << std::endl;
	//std::cout << "tps(X) - Y = \n" << (tps.evaluate(X1) - Y1) << std::endl;

	Matrix X1_more(801, 1);
	for (int i = 0; i < X1_more.size(); i++) X1_more(i) = 0.0025f * i - 1.0f;
	Matrix Y1_more = tps.evaluate(X1_more);

	color.val[0] = 0.0;
	color.val[1] = 0.0;
	color.val[2] = 255.0;
	color.val[3] = 0.0;
	drawLine(X1_more, Y1_more, color, img, true);

	Matrix X2 = X1 + 0.05f * Matrix::Random(11, 1);
	Matrix Y2 = Y1 + 0.2f * Matrix::Random(11, 1) + 0.3f * Matrix::Ones(11, 1);
	//Matrix Y2 = Y1;
	Y2(0) = Y1(0) - 0.05f;
	Y2(1) = Y1(1) - 0.05f;
	Y2(9) = Y1(9) - 0.04f;
	Y2(10) = Y1(10) - 0.04f;

	color.val[0] = 255.0;
	color.val[1] = 0.0;
	color.val[2] = 0.0;
	color.val[3] = 0.0;
	drawCircles(X2, Y2, color, img, true);

 	Matrix X1Y1(11,2), X2Y2(11, 2);
	X1Y1.col(0) = X1;
	X1Y1.col(1) = Y1;
	X2Y2.col(0) = X2;
	X2Y2.col(1) = Y2;

  	TPSFunction tps1(X1Y1, X2Y2, 0.00001f);
	Matrix X1Y1_more(801, 2), X2Y2_more(801, 2);
	X1Y1_more.col(0) = X1_more;
	X1Y1_more.col(1) = Y1_more;
	X2Y2_more = tps1.evaluate(X1Y1_more);

	color.val[0] = 255.0;
	color.val[1] = 0.0;
	color.val[2] = 0.0;
	color.val[3] = 0.0;
	drawLine(X2Y2_more.col(0), X2Y2_more.col(1), color, img, true);

	while( (unsigned char)cvWaitKey(15) != 27 );

	return 0;
}