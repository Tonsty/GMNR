#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>

#include <GMNR/ThinPlateSplines.h>
#include <GMNR/DualTPS.h>
#include <GMNR/MultiTPS.h>

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
			cvShowImage("MultiTPS_test", img);
			unsigned char c = cvWaitKey(150);
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
			cvShowImage("MultiTPS_test", img);
			unsigned char c = cvWaitKey(150);
			if (c == 32) break; 
			else if (c == 27) exit(0); 
		}
	}
}

int main(int argc, char** argv){

	std::cout << "float epsilon = " << Eigen::NumTraits<float>::epsilon() << std::endl; 
	std::cout << "double epsilon = " << Eigen::NumTraits<double>::epsilon() << std::endl; 
	std::cout << "float lowest = " << Eigen::NumTraits<float>::lowest() << std::endl; 
	std::cout << "double lowest = " << Eigen::NumTraits<double>::lowest() << std::endl; 

	int kappa1 = 100, kappa2 = 100, kappa3 = 100;
	int lambda1 = 10, lambda2 = 10, lambda3 = 10;
	int kappa1_last = 100, kappa2_last = 100, kappa3_last = 100;
	int lambda1_last = 10, lambda2_last = 10, lambda3_last = 10;

	cvNamedWindow("MultiTPS_test", CV_WINDOW_AUTOSIZE);
	cvCreateTrackbar("kappa1", "MultiTPS_test", &kappa1, 10000);
	cvCreateTrackbar("kappa2", "MultiTPS_test", &kappa2, 10000);
	cvCreateTrackbar("kappa3", "MultiTPS_test", &kappa3, 10000);
	cvCreateTrackbar("lambda1", "MultiTPS_test", &lambda1, 10000);
	cvCreateTrackbar("lambda2", "MultiTPS_test", &lambda2, 10000);
	cvCreateTrackbar("lambda3", "MultiTPS_test", &lambda3, 10000);

	cvSetMouseCallback("MultiTPS_test", my_mouse_callback);

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

	bool draw_1, draw_2, draw_3;
	bool draw_1_more, draw_2_more, draw_3_more;
	bool draw_1_new, draw_2_new, draw_3_new;
	bool draw_1_more_new, draw_2_more_new, draw_3_more_new;

	draw_1 = true;
	draw_2 = true;
	draw_3 = true;

	draw_1_more = true;
	draw_2_more = true;
	draw_3_more = true;

	//draw_1 = false; 
	//draw_2 = false; 
	//draw_3 = false; 

	//draw_1_more = false;
	//draw_2_more = false;
	//draw_3_more = false;

	draw_1_new = true; 
	draw_2_new = true; 
	draw_3_new = true; 

	draw_1_more_new = true;
	draw_2_more_new = true;
	draw_3_more_new = true;

	//draw_1_new = false; 
	//draw_2_new = false; 
	//draw_3_new = false; 

	//draw_1_more_new = false;
	//draw_2_more_new = false;
	//draw_3_more_new = false;


	//------------------------------
	Matrix X1(11, 1);
	for (int i = 0; i < X1.size(); i++) X1(i) = 0.2f * i - 1.0f;
	Matrix Y1 = 0.1f * Matrix::Random(11, 1);
	if(draw_1){
		color.val[0] = 0.0;
		color.val[1] = 0.0;
		color.val[2] = 255.0;
		color.val[3] = 0.0;
		drawCircles(X1, Y1, color, img, true);
	}

	TPSFunction tps(X1, Y1, 1.e-6f);

	Matrix X1_more(801, 1);
	for (int i = 0; i < X1_more.size(); i++) X1_more(i) = 0.0025f * i - 1.0f;
	Matrix Y1_more = tps.evaluate(X1_more);
	if ( draw_1_more )
	{
		color.val[0] = 0.0;
		color.val[1] = 0.0;
		color.val[2] = 255.0;
		color.val[3] = 0.0;
		drawLine(X1_more, Y1_more, color, img, true);
	}

	Matrix X2 = X1 + 0.05f * Matrix::Random(11, 1);
	Matrix Y2 = Y1 + 0.2f * Matrix::Random(11, 1) + 0.3f * Matrix::Ones(11, 1);
	//Matrix Y2 = Y1;
	Y2(0) = Y1(0) - 0.05f;
	Y2(1) = Y1(1) - 0.05f;
	Y2(9) = Y1(9) - 0.04f;
	Y2(10) = Y1(10) - 0.04f;
	if (draw_2)
	{
		color.val[0] = 255.0;
		color.val[1] = 0.0;
		color.val[2] = 0.0;
		color.val[3] = 0.0;
		drawCircles(X2, Y2, color, img, true);
	}

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
	if (draw_2_more)
	{
		color.val[0] = 255.0;
		color.val[1] = 0.0;
		color.val[2] = 0.0;
		color.val[3] = 0.0;
		drawLine(X2Y2_more.col(0), X2Y2_more.col(1), color, img, true);
	}

	Matrix X3 = X1 + 0.05f * Matrix::Random(11, 1);
	Matrix Y3 = Y1 + 0.05f * Matrix::Random(11, 1) + 0.2f * Matrix::Ones(11, 1);
	//Matrix Y3 = Y1;
	//Y3(0) = Y1(0) - 0.01f;
	//Y3(1) = Y1(1) - 0.02f;
	//Y3(9) = Y1(9) - 0.01f;
	//Y3(10) = Y1(10) - 0.02f;
	if (draw_3)
	{
		color.val[0] = 125.0;
		color.val[1] = 125.0;
		color.val[2] = 0.0;
		color.val[3] = 0.0;
		drawCircles(X3, Y3, color, img, true);
	}

	Matrix X3Y3(11, 2);
	X3Y3.col(0) = X3;
	X3Y3.col(1) = Y3;
	TPSFunction tps2(X1Y1, X3Y3, 0.00001f);
	Matrix X3Y3_more(801, 2);
	X3Y3_more = tps2.evaluate(X1Y1_more);
	if (draw_3_more)
	{
		color.val[0] = 125.0;
		color.val[1] = 125.0;
		color.val[2] = 0.0;
		color.val[3] = 0.0;
		drawLine(X3Y3_more.col(0), X3Y3_more.col(1), color, img, true);
	}

	//{
	//	//DualTPS dtps(X1Y1, X2Y2, 0.01, 0.01, 0.001, 0.001);
	//	//std::cout << "dtps = \n" << dtps << std::endl;
	//	//TPSFunction f1 = dtps.getf1(), f2 = dtps.getf2();

	//	Matrix X = X1Y1, Y = X2Y2;
	//	std::vector<int> m(1, X.rows());
	//	std::vector<int> alpha(1, 0);
	//	std::vector<int> beta(1, 1);
	//	Vector kappa(2);
	//	kappa << 0.01, 0.01;
	//	Vector lambda(2);
	//	lambda << 0.001, 0.001;
	//	MultiTPS mtps(X, Y, m, alpha, beta, kappa, lambda);
	//	TPSFunction f1 = mtps.getfs()[0], f2 = mtps.getfs()[1];

	//	CvScalar white;
	//	white.val[0] = 255.0;
	//	white.val[1] = 255.0;
	//	white.val[2] = 255.0;
	//	white.val[3] = 0.0;
	//	cvSet(img, white);
	//	cvShowImage("MultiTPS_test", img);

	//	Matrix new1 = f1.evaluate(X1Y1);
	//	if (draw_1_new)
	//	{
	//		color.val[0] = 0.0;
	//		color.val[1] = 0.0;
	//		color.val[2] = 250.0;
	//		color.val[3] = 0.0;
	//		drawCircles(new1.col(0), new1.col(1), color, img, true);
	//	}

	//	Matrix new2 = f2.evaluate(X2Y2);
	//	if (draw_2_new)
	//	{
	//		color.val[0] = 255.0;
	//		color.val[1] = 0.0;
	//		color.val[2] = 0.0;
	//		color.val[3] = 0.0;
	//		drawCircles(new2.col(0), new2.col(1), color, img, true);
	//	}

	//	std::cout << "tps_energy1_base =" << (new1 - new2).squaredNorm() << std::endl;
	//	std::cout << "tps_energy2_base =" << kappa(0) * (new1 - X1Y1).squaredNorm() << std::endl;
	//	std::cout << "tps_energy3_base =" << kappa(1) * (new2 - X2Y2).squaredNorm() << std::endl;
	//	std::cout << "tps_energy4_base =" << lambda(0) * f1.getA().col(0).transpose() * greenFunc(f1.getX(), f1.getX()) * f1.getA().col(0) +
	//		lambda(0) * f1.getA().col(1).transpose() * greenFunc(f1.getX(), f1.getX()) * f1.getA().col(1) +
	//		lambda(1) * f2.getA().col(0).transpose() * greenFunc(f2.getX(), f2.getX()) * f2.getA().col(0) + 
	//		lambda(1) * f2.getA().col(1).transpose() * greenFunc(f2.getX(), f2.getX()) * f2.getA().col(1) << std::endl;

	//	Matrix new1_more = f1.evaluate(X1Y1_more);
	//	if (draw_1_more_new)
	//	{
	//		color.val[0] = 0.0;
	//		color.val[1] = 0.0;
	//		color.val[2] = 255.0;
	//		color.val[3] = 0.0;
	//		drawLine(new1_more.col(0), new1_more.col(1), color, img, true);
	//	}

	//	Matrix new2_more = f2.evaluate(X2Y2_more);
	//	if (draw_2_more_new)
	//	{
	//		color.val[0] = 255.0;
	//		color.val[1] = 0.0;
	//		color.val[2] = 0.0;
	//		color.val[3] = 0.0;
	//		drawLine(new2_more.col(0), new2_more.col(1), color, img, true);		
	//	}
	//}

	//-------------------------------------------------------------------------------------------------

	Matrix X(X1Y1.rows() + X2Y2.rows(), X1Y1.cols()), Y(X2Y2.rows() + X2Y2.rows(), X2Y2.cols());
	X << X1Y1, X2Y2;
	Y << X2Y2, X3Y3;
	std::vector<int> m(2);
	m[0] = X1Y1.rows();
	m[1] = X2Y2.rows();
	std::vector<int> alpha(2);
	alpha[0] = 0;
	alpha[1] = 1;
	std::vector<int> beta(2);
	beta[0] = 1;
	beta[1] = 2;
	Vector kappa(3);
	kappa << 0.01, 0.01, 0.01;
	Vector lambda(3);
	lambda << 0.001, 0.001, 0.001;

	MultiTPS mtps(X, Y, m, alpha, beta, kappa, lambda);
	//std::cout << "mtps = \n" << mtps << std::endl;
	TPSFunction f1 = mtps.getfs()[0], f2 = mtps.getfs()[1], f3 = mtps.getfs()[2];

	Scalar XtA_norm_1 = (f1.getX().transpose() * f1.getA()).norm();
	std::cout << "XtA_norm_1 is: " << XtA_norm_1 << std::endl;

	Scalar XtA_norm_2 = (f2.getX().transpose() * f2.getA()).norm();
	std::cout << "XtA_norm_2 is: " << XtA_norm_2 << std::endl;

	Scalar XtA_norm_3 = (f1.getX().transpose() * f3.getA()).norm();
	std::cout << "XtA_norm_3 is: " << XtA_norm_1 << std::endl;

	CvScalar white;
	white.val[0] = 255.0;
	white.val[1] = 255.0;
	white.val[2] = 255.0;
	white.val[3] = 0.0;
	cvSet(img, white);
	cvShowImage("MultiTPS_test", img);

	Matrix new1 = f1.evaluate(X1Y1);
	if (draw_1_new)
	{
		color.val[0] = 0.0;
		color.val[1] = 0.0;
		color.val[2] = 250.0;
		color.val[3] = 0.0;
		drawCircles(new1.col(0), new1.col(1), color, img, true);
	}

	Matrix new2 = f2.evaluate(X2Y2);
	if (draw_2_new)
	{
		color.val[0] = 255.0;
		color.val[1] = 0.0;
		color.val[2] = 0.0;
		color.val[3] = 0.0;
		drawCircles(new2.col(0), new2.col(1), color, img, true);
	}

	Matrix new3 = f3.evaluate(X3Y3);
	if (draw_3_new)
	{
		color.val[0] = 125.0;
		color.val[1] = 125.0;
		color.val[2] = 0.0;
		color.val[3] = 0.0;
		drawCircles(new3.col(0), new3.col(1), color, img, true);
	}

	Matrix new1_more = f1.evaluate(X1Y1_more);
	if (draw_1_more_new)
	{
		color.val[0] = 0.0;
		color.val[1] = 0.0;
		color.val[2] = 255.0;
		color.val[3] = 0.0;
		drawLine(new1_more.col(0), new1_more.col(1), color, img, true);
	}

	Matrix new2_more = f2.evaluate(X2Y2_more);
	if (draw_2_more_new)
	{
		color.val[0] = 255.0;
		color.val[1] = 0.0;
		color.val[2] = 0.0;
		color.val[3] = 0.0;
		drawLine(new2_more.col(0), new2_more.col(1), color, img, true);		
	}

	Matrix new3_more = f3.evaluate(X3Y3_more);
	if (draw_3_more_new)
	{
		color.val[0] = 125.0;
		color.val[1] = 125.0;
		color.val[2] = 0.0;
		color.val[3] = 0.0;
		drawLine(new3_more.col(0), new3_more.col(1), color, img, true);
	}

	color.val[0] = 0.0;
	color.val[1] = 0.0;
	color.val[2] = 250.0;
	color.val[3] = 0.0;
	drawCircles(X1, Y1, color, img);
	drawLine(X1Y1_more.col(0), X1Y1_more.col(1), color, img);

	color.val[0] = 255.0;
	color.val[1] = 0.0;
	color.val[2] = 0.0;
	color.val[3] = 0.0;
	drawCircles(X2, Y2, color, img);
	drawLine(X2Y2_more.col(0), X2Y2_more.col(1), color, img);

	color.val[0] = 125.0;
	color.val[1] = 125.0;
	color.val[2] = 0.0;
	color.val[3] = 0.0;
	drawCircles(X3, Y3, color, img);
	drawLine(X3Y3_more.col(0), X3Y3_more.col(1), color, img);

	while(1){
		if ( kappa1 != kappa1_last  ||
			kappa2 != kappa2_last ||
			kappa3 != kappa3_last ||
			lambda1 != lambda1_last ||
			lambda2 != lambda2_last ||
			lambda3 != lambda3_last){

				kappa1_last = kappa1;
				kappa2_last = kappa2;
				kappa3_last = kappa3;
				lambda1_last = lambda1;
				lambda2_last = lambda2;
				lambda3_last = lambda3;

				CvScalar color;
				color.val[0] = 255.0;
				color.val[1] = 255.0;
				color.val[2] = 255.0;
				color.val[3] = 0.0;
				cvSet(img, color);

				color.val[0] = 0.0;
				color.val[1] = 0.0;
				color.val[2] = 250.0;
				color.val[3] = 0.0;
				drawCircles(X1, Y1, color, img);
				drawLine(X1Y1_more.col(0), X1Y1_more.col(1), color, img);

				color.val[0] = 255.0;
				color.val[1] = 0.0;
				color.val[2] = 0.0;
				color.val[3] = 0.0;
				drawCircles(X2, Y2, color, img);
				drawLine(X2Y2_more.col(0), X2Y2_more.col(1), color, img);

				color.val[0] = 125.0;
				color.val[1] = 125.0;
				color.val[2] = 0.0;
				color.val[3] = 0.0;
				drawCircles(X3, Y3, color, img);
				drawLine(X3Y3_more.col(0), X3Y3_more.col(1), color, img);

				Vector kappa(3);
				kappa << (kappa1+1.0f)/10000, (kappa2+1.0f)/10000, (kappa3+1.0f)/10000;
				Vector lambda(3);
				lambda << (lambda1+1.0f)/10000, (lambda1+1.0f)/10000, (lambda1+1.0f)/10000;

				MultiTPS mtps(X, Y, m, alpha, beta, kappa, lambda);

				new1 = mtps.getfs()[0].evaluate(X1Y1);
				color.val[0] = 0.0;
				color.val[1] = 0.0;
				color.val[2] = 250.0;
				color.val[3] = 0.0;
				drawCircles(new1.col(0), new1.col(1), color, img);

				new2 = mtps.getfs()[1].evaluate(X2Y2);
				color.val[0] = 255.0;
				color.val[1] = 0.0;
				color.val[2] = 0.0;
				color.val[3] = 0.0;
				drawCircles(new2.col(0), new2.col(1), color, img);

				new3 = mtps.getfs()[2].evaluate(X3Y3);
				color.val[0] = 125.0;
				color.val[1] = 125.0;
				color.val[2] = 0.0;
				color.val[3] = 0.0;
				drawCircles(new2.col(0), new2.col(1), color, img);

				new1_more = mtps.getfs()[0].evaluate(X1Y1_more);
				color.val[0] = 0.0;
				color.val[1] = 0.0;
				color.val[2] = 250.0;
				color.val[3] = 0.0;
				drawLine(new1_more.col(0), new1_more.col(1), color, img);

				new2_more = mtps.getfs()[1].evaluate(X2Y2_more);
				color.val[0] = 255.0;
				color.val[1] = 0.0;
				color.val[2] = 0.0;
				color.val[3] = 0.0;
				drawLine(new2_more.col(0), new2_more.col(1), color, img);

				new3_more = mtps.getfs()[2].evaluate(X3Y3_more);
				color.val[0] = 125.0;
				color.val[1] = 125.0;
				color.val[2] = 0.0;
				color.val[3] = 0.0;
				drawLine(new3_more.col(0), new3_more.col(1), color, img);
		}
		cvShowImage("MultiTPS_test", img);
		unsigned char c = cvWaitKey(150);
		if (c == 27) return 0;
	}

	while( (unsigned char)cvWaitKey(150) != 27 );

 	return 0;
}