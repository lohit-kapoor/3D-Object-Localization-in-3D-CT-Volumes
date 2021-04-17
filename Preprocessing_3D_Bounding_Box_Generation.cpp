#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <dcmtk\config\osconfig.h>
#include <dcmtk\dcmdata\dctk.h>
#include <dcmtk\dcmimgle\dcmimage.h>
#include <new>
#include <iostream>
#include <stdio.h>
#include <string>

int main()
{
	int xr[85];
	int xl[85];
	int ya[85];
	int yp[85];
	for (int q = 0; q < 85; q++)
	{
		int z = q + 1;
		char filename[100];
		sprintf_s(filename,"E:\\Patient12_6_2_2_Bin\\BW (%d).dcm", z);
		DicomImage* image_in = new DicomImage(filename);
		int nrows = image_in->getHeight();
		int ncols = image_in->getWidth();
		int nchannels = image_in->getFrameCount();
		int type = CV_MAKETYPE(CV_16U, nchannels);
		ushort* pdata = (ushort*)image_in->getOutputData(16);

		cv::Mat dicom_in;
		dicom_in.create(nrows, ncols, type);
		for (int j = 0; j < nrows; j++)
		{
			ushort* dicom_data = dicom_in.ptr<ushort>(j);
			for (int i = 0; i < ncols; i++)
			{
				dicom_data[i] = pdata[(j * ncols) + i];
			}
		}
		std::cout << *(dicom_in.ptr<ushort>(0)) << std::endl;
		cv::namedWindow("dicom", CV_WINDOW_AUTOSIZE);
		cv::imshow("dicom", dicom_in);
		cv::waitKey(0);
		cv::Mat points;
		cv::Mat dicom_in_mod;
		cv::Mat dicom_in_2 = dicom_in;
		dicom_in_2.convertTo(dicom_in_mod, CV_8UC1);
		cv::findNonZero(dicom_in_mod, points);
		cv::Rect Min_Rect = cv::boundingRect(points);
		cv::rectangle(dicom_in_mod, Min_Rect.tl(), Min_Rect.br(), cv::Scalar(255, 255, 255, 2));
		cv::imshow("Result", dicom_in_mod);
		cv::waitKey(0);
		std::cout << Min_Rect.tl() << std::endl << Min_Rect.br() << std::endl << Min_Rect.x << std::endl << Min_Rect.y << std::endl << Min_Rect.x + Min_Rect.width << std::endl << Min_Rect.y + Min_Rect.height << std::endl;
		xr[q] = Min_Rect.x;
		xl[q] = Min_Rect.x + Min_Rect.width;
		ya[q] = Min_Rect.y;
		yp[q] = Min_Rect.y + Min_Rect.height;
				 
	}
	int tempxl = 0;
	int tempyp = 0;
	int tempxr = xr[0];
	int tempya = ya[0];
	for (int q = 0; q < 85; q++)
		{
			if (xl[q] > tempxl)
				tempxl = xl[q];
			if (yp[q] > tempyp)
				tempyp = yp[q];
			if (xr[q] < tempxr)
				tempxr = xr[q];
			if (ya[q] < tempya)
				tempya = ya[q];

		}
	std::cout << "Largest value of xl is ::" << tempxl << std::endl;
	std::cout << "Largest value of yp is ::" << tempyp << std::endl;
	std::cout << "Minimum value of xr is ::" << tempxr << std::endl;
	std::cout << "Minimum value of ya is ::" << tempya << std::endl;
	
		return 0;
}