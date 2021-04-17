#include <opencv2\core\core.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\imgproc\imgproc.hpp>
#include <dcmtk\config\osconfig.h>
#include <dcmtk\dcmdata\dctk.h>
#include <dcmtk\dcmimgle\dcmimage.h>
#include <new>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
#include <array>
#include <omp.h>


std::array<long double, 10>  calcftheta(int vtheta[10][12], int v[3], int ud)
{

	int vx = v[0];
	int vy = v[1];
	int vz = v[2];
	std::array<long double, 10> ftheta = { 0 };

	for (int i = 0; i < 10; i++)
	{

		int offsetx_bb1 = vx + vtheta[i][0];
		int offsety_bb1 = vy + vtheta[i][1];
		int offsetz_bb1 = vz + vtheta[i][2];
		int offsetx_bb2 = vx + vtheta[i][6];
		int offsety_bb2 = vy + vtheta[i][7];
		int offsetz_bb2 = vz + vtheta[i][8];
		int extentx_bb1 = vtheta[i][3];
		int extenty_bb1 = vtheta[i][4];
		int extentz_bb1 = vtheta[i][5];
		int extentx_bb2 = vtheta[i][9];
		int extenty_bb2 = vtheta[i][10];
		int extentz_bb2 = vtheta[i][11];
		if ((offsetx_bb1 >= 0 && offsetx_bb1 <= 249) && (offsety_bb1 >= 0 && offsety_bb1 <= 249) && (offsetz_bb1 >= 0 && offsetz_bb1 <= 107) && ((extentx_bb1 + offsetx_bb1 - 1) >= 0 && (extentx_bb1 + offsetx_bb1 - 1) <= 249) && ((extenty_bb1 + offsety_bb1 - 1) >= 0 && (extenty_bb1 + offsety_bb1 - 1) <= 249) && ((extentz_bb1 + offsetz_bb1 - 1) >= 0 && (extentz_bb1 + offsetz_bb1 - 1) <= 107) && (offsetx_bb2 >= 0 && offsetx_bb2 <= 249) && (offsety_bb2 >= 0 && offsety_bb2 <= 249) && (offsetz_bb2 >= 0 && offsetz_bb2 <= 107) && ((extentx_bb2 + offsetx_bb2 - 1) >= 0 && (extentx_bb2 + offsetx_bb2 - 1) <= 249) && ((extenty_bb2 + offsety_bb2 - 1) >= 0 && (extenty_bb2 + offsety_bb2 - 1) <= 249) && ((extentz_bb2 + offsetz_bb2 - 1) >= 0 && (extentz_bb2 + offsetz_bb2 - 1) <= 107))
		{
			unsigned long int sum_bb1 = 0;

			for (int j = 0; j < extentz_bb1; j++)
			{

				char filename_bb1[100];
				sprintf_s(filename_bb1, "E:\\Patient12_6_2_2_%d\\pt (%d).dcm", ud, offsetz_bb1 + j + 1);

				DicomImage* mem_bb1 = (DicomImage*)kmp_malloc(50 * sizeof(DicomImage(filename_bb1)));
				DicomImage* image_in = new (mem_bb1)DicomImage(filename_bb1);
				int nrows = image_in->getHeight();
				int ncols = image_in->getWidth();
				int nchannels = image_in->getFrameCount();
				int type = CV_MAKETYPE(CV_16U, nchannels);
				ushort* pdata = (ushort*)image_in->getOutputData(16);

				cv::Mat input_bb1;
				input_bb1.create(nrows, ncols, type);
				for (int f = 0; f < nrows; f++)
				{
					ushort* dicom_data = input_bb1.ptr<ushort>(f);
					for (int i = 0; i < ncols; i++)
					{
						dicom_data[i] = pdata[(f * ncols) + i];
					}
				}
				image_in->~DicomImage();
				kmp_free(mem_bb1);

				cv::Mat imageROI_bb1;
				imageROI_bb1 = input_bb1(cv::Rect(offsetx_bb1, offsety_bb1, extentx_bb1, extenty_bb1));
				cv::Scalar total_bb1 = cv::sum(imageROI_bb1);
				sum_bb1 += (unsigned long int)total_bb1[0] + (unsigned long int)total_bb1[1] + (unsigned long int)total_bb1[2];
				input_bb1.release();
				imageROI_bb1.release();
			}

			int pixels_bb1 = extentx_bb1*extenty_bb1*extentz_bb1;
			long double mean_bb1 = (long double)sum_bb1 / (long double)pixels_bb1;

			unsigned long int sum_bb2 = 0;

			for (int k = 0; k < extentz_bb2; k++)
			{

				char filename_bb2[100];
				sprintf_s(filename_bb2, "E:\\Patient12_6_2_2_%d\\pt (%d).dcm", ud, offsetz_bb2 + k + 1);

				DicomImage* mem_bb2 = (DicomImage*)kmp_malloc(50 * sizeof(DicomImage(filename_bb2)));
				DicomImage* image_in_bb2 = new (mem_bb2)DicomImage(filename_bb2);
				//DicomImage* image_in_bb2 = new DicomImage(filename_bb2);
				int nrows_bb2 = image_in_bb2->getHeight();
				int ncols_bb2 = image_in_bb2->getWidth();
				int nchannels_bb2 = image_in_bb2->getFrameCount();
				int type_bb2 = CV_MAKETYPE(CV_16U, nchannels_bb2);
				ushort* pdata_bb2 = (ushort*)image_in_bb2->getOutputData(16);

				cv::Mat input_bb2;
				input_bb2.create(nrows_bb2, ncols_bb2, type_bb2);
				for (int j = 0; j < nrows_bb2; j++)
				{
					ushort* dicom_data_bb2 = input_bb2.ptr<ushort>(j);
					for (int i = 0; i < ncols_bb2; i++)
					{
						dicom_data_bb2[i] = pdata_bb2[(j * ncols_bb2) + i];
					}
				}
				image_in_bb2->~DicomImage();
				kmp_free(mem_bb2);

				cv::Mat imageROI_bb2;
				imageROI_bb2 = input_bb2(cv::Rect(offsetx_bb2, offsety_bb2, extentx_bb2, extenty_bb2));
				cv::Scalar total_bb2 = cv::sum(imageROI_bb2);
				sum_bb2 += (unsigned long int)total_bb2[0] + (unsigned long int)total_bb2[1] + (unsigned long int)total_bb2[2];
				input_bb2.release();
				imageROI_bb2.release();
			}


			int pixels_bb2 = extentx_bb2*extenty_bb2*extentz_bb2;
			long double mean_bb2 = (long double)sum_bb2 / (long double)pixels_bb2;

			ftheta[i] = mean_bb1 - mean_bb2;
		}
		else
		{

			ftheta[i] = 191.919;
		}

	}

	return ftheta;

}


int main()
{

	std::ofstream MyExcelFile0, MyExcelFile1, MyExcelFile2, MyExcelFile3;
	MyExcelFile0.open("E:\\featuretheta_12_0.csv");
	MyExcelFile1.open("E:\\featuretheta_12_1.csv");
	MyExcelFile2.open("E:\\featuretheta_12_2.csv");
	MyExcelFile3.open("E:\\featuretheta_12_3.csv");

	for (int slice = 0; slice < 108; slice = slice + 2)
	{
		for (int rows = 40; rows < 190; rows = rows + 3)
		{
			omp_set_num_threads(4);
#pragma omp parallel
			{
				int id = omp_get_thread_num();
				int nthreads = omp_get_num_threads();
				for (int cols = (3 * id) + 110; cols < 225; cols = cols + 12)
				{
					int bb[6] = { 199, 131, 98, 144, 94, 10 };
					int vx = cols; int vy = rows; int vz = slice;
					int v[3] = { vx, vy, vz };
					int d[6] = { vx - bb[0], vx - bb[1], vy - bb[2], vy - bb[3], vz - bb[4], vz - bb[5] };

					int theta[10][12] = { { 4, -8, 7, 9, 9, 5, -3, 9, 4, 6, 6, 5 },
					{ 7, 5, -6, 7, 9, 4, -9, 4, -4, 5, 8, 5 },
					{ 3, 6, 4, 5, 8, 6, -3, 4, -5, 7, 9, 5 },
					{ 3, 7, 6, 7, 6, 3, -4, -6, -7, 3, 9, 3 },
					{ -5, -9, -4, 4, 7, 7, 7, -6, -4, 3, 6, 6 },
					{ -4, -7, 6, 9, 5, 7, 6, 8, -7, 5, 4, 4 },
					{ 3, 6, -5, 9, 8, 5, 4, -7, -4, 7, 6, 5 },
					{ 9, 8, -3, 8, 7, 4, -9, 6, 5, 10, 4, 6 },
					{ 8, 3, 6, 9, 5, 4, -8, -7, 5, 10, 5, 7 },
					{ -9, -4, 4, 5, 9, 4, -4, -7, 6, 10, 7, 7 } };

					std::array<long double, 10>  feature = calcftheta(theta, v, id);

					if (id == 0)
					{
						for (int c = 0; c < 6; c++)
						{
							MyExcelFile0 << d[c] << ",";
						}
						for (int e = 0; e < 10; e++)
						{
							if (feature[e] == 191.919)
								MyExcelFile0 << "NA" << ",";

							else
								MyExcelFile0 << feature[e] << ",";

						}
						MyExcelFile0 << std::endl;
					}
					else if (id == 1)
					{
						for (int c = 0; c < 6; c++)
						{
							MyExcelFile1 << d[c] << ",";
						}
						for (int e = 0; e < 10; e++)
						{
							if (feature[e] == 191.919)
								MyExcelFile1 << "NA" << ",";

							else
								MyExcelFile1 << feature[e] << ",";

						}
						MyExcelFile1 << std::endl;
					}

					else if (id == 2)
					{
						for (int c = 0; c < 6; c++)
						{
							MyExcelFile2 << d[c] << ",";
						}
						for (int e = 0; e < 10; e++)
						{
							if (feature[e] == 191.919)
								MyExcelFile2 << "NA" << ",";

							else
								MyExcelFile2 << feature[e] << ",";

						}
						MyExcelFile2 << std::endl;
					}
					else
					{
						for (int c = 0; c < 6; c++)
						{
							MyExcelFile3 << d[c] << ",";
						}
						for (int e = 0; e < 10; e++)
						{
							if (feature[e] == 191.919)
								MyExcelFile3 << "NA" << ",";

							else
								MyExcelFile3 << feature[e] << ",";

						}
						MyExcelFile3 << std::endl;
					}



				}
			}
		}

	}

	MyExcelFile0.close();
	MyExcelFile1.close();
	MyExcelFile2.close();
	MyExcelFile3.close();





	return 0;
}
