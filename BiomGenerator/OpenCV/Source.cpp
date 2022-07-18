#include "opencv2/objdetect.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/videoio.hpp"
#include<thread>
#include<time.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

#include "BaseGenerator.h"

float A_E[257];

template<typename HSVType = uint16_t, typename RGBType = uint8_t, HSVType MAX_H = 360/*MAX_H % 6 == 0*/, HSVType MAX_S = 255, HSVType MAX_V = 255, RGBType MAX_RGB = 255 >
RGBType* hsv2rgb(HSVType* in)
{
	RGBType* out = new RGBType[3];
	RGBType
		& R = out[2],
		& G = out[1],
		& B = out[0];
	HSVType
		& H = in[0],
		& S = in[1],
		& V = in[2];

	HSVType
		H_D = MAX_H / 6;

	uint8_t NH = ((H) / (H_D)) % 6;

	HSVType F = (H % MAX_H) - NH * H_D;

	RGBType
		P = (V * MAX_RGB - (V * MAX_RGB * S) / MAX_S) / MAX_V,
		Q = (V * MAX_RGB - (V * MAX_RGB * S * F) / MAX_S / H_D) / MAX_V,
		T = (V * MAX_RGB - (V * MAX_RGB * S * (H_D - F)) / MAX_S / H_D) / MAX_V;

	switch (NH) {
	case 0:
		R = V;
		G = T;
		B = P;
		break;
	case 1:
		R = Q;
		G = V;
		B = P;
		break;
	case 2:
		R = P;
		G = V;
		B = T;
		break;
	case 3:
		R = P;
		G = Q;
		B = V;
		break;
	case 4:
		R = T;
		G = P;
		B = V;
		break;
	case 5:
		R = V;
		G = P;
		B = Q;
		break;
	}

	return out;
}

void main()
{
	const size_t size = 2017;

	BaseGenerator A = BaseGenerator(size);

	cv::Mat1w test(size, size);

	A.Init(10, 40, 0.01, 50, 0, 20);

	{
		float _max = -1E5, _min = -_max;
		for(int i = 0; i < size; i++) for(int j = 0; j < size; j++)
		{
			_max = MAX(_max, A.map[i][j]);
			_min = MIN(_min, A.map[i][j]);
		}

		float a = _min * 0.8 + _max * 0.2;
		float b = _min * 0.2 + _max * 0.8;

		for(int i = 0; i < size; i++) for(int j = 0; j < size; j++)
		{
			A.map[i][j] = fmax(fmin(A.map[i][j], b), a);
		}
	}

	unsigned int max_i = 10;

	

	for(unsigned int i = 0; i < max_i; i++)
	{
		A.SymulateDrop(1200, 80);

		A.Expoze(fmin(30, i * 10));

		float summ = 0;
		for(int i = 0; i < size; i++) for(int j = 0; j < size; j++)
		{
			summ += A.buff[i][j];
		}

		A.Eroze(fmin(75, i * 20));


		float _max = -1E5, _min = -_max;
		for(int i = 0; i < size; i++) for(int j = 0; j < size; j++)
		{
			_max = MAX(_max, A.map[i][j]);
			_min = MIN(_min, A.map[i][j]);
		}
		for(int i = 0; i < size; i++) for(int j = 0; j < size; j++)
		{
			test[i][j] = 100 + (UINT16_MAX - 200) * (A.map[i][j] - _min) / (_max - _min);
		}

		std::cout << i << " : " << summ << " => " << _min << " -- " << _max << "      \r";

		cv::imshow("teat", test);

		cv::imwrite("D:\\!!!\\T\\I_" + std::to_string(i) + ".png", test, { cv::IMWRITE_PXM_BINARY });

		cv::waitKey(1);



	}
}