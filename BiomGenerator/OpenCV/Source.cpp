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

float A_E[257];

uint16_t xorshift16_xy(uint16_t x, uint16_t y)
{
	y = ((y & 0xff00) >> 8) ^ ((y & 0x00ff) << 8);
	x = ((x & 0xf0f0) >> 4) ^ ((x & 0x0f0f) << 4);
	uint16_t t = x ^ y;
	t = ((t & 0xfff0) >> 4) ^ ((t & 0x000f) << 12);	
	t ^= t >> 3;
	t ^= t << 5;
	uint16_t u = (x ^ y);
	return (y ^ (u << 3)) ^ t;
}

#define E(x) A_E[(unsigned int)(x + 128) % 257] 

const unsigned int randbytelenght = 5333;

uint8_t randbyte[randbytelenght];

uint16_t randdword[randbytelenght];


const int Size = 4;

inline
float _F(float x)
{
	float T = -x / (0.5 + x * x / 200.0) + x;

	return
		T;
}

void main()
{
	for (int i = 0; i < randbytelenght; i++)
	{
		randbyte[i] = rand();
		randdword[i] = rand();
	}
	for (int i = 0; i < 257; i++)
	{
		A_E[i] = expf(i / 256.0 * 3 - 3);
	}
	const int
		Dim = 2,
		State = 3,
		Harmonic = 40;

	float
		A[State][Harmonic],
		B[State][Harmonic],
		F[State][Harmonic][Dim];

	uint16_t __Seed = 0;

	float BaseHZ[State]
	{
		0.02,
		0.05,
		0.005
	};

	for (int s = 0; s < State; s++)
	{
		float S = 0;
	
		for (int h = 0; h < Harmonic; h++)
		{
			float a = rand() / float(RAND_MAX) * 6.28318530718;
			float waveLenght = 4 + (h / 4) * 4;
			float Hz = BaseHZ[s] * (1 + ((h / 4) * 0.25));
			{
				float T = 0;
				for (int i = 0; i < 5; i++)
				{
					T += randdword[__Seed = (__Seed + 1 == randbytelenght ? 0 : __Seed + 1)] / float(UINT16_MAX * 0.5) - 1;
				}
				T = T / 5.0 * 0.5;
				Hz = Hz * (1 + T);
			}
			waveLenght = 1 / Hz;
			A[s][h] = 1.0 * pow(waveLenght, 1.5);
			B[s][h] = rand() / float(RAND_MAX) * 6.3;
			S += A[s][h];

			F[s][h][0] = cos(a) * Hz;
			F[s][h][1] = sin(a) * Hz;
		}
		for (int h = 0; h < Harmonic; h++)
		{
			A[s][h] = A[s][h] / S;
		}
	}

	const int SIZE = 1024;

	unsigned char*** M = new unsigned char** [SIZE];


	for (int x = 0; x < SIZE; x++)
	{
		M[x] = new unsigned char* [SIZE];
		for (int y = 0; y < SIZE; y++)
		{
			M[x][y] = new unsigned char[State];
			for (int s = 0; s < State; s++)
			{
				float v = 0;
				for (int h = 0; h < Harmonic; h++)
				{
					v += A[s][h] * sinf(x * F[s][h][0] + y * F[s][h][1] + B[s][h]);
				}
				M[x][y][s] = powf(4, v - 1) * 255 * 1.0;
			}

		}
	}

	cv::Mat3b m(SIZE, SIZE);

	for (int x = 0; x < SIZE; x++)
	{
		for (int y = 0; y < SIZE; y++)
		{
			for (int s = 0; s < State; s++)
			{
				m[x][y][s] = M[x][y][s];

			}
		}
	}

	cv::imwrite("_t.jpg", m);
	//cv::imshow("t.jpg", m);
	//cv::waitKey();

	const int
		_S1 = 8, 
		_S2 = 2,
		_S0 = _S1 * _S2;

	float _SM[_S0][2];

	for (int i = 0; i < _S2; i++)
	{
		for (int j = 0; j < _S1; j++)
		{
			float a = j / float(_S1) * CV_PI * 2;
			_SM[i * _S1 + j][0] = cos(a) * (i + 1);
			_SM[i * _S1 + j][1] = sin(a) * (i + 1);
		}
	}

	float TG = 4;

	float**** buff = new float***[2]{ new float** [SIZE], new float** [SIZE] };


	for (int x = 0; x < SIZE; x++)
	{
		buff[0][x] = new float* [SIZE];
		buff[1][x] = new float* [SIZE];
		for (int y = 0; y < SIZE; y++)
		{
			buff[0][x][y] = new float [State];
			buff[1][x][y] = new float [State];
			for (int s = 0; s < State; s++)
			{
				buff[0][x][y][s] = M[x][y][s];
				buff[1][x][y][s] = M[x][y][s];
			}
		}
	}
	float MaxS[3]
	{
		255,
		255,
		255
	};
	for (int I = 0; I < 10; I++)
	{
		float MaxSL[3]
		{
			0,
			0,
			0
		};
		{
			auto BT = buff[0];
			buff[0] = buff[1];
			buff[1] = BT;
		}
		float
			D = 1
			;

		unsigned char Mask[3][3];

		const float 
			__V0__ = -1,
			__V1__ = 4 + 4 / sqrt(2),
			__V2__ = 1 / __V1__,
			__V3__ = 1 / sqrt(2) / __V1__;

		float
			MG[3][3]
		{
			-0.25, -0.50, -0.25,
			 0.00,  0.00,  0.00,
			+0.25, +0.50, +0.25
		},
			ML[3][3]
		{
			__V3__, __V2__, __V3__,
			__V2__, __V0__, __V2__,
			__V3__, __V2__, __V3__
		};

		for (int x = 0; x < SIZE; x++)
		{
			for (int y = 0; y < SIZE; y++)
			{
				for (int s = 0; s < State; s++)
				{
					int
						X[3]
					{
						MAX(0, x - 1),
						x,
						MIN(SIZE - 1, x + 1)
					},
						Y[3]
					{
						MAX(0, y - 1),
						y,
						MIN(SIZE - 1, y + 1)
					};

					unsigned char Max = buff[0][x][y][s];

					float GX = 0, GY = 0, L = 0;

					for(int i = 0; i < 3; i++)
						for (int j = 0; j < 3; j++)
						{
							float& T = buff[0][X[i]][Y[j]][s];
							GX += MG[i][j] * T;
							GY += MG[j][i] * T;
							L += ML[i][j] * T;
							Max = MAX(Max, T);
						}

					float G = 0.25 * sqrt(GX * GX + GY * GY);

					if (abs(L) < 1)
						L = 0;

					MaxSL[s] = MAX(Max, MaxSL[s]);

					if (Max > buff[0][x][y][s] || Max < MaxS[s] * 0.0)
					{
						float T = buff[0][x][y][s] - 0 * (TG - G) + 0.5 * L;
						buff[1][x][y][s] =
							MAX(0, MIN(255, T));
						;
					}
					else
						buff[1][x][y][s] = buff[0][x][y][s];
				}
			}
		}
		for (int s = 0; s < State; s++)
			MaxS[s] = MaxSL[s];
		for (int x = 0; x < SIZE; x++)
		{
			for (int y = 0; y < SIZE; y++)
			{
				for (int s = 0; s < State; s++)
				{;
					m[x][y][s] = buff[1][x][y][s];
				}
			}
		}

		cv::imshow("d.jpg", m);
		cv::waitKey(1);
	}

	cv::imwrite("_d.jpg", m);
	cv::imshow("d.jpg", m);
	cv::waitKey(1);

}