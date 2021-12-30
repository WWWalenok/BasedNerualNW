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
		State = 1,
		Harmonic = 40;

	float
		A[State][Harmonic],
		B[State][Harmonic],
		F[State][Harmonic][Dim];

	uint16_t __Seed = 0;

	float BaseHZ[State]
	{
		0.15
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

	const int SIZE = 256;

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
				M[x][y][s] = (1 / (1 + exp(v * 8))) * 255 * 1.0;
				//M[x][y][s] = 0 + round(tanh((x - SIZE * 0.5) * (x - SIZE * 0.5) * 0.001 + (y - SIZE * 0.5) * (y - SIZE * 0.5) * 0.001) * 255);
			}

		}
	}

	cv::Mat3b m(SIZE * 4, SIZE * 4);

	for (int x = 0; x < SIZE * 4; x++)
	{
		for (int y = 0; y < SIZE * 4; y++)
		{
			for (int s = 0; s < State; s++)
			{
				m[x][y][s] = M[x / 4][y / 4][0];

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

	float**** buff = new float*** [2]{ new float** [SIZE], new float** [SIZE] };


	for (int x = 0; x < SIZE; x++)
	{
		buff[0][x] = new float* [SIZE];
		buff[1][x] = new float* [SIZE];
		for (int y = 0; y < SIZE; y++)
		{
			buff[0][x][y] = new float[State];
			buff[1][x][y] = new float[State];
			for (int s = 0; s < State; s++)
			{
				buff[0][x][y][s] = M[x][y][0];
				buff[1][x][y][s] = M[x][y][0];
			}
		}
	}
	float MaxS[3]
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
		MGX[3][3]
	{
		-0.25, -0.50, -0.25,
		 0.00,  0.00,  0.00,
		+0.25, +0.50, +0.25
	}, 
		MGY[3][3]
	{
		-0.25,  0.00, +0.25,
		-0.50,  0.00, +0.50,
		-0.25,  0.00, +0.25
	},
		MG1[3][3]
	{
		-0.25, -0.50, -0.25,
		+0.25, +0.25, +0.25,
		+0.00, +0.00, +0.00
	},
		MG2[3][3]
	{
		-0.00, -0.00, -0.00,
		-0.25, -0.50, -0.25,
		+0.25, +0.50, +0.25
	},

		ML[3][3]
	{
		__V3__, __V2__, __V3__,
		__V2__, __V0__, __V2__,
		__V3__, __V2__, __V3__
	},
		MF[3][3]
	{
		__V3__, __V2__, __V3__,
		__V2__, -__V0__, __V2__,
		__V3__, __V2__, __V3__
	};

	for (int I = 0; I < 100; I++)
	{
		float MaxSL[3]
		{
			0,
			0,
			0
		};

		cv::Mat3b Test(SIZE * 4, SIZE * 4);
		for (int x = 0; x < SIZE; x++)
		{
			for (int y = 0; y < SIZE; y++)
			{
				int C = 1;

				int
					X[3]
				{
					MAX(0, x - C),
					x,
					MIN(SIZE - 1, x + C)
				},
					Y[3]
				{
					MAX(0, y - C),
					y,
					MIN(SIZE - 1, y + C)
				};

				float GX = 0, GY = 0, L = 0;

				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
					{
						float& T = buff[0][X[i]][Y[j]][0];
						GX += MGX[i][j] * T;
						GY += MGY[i][j] * T;
						L += ML[i][j] * T;
					}

				float G = sqrt(GX * GX + GY * GY);

				float A = acos((G != 0 ? GX / G : 0));

				A = (GY > 0 ? A : 2 * CV_PI - A);

				A = A / CV_PI * 180;

				uint16_t __C__[3]
				{
					A, tanh(G * 1) * 20, buff[0][x][y][0]
				};

				auto Color = hsv2rgb(__C__);
				if(G > 0)
					A = A / CV_PI * 180;

				for (int _x = 0; _x < 4; _x++)
					for (int _y = 0; _y < 4; _y++)
						for (int s = 0; s < 3; s++)
							Test[x * 4 + _x][y * 4 + _y][s] = Color[s];

			}
		}
		float
			Pos[1000][2],
			V[2]
		{
			0, 0
		};
		Pos[0][0] = (rand() * SIZE) / RAND_MAX;
		Pos[0][1] = (rand() * SIZE) / RAND_MAX;
		{
			int 
				x = round(Pos[0][0]),
				y = round(Pos[0][1]);
			while (buff[0][x][y][0] > 50)
			{
				Pos[0][0] = (rand() * SIZE) / RAND_MAX;
				Pos[0][1] = (rand() * SIZE) / RAND_MAX;
				x = round(Pos[0][0]);
				y = round(Pos[0][1]);
			}
		}
		float Map[3][3];

		for (int I = 1; I < 100; I++)
		{
			int C = 1;
			int
				x = round(Pos[I - 1][0]),
				y = round(Pos[I - 1][1]);

			x = MAX(0, MIN(SIZE - 1, x));
			y = MAX(0, MIN(SIZE - 1, y));

			int
				X[3]
			{
				MAX(0, x - C),
				x,
				MIN(SIZE - 1, x + C)
			},
				Y[3]
			{
				MAX(0, y - C),
				y,
				MIN(SIZE - 1, y + C)
			};

			unsigned char Max = buff[0][x][y][0];

			float GX = 0, GY = 0, L = 0;

			float Tx = 0, Ty = 0, S = 0;
			float Min = Max;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					float& T = buff[0][X[i]][Y[j]][0];
					GX += MGX[i][j] * T;
					GY += MGY[i][j] * T;
					L += ML[i][j] * T;				
					Min = MIN(T, Min);
				}

			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
				{
					float T = buff[0][X[i]][Y[j]][0] - Min;
					Tx += X[i] * T;
					Ty += Y[j] * T;
					S += T;
					Map[i][j] = buff[0][X[i]][Y[j]][0];
				}


			Tx = Tx / S - Pos[I - 1][0];
			Ty = Ty / S - Pos[I - 1][1];

			//float TL = sqrt(Tx * Tx + Ty * Ty);

			//Tx /= TL;
			//Ty /= TL;

			V[0] = GX + 0.1 * Tx * 0.5;
			V[1] = GY + 0.1 * Ty * 0.5;

			float VL = sqrt(V[0] * V[0] + V[1] * V[1]);
			if (VL > 1)
			{
				V[0] /= VL;
				V[1] /= VL;
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
					{
						buff[1][X[i]][Y[j]][0] = MAX(MAX(buff[0][X[i]][Y[j]][0] - 20, 0), buff[1][X[i]][Y[j]][0] - 3 * MF[i][j]);;
					}

			}
			else if(VL < 0.01)
			{
				break;
			}

			Pos[I][0] = Pos[I - 1][0] + V[0];
			Pos[I][1] = Pos[I - 1][1] + V[1];

			Pos[I][0] = MAX(0, MIN(SIZE - 1, Pos[I][0]));
			Pos[I][1] = MAX(0, MIN(SIZE - 1, Pos[I][1]));

			uint16_t __C__[3]
			{
				I * 1, 255, 255
			};

			auto Color = hsv2rgb(__C__);

			cv::line(Test, 
				cv::Point((Pos[I - 1][1]) * 4 + 2, (Pos[I - 1][0]) * 4 + 2),
				cv::Point((Pos[I][1]) * 4 + 2, (Pos[I][0]) * 4 + 2),
				cv::Scalar(Color[0], Color[1], Color[2], 0.5));

		}

		for (int s = 0; s < State; s++)
			MaxS[s] = MaxSL[s];
		for (int x = 0; x < SIZE * 4; x++)
		{
			for (int y = 0; y < SIZE * 4; y++)
			{
				m[x][y][0] = buff[1][x / 4][y / 4][0]; // B
				m[x][y][1] = buff[1][x / 4][y / 4][0]; // G
				m[x][y][2] = buff[1][x / 4][y / 4][0]; // R
			}
		}

		cv::imshow("d.jpg", m);
		cv::imshow("Test", Test);
		cv::waitKey(1);






		cv::waitKey(1);
	}

	cv::imwrite("_d.jpg", m);
	cv::imshow("d.jpg", m);
	cv::waitKey(1);

}