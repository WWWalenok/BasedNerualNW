#include <iostream>
#include <SFML/Graphics.hpp>
#include<thread>
#include<time.h>
#include <ctime>
#include <ratio>
#include <chrono>
#include <string>
#include <fstream>

float A_E[257];

#define E(x) A_E[(unsigned int)(x + 128) % 257] 

void main()
{
	for (int i = 0; i < 257; i++)
	{
		A_E[i] = expf(i / 256.0 * 4 - 4);
	}
	const int
		Dim = 2,
		State = 3,
		Harmonic = 10;

	float
		A[State][Harmonic],
		K[State][Harmonic][Dim];

	for (int s = 0; s < State; s++)
	{
		float S = 0;
		for (int h = 0; h < Harmonic; h++)
		{
			A[s][h] = rand() / float(RAND_MAX);
			S += A[s][h];
			for (int i = 0; i < Dim; i++)
			{
				K[s][h][i] = (rand() * 0.80 + RAND_MAX * 0.10) / (rand() * 0.1 + RAND_MAX * 0.1);
			}
		}
		for (int h = 0; h < Harmonic; h++)
		{
			A[s][h] = A[s][h] / S * 128.0;
		}
	}

	unsigned char ***M = new unsigned char**[1024];

	std::ofstream out("out.txt");
	std::string str = "           ";
	out << "{\n\t";
	for (int x = 0; x < 1024; x++)
	{
		M[x] = new unsigned char* [1024];
		for (int y = 0; y < 1024; y++)
		{
			M[x][y] = new unsigned char[State];
			out << '{';
			for (int s = 0; s < State; s++)
			{
				float v = 0;
				for (int h = 0; h < Harmonic; h++)
				{
					v += A[s][h] * sinf(x * 0.125 * K[s][h][0] + y * 0.125 * K[s][h][1]);
				}
				M[x][y][s] = 255 * E(v);
			}
			str[3] = str[7] = ',';

			str[0] = M[x][y][0] / 100 + '0';
			str[1] = (M[x][y][0] % 100) / 10 + '0';
			str[2] = M[x][y][0] % 10 + '0';

			str[4] = M[x][y][1] / 100 + '0';
			str[5] = (M[x][y][1] % 100) / 10 + '0';
			str[6] = M[x][y][1] % 10 + '0';

			str[8] = M[x][y][2] / 100 + '0';
			str[9] = (M[x][y][2] % 100) / 10 + '0';
			str[10] = M[x][y][2] % 10 + '0';

			out << str;

			if(y == 1023)
				out << "}";
			else
				out << "},";
		}
		if (x == 1023)
			out << "\n";
		else
			out << ",\n\t";
	}
	out << "}";



	out.close();
}