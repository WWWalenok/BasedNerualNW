﻿// SFML.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <SFML/Graphics.hpp>
#include<thread>
#include<time.h>
#include <ctime>
#include <ratio>
#include <chrono>

void TH_Delay(uint32_t ms)
{
	std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

void Delay(int32_t ms)
{
	if (ms <= 0) return;
	std::thread t(TH_Delay, ms);
	t.join();
}

double Softsign(double x) { return (x > 0 ? x / (1 + x): x / (1 - x)); }
double DSoftsign(double x) { return (x > 0 ? 1 / (1 + 2 * x + x * x) : 1 / (1 - 2 * x + x * x)); }

double ReLU(double x) { return (x > 0 ? x : 0.01 * x); }
double DReLU(double x) { return (x > 0 ? 1 : 0.01); }


double F(double x) { return ReLU(x); }
double DF(double x) { return ReLU(x); }


struct NW
{
	double (*F)(double x) = Softsign;
	double (*DF)(double x) =  DSoftsign;

	double **N, **B, **E, **D, A = 0.0002;

	double ***Links;

	int* S;
	int* M;

	int L;

	int NS;

	NW(int *_S, int _L)
	{
		L = _L;
		S = new int[L];
		M = new int[L];
		NS = 0;
		M[0] = 0;
		for (int l = 0; l < L; NS += S[l], l++, M[l] = NS)
			S[l] = _S[l];

		Links = new double**[L];
		N = new double*[L];
		B = new double*[L];
		E = new double*[L];
		D = new double*[L];

		for (int l = 0; l < L; l++)
		{
			N[l] = new double[S[l]];
			B[l] = new double[S[l]];
			E[l] = new double[S[l]];
			D[l] = new double[S[l]];
			if (l + 1 == L) break;
			Links[l] = new double*[S[l]];
			for (int i = 0; i < S[l]; i++)
			{
				Links[l][i] = new double[S[l + 1]];
				for (int o = 0; o < S[l + 1]; o++)
				{
					Links[l][i][o] = (rand() / float(RAND_MAX) * 2 - 1) * 1;
				}
			}
		}


	}

	void Upd()
	{
		for (int l = 0; l < L - 1; l++)
		{
			for (int o = 0; o < S[l + 1]; o++)
			{
				B[l + 1][o] = 0;
				for (int i = 0; i < S[l]; i++)
					B[l + 1][o] += Links[l][i][o] * N[l][i];
				N[l + 1][o] = F(B[l + 1][o]);
				D[l + 1][o] = DF(B[l + 1][o]);
				E[l + 1][o] = 0;
			}
		}
	}

	void Train()
	{
		for (int l = L - 2; l >= 0; l--)
			for (int i = 0; i < S[l]; i++)			
			{
				E[l][i] = 0;
				for (int o = 0; o < S[l + 1]; o++)
					E[l][i] += Links[l][i][o] * D[l + 1][o] * E[l + 1][o];
			}

		for (int l = L - 2; l >= 0; l--)
			for (int i = 0; i < S[l]; i++)			
				for (int o = 0; o < S[l + 1]; o++)
					Links[l][i][o] += N[l][i] * D[l + 1][o] * E[l + 1][o] * A;

	}
};

uint32_t size = 1000;

double x = 0, y = 0, a = 0, dx = 0, dy = 0, da = 0, a1 = 0, a2 = 0, p1 = 1, p2 = 1;

const int Deep = 5;

int __SS[Deep]
{
	4,
	4,
	3,
	2,
	2
};

sf::CircleShape **Layers;
sf::Vertex ****Links;

NW A(__SS, Deep);

double tx = 0, ty = 1;
uint64_t t = 0;

sf::CircleShape shape(6.f), unit(6.f), Mid(4.f);

sf::Text Outs[4];

double fps = 60;

void Train()
{
	tx = cos(t * 0.00999216122);
	ty = sin(t * 0.0198564215366);
	//x = 0;
	//y = 0;
	A.N[0][0] = tx - x;
	A.N[0][1] = ty - y;
	A.N[0][2] = dx;
	A.N[0][3] = dy;


	A.Upd();

	p1 = A.N[A.L - 1][0];
	p2 = A.N[A.L - 1][1];
	dx += p1 * .1;
	dy += p2 * .1;
	x += dx * 0.1;
	y += dy * 0.1;
	A.E[A.L - 1][0] = (tx - x) - dx * dx * dx * 0.1;
	A.E[A.L - 1][1] = (ty - y) - dy * dy * dy * 0.1;

	A.Train();

	if (abs(y) > 2 || abs(x) > 2)
	{
		x = y = dx = dy = 0;
	}
	t++;
}

void MDraw()
{
	Train();

	for (int i = 0; i < Deep; i++)
		for (int j = 0; j < __SS[i]; j++)
		{
			float T = A.N[i][j];
			Layers[i][j].setFillColor(
				sf::Color(
					(T > 0 ? 55+200 * T / (1 + T) : 55), 
					(T < 0 ? 55-200 * T / (1 - T) : 55), 
					(T == 0 ? 255 : 55)
				) 
			);
		}


	for (int l = 0; l < Deep - 1; l++)
	{
		for (int i = 0; i < __SS[l]; i++)
		{
			for (int o = 0; o < __SS[l + 1]; o++)
			{
				float T = A.Links[l][i][o];
				Links[l][i][o][0].color = Links[l][i][o][1].color = sf::Color(
					(T > 0 ? 200: 0), 
					(T < 0 ? 200 : 0), 
					(T == 0 ? 0 : 0),
					abs(200 * T / (1 + abs(T)))
				);
			}
		}
	}

	Outs[0].setString(std::to_string(p1).c_str());
	Outs[1].setString(std::to_string(p2).c_str());

	shape.setPosition(3 * size / 4.0 - ty * 100, 3 * size / 4.0 + tx * 100);
	unit.setPosition(3 * size / 4.0 - y * 100, 3 * size / 4.0 + x * 100);

}

bool isWindowOpen = true;

uint32_t delay = 1;

void TH_Draw()
{
	while (isWindowOpen)
	{
		//auto t1 = std::chrono::high_resolution_clock::now();
		MDraw();
		std::this_thread::sleep_for(std::chrono::microseconds(delay));
	}
}

int main()
{

	sf::Font font_arial;
	//if (!font_arial.loadFromFile("arial.ttf"));
	//{
		//return 1;
		// error...
	//}

	{
		int Max = 0;
		Layers = new sf::CircleShape*[Deep];
		for (int i = 0; i < Deep; i++)
		{
			Layers[i] = new sf::CircleShape[__SS[i]];
			Max = std::max(Max, __SS[i]);
			for (int j = 0; j < __SS[i]; j++)
			{
				Layers[i][j].setRadius(15);
				Layers[i][j].setFillColor(sf::Color::White);
			}
		}

		for (int i = 0; i < Deep; i++)
		{
			for (int j = 0; j < __SS[i]; j++)
			{
				Layers[i][j].setPosition(i * 80 + 20, 20 + Max * 20 - __SS[i] * 20 + j * 40);
			}
		}

		Links = new sf::Vertex***[Deep];
		for (int l = 0; l < Deep - 1; l++)
		{
			Links[l] = new sf::Vertex**[__SS[l]];
			for (int i = 0; i < __SS[l]; i++)
			{
				Links[l][i] = new sf::Vertex*[__SS[l + 1]];
				for (int o = 0; o < __SS[l + 1]; o++)
				{
					Links[l][i][o] = new sf::Vertex[2];
					Links[l][i][o][0].position = Layers[l][i].getPosition() + sf::Vector2f(15,15);
					Links[l][i][o][1].position = Layers[l + 1][o].getPosition() + sf::Vector2f(15,15);
					Links[l][i][o][0].color = sf::Color::White;
					Links[l][i][o][1].color = sf::Color::White;
				}
			}
		}


	}

    sf::ContextSettings setting;
    setting.antialiasingLevel = 4;
	
	sf::RenderWindow window(sf::VideoMode(size, size), "SFML works!",sf::Style::Default, setting);
   
    shape.setFillColor(sf::Color::Green);
	unit.setFillColor(sf::Color::Red);
	Mid.setFillColor(sf::Color::Blue);

	Mid.setPosition(3 * size / 4.0, 3 * size / 4.0);

	Outs[0].setPosition(0, 0);
	Outs[1].setPosition(0, 40);
	Outs[0].setFillColor(sf::Color::White);
	Outs[1].setFillColor(sf::Color::White);
	Outs[0].setFont(font_arial);
	Outs[0].setString("Asdsds");
	Outs[1].setFont(font_arial);
	Outs[1].setString("Asdsds");

	std::thread draw(TH_Draw);
	double dt = 0;
	auto t1 = std::chrono::high_resolution_clock::now();
	auto t2 = std::chrono::high_resolution_clock::now();
	double adt = -1;
	int __GG = 0;
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();
        window.draw(shape);
		window.draw(unit);
		window.draw(Mid);
		//window.draw(Outs[0]);
		//window.draw(Outs[1]);

		for (int l = 0; l < Deep - 1; l++)
		{
			for (int i = 0; i < __SS[l]; i++)
			{
				for (int o = 0; o < __SS[l + 1]; o++)
				{
					window.draw(Links[l][i][o], 2, sf::Lines, sf::RenderStates::Default);
				}
			}
		}

		for (int i = 0; i < Deep; i++)
		{
			for (int j = 0; j < __SS[i]; j++)
			{
				window.draw(Layers[i][j]);
			}
		}


		t2 = std::chrono::high_resolution_clock::now();
		double dt = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() * 1000;
		t1 = std::chrono::high_resolution_clock::now();
		adt = (adt < 0 ? 4 : adt * 0.9995 + dt * 0.0005);
		std::cout << int(dt * 1000) << " " << int(adt * 1000) << " " << int(1000.0 / adt) << "        \r";
		
		Delay(std::max(1000 / fps - dt, 1.0));
		window.display();
		__GG++;
		if (__GG == 50000)
			delay = 10000;
    }
	isWindowOpen = false;
    return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.