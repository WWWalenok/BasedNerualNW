#include <iostream>
#include <SFML/Graphics.hpp>
#include<thread>
#include<time.h>
#include <ctime>
#include <ratio>
#include <chrono>

const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

const double A092676[]
{
	1,1,7,127,4369,34807,20036983,2280356863,49020204823,65967241200001
};

const double A092677[]
{
	1, 3, 30, 630, 22680, 178200, 97297200, 10216206000, 198486288000, 237588086736000
};

double inverf(double x)
{
	double tx = x, r = x;
	for (int i = 1; i < 10; ++i)
	{
		tx *= x * x * PI;
		r += tx * A092676[i] / A092677[i];
	}
	return r;
}

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

double Tanh(double x) { return tanh(x); }
double DTanh(double x) { double t = 1 / cosh(x); return t * t; }

double ReLU(double x) { return (x > 0 ? x : 0.01 * x); }
double DReLU(double x) { return (x > 0 ? 1 : 0.01); }


struct NW
{
	double (*F)(double x) = Tanh;
	double (*DF)(double x) =  DTanh;

	double **N, **B, **E, **D, A = 0.002;

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
					//if (rand() % 4 == 0) Links[l][i][o] = 0;
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

double x = 0, y = 0, a = 0, dx = 0, dy = 0, da = 0, ix = 0, iy = 0, a1 = 0, a2 = 0, p1 = 1, p2 = 1;

const int Deep = 4;

int __SS[10]
{
	9,
	4,
	4,
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


double ex = inverf(rand() / float(RAND_MAX) * 2 - 1);
double ey = inverf(rand() / float(RAND_MAX) * 2 - 1) ;


void Train()
{

	if (abs(y) > 3 || abs(x) > 3 || t % 500 == 0)
	{
		x = y = dx = dy = ix = iy = 0;
		A.N[A.L - 2][0] = A.N[A.L - 2][1] = A.N[A.L - 2][2] = A.N[A.L - 2][3] = 0;
	}

	tx = cos(t * 0.00999216122) + cos((t + 12) * 0.0066228562) * 0.5;
	ty = sin(t * 0.0098564215366) + cos((t -651) * 0.0077251612165) * 0.5;
	ex = ex * 0.9 + inverf(rand() / float(RAND_MAX) * 2 - 1) * 0.1;
	ey = ey * 0.9 + inverf(rand() / float(RAND_MAX) * 2 - 1) * 0.1;
	tx += ex * 0.0000333;
	ty += ey * 0.0000333;
	//x = 0;
	//y = 0;
	//A.N[A.L - 2][2] = A.N[A.L - 2][3] = A.N[A.L - 2][2] = A.N[A.L - 2][3] = 0;
	A.N[0][0] = 1;
	A.N[0][1] = tx - x;
	A.N[0][2] = ty - y;
	A.N[0][3] = 0;
	A.N[0][4] = 0;
	A.N[0][5] = A.N[A.L - 2][0];
	A.N[0][6] = A.N[A.L - 2][1];
	A.N[0][7] = A.N[A.L - 2][2];
	A.N[0][8] = A.N[A.L - 2][3];

	A.Upd();

	p1 = A.N[A.L - 1][0];
	p2 = A.N[A.L - 1][1];
	dx += p1 * 0.1;
	dy += p2 * 0.1;
	x += dx * 0.1;
	y += dy * 0.1;
	ix = ix * 0.9 + (tx - x) * 0.1;
	iy = iy * 0.9 + (ty - y) * 0.1;
	A.E[A.L - 1][0] = (tx - x) - dx * dx * dx;
	A.E[A.L - 1][1] = (ty - y) - dy * dy * dy;

	A.Train();

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
					(T < 0 ? 55+200 * abs(T) / (1 + abs(T)) : 55), 
					(T > 0 ? 55+200 * abs(T) / (1 + abs(T)) : 55), 
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
					(T < 0 ? 200: 0), 
					(T > 0 ? 200 : 0), 
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
		if (__GG == 500000)
			delay = 10000;
	}
	isWindowOpen = false;
	return 0;
}
