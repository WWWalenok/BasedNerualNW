﻿// TEster.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>


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
	return r * 0.000149832;
}

 
#define RandGaus inverf(2 * (rand() / float(RAND_MAX)) - 1)


int main()
{
    
	uint32_t M[20003];
	for (int i = 0; i < 20003; i++)
	{
		M[i] = 0;
	}
	int FF = 100000000;
	for (int i = 0; i < FF; i++)
	{
		M[
			uint32_t(std::min(20002.0, std::max(0.0, round(RandGaus * 5000) + 10001)))
		]++;
	}


	int i = 1;
	int C = M[10001];
	for (; i < 10000; i++)
	{
		if (C + M[10001 + i] + M[10001 - i] < FF * 0.97)
		{
			C = C + M[10001 + i] + M[10001 - i];
		}
		else
			break;
	}
	std::cout << i;
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