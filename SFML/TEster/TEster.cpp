// TEster.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <chrono>
#include<string>

const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

const double IerfCC[]
{
	1.0, 
	0.0, 
	1.0 / 3.0, 
	0.0, 
	7.0 / 30.0, 
	0.0, 
	127.0 / 630.0, 
	0.0, 
	4369.0 / 22680.0, 
	0.0,
	34807.0 / 178200.0, 
	0.0, 
	20036983.0 / 97297200.0, 
	0.0, 
	2280356863.0 / 10216206000.0, 
	0.0,
	49020204823.0 / 198486288000.0, 
	0.0, 
	65967241200001.0 / 237588086736000.0, 
	0.0,
	15773461423793767.0 / 49893498214560000.0, 
	0.0, 
	655889589032992201.0 /1803293578326240000.0, 
	0.0, 
	94020690191035873697.0 / 222759794969712000000.0, 
	0.0,
	655782249799531714375489.0 / 1329207696584271504000000.0, 
	0.0,
	44737200694996264619809969.0 / 77094046401887747232000000.0, 
	0.0,
	10129509912509255673830968079.0 / 14761242414008506896480000000.0, 
	0.0,
	108026349476762041127839800617281.0 / 132496911908140357902804480000000.0, 
	0.0,
	10954814567103825758202995557819063.0 / 11262237512191930421738380800000000.0, 
	0.0,
	61154674195324330125295778531172438727.0 / 52504551281838779626144331289600000000.0, 
	0.0,
	54441029530574028687402753586278549396607.0 / 38905872499842535702972949485593600000000.0, 
	0.0,
	452015832786609665624579410056180824562551.0 / 268090886133368733415443853598208000000000.0, 
	0.0,
	2551405765475004343830620568825540664310892263.0 / 1252532276140582782027102181569679872000000000.0, 
	0.0,
	70358041406630998834159902148730577164631303295543.0 / 28520159927721069946757116674341610685440000000000.0, 
	0.0,
	775752883029173334450858052496704319194646607263417.0 / 259078091444256105986928093487086396226560000000000.0, 
	0.0,
	132034545522738294934559794712527229683368402215775110881.0 / 36256424429074976496234665114956818633529712640000000000.0
};

double gausrand(double S = 1, double U = 0)
{
	return sqrt(-2 * log((1 + rand()) / float(RAND_MAX + 2))) * cos(2 * PI * (rand() / float(RAND_MAX + 1))) * S + U;
}

struct xorshift32_state {
	uint32_t a;
};

#define xorshift(x) (x = (((x >> 13) ^ x) ^ ((x << 7) ^ (x >> 5)) ^ ((x >> 7) ^ x) ^ ((x << 11) ^ (x >> 13))))

/* The state word must be initialized to non-zero */
uint64_t xorshift64(uint64_t &x)
{
	/* Algorithm "xor" from p. 4 of Marsaglia, "Xorshift RNGs" */
	x ^= x << 13;
	x ^= x >> 7;
	x ^= x << 23;
	return x;
}

const int K = 1000;

float *M = new float[K];

const float KR16 = 1 / float(UINT16_MAX);
const float KR32 = 1 / float(UINT32_MAX);
const float KR64 = 1 / float(UINT64_MAX);

int main()
{
	float T[]
	{
		190.514496, 380.536713, 341.558716, 283.494476, 288.809174, 171.687790, 81.069801, 288.828705, 126.621162, 339.867828, 340.049225, 337.240143, 361.827393, 195.700806, 94.907021, 202.122375, 402.253906, 328.798096, 214.877075, 182.136581, 114.141502, 314.352966, 292.432373, 426.835541, 160.651031, 178.609024, 199.406479, 393.364990, 351.652527, 387.204071, 346.910675, 157.349564, 261.236145, 205.784317, 109.003815, 382.713196, 239.314224, 470.103973, 402.226532, 326.194458
	};

	float A[40];

	std::string s = "";

	for (int u = 0; u < 5; u++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				A[u * 4 * 2 + (3 - i) * 2 + (1 - j)] = T[u * 4 * 2 + i * 2 + j];
			}
		}
	}

	for (int i = 0; i < 40; i++)
	{
		s = s + std::to_string(A[i]) + ' ';
	}

	std::cout << s;






	return 0;
	uint32_t a = 0b1000101001110110100010100111011010001010011101101000101001110110;
	uint64_t b = 0b1000101001110110100010100111011010001010011101101000101001110110;
	uint64_t c = 0b1000101001110110100010100111011010001010011101101000101001110110;
	uint64_t C = 1000000000;
	float _a, _b, _c;
	for (int i = 0; i < K; i++)
	{
		M[i] = xorshift(c) * KR64;
	}
	auto t_start = std::chrono::high_resolution_clock::now();
	auto t_end = std::chrono::high_resolution_clock::now();
	t_start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < C; i++)
		_a = xorshift(a) * KR16;
	t_end = std::chrono::high_resolution_clock::now();
	float T1 = std::chrono::duration<double, std::milli>(t_end - t_start).count() / float(C);
	std::cout << "xorshift : " << T1 << " s." << std::endl;
	t_start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < C; i++)
		_b = xorshift64(b) * KR64;
	t_end = std::chrono::high_resolution_clock::now();
	float T2 = std::chrono::duration<double, std::milli>(t_end - t_start).count() / float(C);
	std::cout << "xorshift64 : " << T2 / T1 << " s." << std::endl;
	t_start = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < C; i++)
		_c = M[i % K];
	t_end = std::chrono::high_resolution_clock::now();
	float T3 = std::chrono::duration<double, std::milli>(t_end - t_start).count() / float(C);
	std::cout << "M[i % K] : " << T3 / T1 << " s." << std::endl;

	return _a + _b + _c;

	float kriteries[9];
	kriteries[0] = 20011.7;
	kriteries[1] = 20025;
	kriteries[2] = 20025.2;
	kriteries[3] = 20018.6;
	kriteries[4] = 20020;
	kriteries[5] = 20016;
	kriteries[6] = 20012.3;
	kriteries[7] = 20043.4;
	kriteries[8] = 20006.7;
	std::string ch = "012345678";


	float maxkrit = kriteries[0], minkrit = kriteries[0];

	for (int n = 1; n < ch.size(); n++)
	{
		if (maxkrit < kriteries[n])maxkrit = kriteries[n];
		if (minkrit > kriteries[n])minkrit = kriteries[n];
	}
	double summ = 0;
	for (int n = 0; n < ch.size(); n++)
	{
		kriteries[n] = (kriteries[n] - minkrit);
		kriteries[n] = kriteries[n] * kriteries[n] * kriteries[n];
		summ += kriteries[n];
	}
	for (int n = 0; n < ch.size(); n++)
		kriteries[n] /= summ;
	summ = 0;
	for (int n = 0; n < ch.size(); n++)
	{
		;
		//if (kriteries[n] < 0.1) kriteries[n] = 0;
		summ += kriteries[n];
	}
	for (int n = 0; n < ch.size(); n++)
		kriteries[n] /= summ;
	int maxn0 = 0;
	for (int n = 0; n < ch.size(); n++)
		if (kriteries[n] > kriteries[maxn0]) maxn0 = n;
	int maxn1 = 0;
	maxkrit = 1e10;
	for (int n = 0; n < ch.size(); n++)
		if (kriteries[n] > kriteries[maxn1] && kriteries[n] < kriteries[maxn0]) maxn1 = n;
	int maxn2 = 0;
	maxkrit = 1e10;
	for (int n = 0; n < ch.size(); n++)
		if (kriteries[n] > kriteries[maxn2] && kriteries[n] < kriteries[maxn1]) maxn2 = n;
	for (int n = 0; n < ch.size(); n++)
		if (kriteries[n] < kriteries[maxn2]) kriteries[n] = 0;
	kriteries[maxn0] = 10/16.0;
	kriteries[maxn1] = 4 / 16.0;
	kriteries[maxn2] = 2 / 16.0;
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
