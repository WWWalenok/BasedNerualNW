#pragma once
#include<cmath>
//
//#undef MAX;
//#undef MIN;
//
//template<typename T>
//inline T MAX(T a, T b)
//{
//	return b + (a > b) * (b - a);
//}
//
//template<typename T>
//inline T MIN(T a, T b)
//{
//	return b + (a < b) * (b - a);
//}


#define CHECK_VALID(x) if(isnan(x) || isinf(x) || isinf(-x)) throw;


double gausrand(double S = 1, double U = 0)
{
	const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

	return sqrt(-2 * log((1 + rand()) / float(RAND_MAX + 1))) * cos(2 * PI * (rand() / float(RAND_MAX))) * S + U;
}

int seed = 0;
const int arS = 5059;
float *gaus_rand_arr = new float[arS];

float arrgsrnd()
{
	return gaus_rand_arr[seed = (seed + 1) % arS];
}


struct BaseGenerator
{

	int size = 100;
	float **map;
	float **buff;
	float **buff_dop;

	float dt = 0.01;

	BaseGenerator(int _size = 100)
	{
		for(int i = 0; i < arS; i++)
		{
			gaus_rand_arr[i] = gausrand();
		}


		size = _size;
		map = new float *[size];
		for(int i = 0; i < size; i++)
		{
			map[i] = new float[size];
			for(int j = 0; j < size; j++)
				map[i][j] = 0;
		}

		buff = new float *[size];
		for(int i = 0; i < size; i++)
		{
			buff[i] = new float[size];
			for(int j = 0; j < size; j++)
				buff[i][j] = 0;
		}

		buff_dop = new float *[size];
		for(int i = 0; i < size; i++)
		{
			buff_dop[i] = new float[size];
			for(int j = 0; j < size; j++)
				buff_dop[i][j] = 0;
		}
	}

	inline float GetVal(float **map, float x, float y)
	{
		x = MAX(0, MIN(size - 1, x));
		y = MAX(0, MIN(size - 1, y));

		return map[int(x)][int(y)];
	}

	inline void AddVal(float **map, float x, float y, float val)
	{

		x = MAX(0, MIN(size - 1, x));
		y = MAX(0, MIN(size - 1, y));

		map[int(x)][int(y)] += val;
		return;
	}

	inline float GetDX(float **map, float x, float y)
	{
		float ret = 0.5 * (GetVal(map, x + 1, y) - GetVal(map, x - 1, y));

		return ret;
	}

	inline float GetDY(float **map, float x, float y)
	{
		float ret = 0.5 * (GetVal(map, x, y + 1) - GetVal(map, x, y - 1));

		return ret;
	}

	void Init(size_t seed, float max_len, float min_len, int count = 50, float ofset = 0, float ampl = 1, float color = 1.75, int h_count = 9)
	{
		float *as = new float[count];
		float *bs = new float[count];
		float *fs = new float[count];
		float *k1_s = new float[count];
		float *k2_s = new float[count];
		float *hs = new float[count];
		if(max_len < min_len)
		{
			max_len = max_len + min_len;
			min_len = max_len - min_len;
			max_len = max_len - min_len;
		}
		if(min_len <= 0)
		{
			min_len = max_len * 0.01;
		}
		{
			float step = (max_len - min_len) / MAX(count / h_count, 1.0);
			srand(seed);
			float summ = 0;

			float obs = 0;

			for(int i = 0; i < count; i++)
			{
				as[i] = rand() / float(RAND_MAX) * 3.1416;
				bs[i] = rand() / float(RAND_MAX) * 31.416;
				fs[i] = (0.33 + rand() / float(RAND_MAX) * 0.66) / (max_len - step * (i / h_count));

				k1_s[i] = (cos(as[i]) * fs[i]);
				k2_s[i] = (sin(as[i]) * fs[i]);

				hs[i] = (rand() / float(RAND_MAX) * 0.2 + 0.8) * pow(fs[i], -color);
				summ += hs[i];
			}

			for(int i = 0; i < count; i++)
			{
				hs[i] = hs[i] * ampl / summ;
			}
		}

		for(int i = 0; i < size; i++)
			for(int j = 0; j < size; j++)
			{
				float &A = map[i][j] = ofset;
				for(int k = 0; k < count; k++)
					A += hs[k] * sinf(k1_s[k] * i + k2_s[k] * j + bs[k]);
			}
	}

	struct DropElement
	{
		float
			flow = 0.02,
			base_c_max = 0.5,
			c_max_v_multy = 2,
			h_multy = 50,
			flow_v_multy = -0.01,
			fric = 0.4
			;

		float
			x = 0,
			y = 0,
			vx = 0,
			vy = 0,
			c_max = 0,
			c = 0,
			m = 1,
			t = 0
			;

		void Ubdate(BaseGenerator *parrent)
		{
			float **map = parrent->map;
			float **buff = parrent->buff;
			float &dt = parrent->dt;

			if(t == 0 && parrent->GetDX(map, x, y) == 0 && parrent->GetDY(map, x, y))
			{
				m = -10000000;
			}

			t += dt;

			if(t > 50)
				m = -10000000;

			float v = sqrtf(vx * vx + vy * vy);

			m =  fmin(m, m - dt * m * (flow + v * flow_v_multy));

			if(m < 0.01)
				m = -1;


			if(m > 0)
				c_max = c_max + dt * ((base_c_max + c_max_v_multy * v) * m - c_max);
			else
				c_max = 0;
			c_max = MAX(0.0f, c_max);

			float delta_c = dt * (c_max - c) * 0.1;
			if(m <= 0)
				delta_c = -c;
			c = c + delta_c;

			vx -= dt * parrent->GetDX(map, x, y) * h_multy;
			vy -=  dt * parrent->GetDY(map, x, y) * h_multy;


			vx -=  dt * vx * fric;
			vy -=  dt * vy * fric;

			vx += arrgsrnd() * dt * 10;
			vy += arrgsrnd() * dt * 10;


			x = x + dt * vx;
			y = y + dt * vy;

			x = MAX(MIN(parrent->size - 1.0f, x), 0.0f);
			y = MAX(MIN(parrent->size - 1.0f, y), 0.0f);

			parrent->AddVal(buff, x, y, -delta_c);

			if(m > 0)if(x == 0 || y == 0 || x == parrent->size - 1.0f || y == parrent->size - 1.0f)
			{
				m = 0.0001;
			}
		}
	};

	void SymulateDrop(int count, int step, DropElement *base = 0)
	{
		const static float
			__K__ = 0.1,
			__V0__ = -1 * __K__,
			__V1__ = 1 / 4.0f * __K__,
			__V2__ = 0 * __K__;

		const static float MF[3][3]
		{
			__V2__, __V1__, __V2__,
			__V1__, __V0__, __V1__,
			__V2__, __V1__, __V2__
		};

		bool is_need_delete_base = base == 0;
		if(is_need_delete_base)
		{
			base = new DropElement;
		}

		DropElement *stack = new DropElement[count];



		for(int s = 0; s < step; s++)
		{
			for(int i = 0; i < count; i++)
			{
				stack[i] = *base;
				stack[i].x = 1 + rand() / float(RAND_MAX) * (size - 3);
				stack[i].y = 1 + rand() / float(RAND_MAX) * (size - 3);
			}

			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++)
					buff[i][j] = buff_dop[i][j] = 0;

			int c = 1;
			while(c > 0)
			{
				c = 0;
				for(int i = 0; i < count; i++) if(stack[i].m > 0)
				{
					c++;
					stack[i].Ubdate(this);
				}
			}
			float _max = 0;
			for(int k = 0; k < 20 || _max > 5; k++)
			{
				for(int i = 0; i < size; i++)
					for(int j = 0; j < size; j++)
						buff_dop[i][j] = buff[i][j];

				_max = 0;
				for(int i = 0; i < size; i++)
					for(int j = 0; j < size; j++)
					{
						float h = fabs(buff[i][j]);
						_max = fmax(h, _max);
						if(h > 0.0001 && !(h < 20 / float(k) && k >= 20))
						{

							int
								X[3]
							{
								MAX(0, i - 1),
								i,
								MIN(size - 1, i + 1)
							},
								Y[3]
							{
								MAX(0, j - 1),
								j,
								MIN(size - 1, j + 1)
							};
							buff_dop[X[1]][Y[1]] += 0.1 * (-4 * buff[X[1]][Y[1]]);
							buff_dop[X[0]][Y[1]] += 0.1 * (buff[X[1]][Y[1]]);
							buff_dop[X[2]][Y[1]] += 0.1 * (buff[X[1]][Y[1]]);
							buff_dop[X[1]][Y[0]] += 0.1 * (buff[X[1]][Y[1]]);
							buff_dop[X[1]][Y[2]] += 0.1 * (buff[X[1]][Y[1]]);

						}
					}

				auto t = buff;
				buff = buff_dop;
				buff_dop = t;
			}
			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++) if(buff[i][j] != 0)
					map[i][j] += buff[i][j] * 1;
		}

		delete[] stack;

		if(is_need_delete_base)
		{
			delete base;
		}
	}

	void Eroze(int steep)
	{
		for(int k = 0; k < steep; k++)
		{
			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++)
					buff_dop[i][j] = map[i][j];

			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++)
				{

					int
						X[3]
					{
						MAX(0, i - 1),
						i,
						MIN(size - 1, i + 1)
					},
						Y[3]
					{
						MAX(0, j - 1),
						j,
						MIN(size - 1, j + 1)
					};
					map[i][j] =
						buff_dop[i][j] + 0.1 *
						(
						-4 * buff_dop[X[1]][Y[1]]
						+ buff_dop[X[0]][Y[1]]
						+ buff_dop[X[2]][Y[1]]
						+ buff_dop[X[1]][Y[0]]
						+ buff_dop[X[1]][Y[2]]
						)
						;
				}
		}
	}

	void Expoze(int steep)
	{
		for(int k = 0; k < steep; k++)
		{
			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++)
					buff_dop[i][j] = map[i][j];

			for(int i = 0; i < size; i++)
				for(int j = 0; j < size; j++)
				{

					int
						X[3]
					{
						MAX(0, i - 1),
						i,
						MIN(size - 1, i + 1)
					},
						Y[3]
					{
						MAX(0, j - 1),
						j,
						MIN(size - 1, j + 1)
					};
					map[i][j] =
						buff_dop[i][j] + 0.1 *
						(
						4 * buff_dop[X[1]][Y[1]]
						- buff_dop[X[0]][Y[1]]
						- buff_dop[X[2]][Y[1]]
						- buff_dop[X[1]][Y[0]]
						- buff_dop[X[1]][Y[2]]
						)
						;
				}
		}
	}
};


