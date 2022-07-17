#pragma once
#include<cmath>


const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

double gausrand(double S = 1, double U = 0)
{
	const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

	return sqrt(-2 * log((1 + rand()) / float(RAND_MAX + 1))) * cos(2 * PI * (rand() / float(RAND_MAX))) * S + U;
}



template<const int map_size>
struct Generator_1
{
	float *map[map_size];

	struct DropObject
	{
		float
			x = 0,
			y = 0,
			dx = 0,
			dy = 0,
			v = 0,
			m = 0;


		DropObject()
		{
		}

	};

	inline float GetVal(float **map, float x, float y)
	{
		int
			ix = x,
			iy = y;

		ix = MAX(0, MIN(map_size - 1, ix));
		iy = MAX(0, MIN(map_size - 1, iy));

		float
			fx = x - ix,
			fy = y - iy;

		int m[2][2]
		{
			{ix, MIN(ix + 1, map_size - 1)},
			{iy, MIN(iy + 1, map_size - 1)}
		};

		return 0
			+ map[m[0][0]][m[1][0]] * (1 - fx) * (1 - fy)
			+ map[m[0][1]][m[1][0]] * (fx) * (1 - fy)
			+ map[m[0][0]][m[1][1]] * (1 - fx) * (fy)
			+ map[m[0][1]][m[1][1]] * (fx) * (fy)
			;
	}

	inline void AddVal(float **map, float x, float y, float val)
	{
		int
			ix = x,
			iy = y;

		ix = MAX(0, MIN(map_size - 1, ix));
		iy = MAX(0, MIN(map_size - 1, iy));

		float
			fx = x - ix,
			fy = y - iy;

		int m[2][2]
		{
			{ix, MIN(ix + 1, map_size - 1)},
			{iy, MIN(iy + 1, map_size - 1)}
		};

		map[m[0][0]][m[1][0]] += val * (1 - fx) * (1 - fy);
		map[m[0][1]][m[1][0]] += val * (fx) * (1 - fy);
		map[m[0][0]][m[1][1]] += val * (1 - fx) * (fy);
		map[m[0][1]][m[1][1]] += val * (fx) * (fy);

	}

	inline float GetDX(float **map, float x, float y)
	{
		return 0.5 * (GetVal(map, x + 1, y) - GetVal(map, x - 1, y));
	}

	inline float GetDY(float **map, float x, float y)
	{
		return 0.5 * (GetVal(map, x, y + 1) - GetVal(map, x, y - 1));
	}

	Generator_1()
	{
		for(int i = 0; i < map_size; i++) map[i] = new float[map_size];
	}

	template<const int count>
	void Init(size_t seed, float ofset, float ampl, float color, float max_len, float min_len, int h_count)
	{
		float as[count];
		float bs[count];
		float bsx[count];
		float bsy[count];
		float cs[count];
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
			srand(seed);
			float summ = 0;

			float obs = 0;

			for(int i = 0; i < count; i++)
			{
				float b = 1 / float(max_len - (max_len - min_len) * (i / h_count) / float(count / h_count));
				obs = b;
				b = gausrand(0.1, 1) * b;
				if(b * b < 0.01 * obs * obs)
				{
					b = 0.1 * obs;
				}
				float a = rand() / float(RAND_MAX) * PI * 2;
				as[i] = gausrand() * powf(fabs(b), -color);
				bs[i] = b;
				bsx[i] = cosf(a) * bs[i];
				bsy[i] = sinf(a) * bs[i];
				cs[i] = gausrand(2 * PI) * bs[i];

				summ += fabs(as[i]);

				obs = bs[i];
			}

			for(int i = 0; i < count; i++)
			{
				as[i] = as[i] * ampl / summ;
			}
		}

		for(int i = 0; i < map_size; i++)
			for(int j = 0; j < map_size; j++)
			{
				float &A = map[i][j] = ofset;
				for(int k = 0; k < count; k++)
					A += as[k] * sinf(bsx[k] * i + bsy[k] * j + cs[k]);
			}


	}

	template<const size_t count>
	void Drop()
	{
		DropObject *obj = new DropObject[count];

		float dt = 0.01;

		float minv = 0.00001;

		float b = 1;

		float k = 1;

		float evap = 0.5;

		float kv = 100;

		int c = 1;

		int __K = 0;

		for(int i = 0; i < count; i++)
		{
			obj[i].x = rand() / float(RAND_MAX) * (map_size - 1);
			obj[i].y = rand() / float(RAND_MAX) * (map_size - 1);

			obj[i].v = 10;
		}

		float ***buffers = new float **[2]
		{
			new float *[map_size],
			new float *[map_size]
		};

		for(int i = 0; i < map_size; i++)
		{
			buffers[0][i] = new float[map_size];
			buffers[1][i] = new float[map_size];
			for(int j = 0; j < map_size; j++)
				buffers[0][i][j] = buffers[1][i][j] = map[i][j];
		}

		while(c > 0)
		{
			c = 0;
			{
				auto _T = buffers[0];
				buffers[0] = buffers[1];
				buffers[1] = _T;
			};

			for(int i = 0; i < count; i++) if(obj[i].v > minv)
			{
				c++;
				obj[i].dx += dt * (GetDX(buffers[0], obj[i].x, obj[i].y) * 10 - b * obj[i].dx);
				obj[i].dy += dt * (GetDY(buffers[0], obj[i].x, obj[i].y) * 10 - b * obj[i].dy);
			}

			for(int i = 0; i < count; i++) if(obj[i].v > minv)
			{
				float op[2]
				{
					obj[i].x,
					obj[i].y
				};

				obj[i].x -= dt * obj[i].dx;
				obj[i].y -= dt * obj[i].dy;

				float baseCap = obj[i].v * sqrtf(obj[i].dx * obj[i].dx + obj[i].dy * obj[i].dy);

				float diff = baseCap - obj[i].m;


				obj[i].m += obj[i].v * dt * k * diff;

				AddVal(buffers[1], op[0], op[1], -obj[i].v * dt * k * diff);

				obj[i].v *= (1 - evap * dt);
			}

			for(int i = 0; i < count; i++) if(obj[i].v <= minv && obj[i].v > 0)
			{
				AddVal(buffers[1], obj[i].x, obj[i].y, obj[i].m);

				obj[i].m = 0;

				obj[i].v = 0;
			}

		}

		for(int i = 0; i < map_size; i++)
		{
			for(int j = 0; j < map_size; j++)
				map[i][j] = buffers[1][i][j];
		}


		for(int i = 0; i < map_size; i++)
		{
			delete[] buffers[0][i];
			delete[] buffers[1][i];
		}

		delete[]  buffers[0];
		delete[]  buffers[1];
		delete[]  buffers;

		delete[] obj;

	}
};

