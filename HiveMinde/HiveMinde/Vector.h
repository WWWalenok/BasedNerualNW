#pragma once
#include <cmath>
#include <math.h>

namespace wm
{
	template <const unsigned char D, typename T>
	struct vector
	{
		T map[D];

		T& operator[](unsigned char i)
		{
			return map[i];
		}

		inline vector<D, T> operator-() const
		{
			wm::vector<D, T> ret;
			for (int i = 0; i < D; i++)
				ret.map[i] = -map[i];
			return ret;
		}

		inline vector<D, T> operator+() const
		{
			return *this;
		}

		inline T GetLenght()
		{
			T ret = T();
			for (int i = 0; i < D; i++)
				ret += map[i] * map[i];
			return sqrt(ret);
		}

		inline void SetLenght(T k)
		{
			T r = GetLenght();
			if (r == 0)
				return;
			for (int i = 0; i < D; i++)
				map[i] = k / r;
		}
	};



	template<unsigned char D, typename T>
	inline T operator*(vector<D, T> l, vector<D, T> r)
	{
		T ret = T();

		for (int i = 0; i < D; i++)
			ret += l.map[i] * r.map[i];
		return ret;
	}

	template<unsigned char D, typename T>
	inline wm::vector<D, T> operator+ (vector<D, T> l, vector<D, T> r)
	{
		wm::vector<D, T> ret;

		for (int i = 0; i < D; i++)
			ret.map[i] = l.map[i] + r.map[i];
		return ret;
	}

	template<unsigned char D, typename T>
	inline wm::vector<D, T> operator- (vector<D, T> l, vector<D, T> r)
	{
		wm::vector<D, T> ret;

		for (int i = 0; i < D; i++)
			ret.map[i] = l.map[i] - r.map[i];
		return ret;
	}

	template<unsigned char D, typename T>
	inline wm::vector<D, T>operator*(vector<D, T> v, T k)
	{
		wm::vector<D, T> ret;
		for (int i = 0; i < D; i++)
			ret.map[i] = v.map[i] * k;
		return ret;
	}

	template<unsigned char D, typename T>
	inline wm::vector<D, T> operator/(vector<D, T> v, T k)
	{
		wm::vector<D, T> ret;
		for (int i = 0; i < D; i++)
			ret.map[i] = v.map[i] / k;
		return ret;
	}

	template<unsigned char D, typename T>
	inline wm::vector<D, T> operator*(T k, vector<D, T> v)
	{
		wm::vector<D, T> ret;
		for (int i = 0; i < D; i++)
			ret.map[i] = v.map[i] * k;
		return ret;
	}

	template<unsigned char D, typename T>
	inline wm::vector<D, T> operator/(T k, vector<D, T> v)
	{
		wm::vector<D, T> ret;
		for (int i = 0; i < D; i++)
			ret.map[i] = v.map[i] / k;
		return ret;
	}

	typedef vector<3, float> vector3f;

	typedef vector<3, double> vector3d;

	typedef vector<3, int> vector3i;
}

