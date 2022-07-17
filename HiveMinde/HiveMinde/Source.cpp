#include "opencv2/objdetect.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/videoio.hpp"
#include "Vector.h"

typedef wm::vector3f vector3f;

void main()
{
	static int S = 1000;

	cv::Mat3b m(S, S);

	vector3f 
		x,
		y, 
		z, 
		c;

	vector3f points[8]
	{
		{10,10,10},
		{-10,10,10},
		{10,-10,10},
		{-10,-10,10},
		{10,10,-10},
		{-10,10,-10},
		{10,-10,-10},
		{-10,-10,-10}
	};

	float
		cx = 0,
		cy = 0,
		cz = 100,
		lx = 0,
		ly = 10,
		gamma,
		focal;

	c[0] = cx;
	c[1] = cy;
	c[2] = cz;

	float
		D = sqrt(S * S + S * S),
		B_D = 43.2666,
		K = D / B_D;

	vector3f l{ lx, ly, 0 };

	z = l - c;



	for(int x = 0; x < S; x++)
		for (int y = 0; y < S; y++)
		{

		}

}