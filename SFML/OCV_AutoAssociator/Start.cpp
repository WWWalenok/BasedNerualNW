#define _CRT_SECURE_NO_WARNINGS
#include "opencv2/objdetect.hpp"
#include "opencv2/highgui.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/videoio.hpp"

#include <iostream>

float norm()
{
	float a = (1 + rand()) / float(RAND_MAX);
	float b = (1 + rand()) / float(RAND_MAX);

	return sinf(6.28318530717958647 * a) * sqrtf(-2 * logf(b));
}

template<const unsigned char M, const unsigned char S>
struct Asociator
{
	float in[M * M];
	float sub[S];
	float links[M * M][S];

	Asociator()
	{
		for(int j = 0; j < S; j++)
			for(int i = 0; i < M * M; i++)
				links[i][j] = norm() * 0.1;
	}

	void set(float *I)
	{
		for(int i = 0; i < M * M; i++)
			in[i] = I[i];
	}

	void set(float **I)
	{
		for(int c = 0, i = 0; c < M; c++)
			for(int r = 0; r < M; r++, i++)
				in[i] = I[c][r];
	}

	float *get()
	{
		for(int j = 0; j < S; j++)
		{
			sub[j] = 0;
			for(int i = 0; i < M * M; i++)
				sub[j] += in[i] * links[i][j];
		}
		return sub;
	}

	float *get(float *I)
	{
		set(I);
		return get();
	}

	float *get(float **I)
	{
		set(I);
		return get();
	}
};

template<const unsigned char M, const unsigned char S>
struct Asociator_trainer
{
private:

	unsigned int __R__ = 0;

public:
	Asociator<M, S> *associator;

	float out[M * M];
	float out_err[M * M];
	float sub_err[S];
	float links_out[M * M][S];

	Asociator_trainer()
	{
		for(int j = 0; j < S; j++)
			for(int i = 0; i < M * M; i++)
				links_out[i][j] = norm() * 0.1;
	}

	void set(float *I)
	{
		associator->set(I);
	}

	void set(float **I)
	{
		associator->set(I);
	}

	float *get()
	{
		Asociator<M, S> &a = *associator;
		a.get();
		for(int i = 0; i < M * M; i++)
		{
			out[i] = 0;
			for(int j = 0; j < S; j++)
				out[i] += a.sub[j] * links_out[i][j];
		}
		return out;
	}

	float *get(float *I)
	{
		set(I);
		return get();
	}

	float *get(float **I)
	{
		set(I);
		return get();
	}

	void set_err(float *I)
	{
		for(int i = 0; i < M * M; i++)
			out_err[i] = I[i];
	}

	void set_err(float **I)
	{
		for(int c = 0, i = 0; c < M; c++)
			for(int r = 0; r < M; r++, i++)
				out_err[i] = I[c][r];
	}

	float *get_err()
	{
		for(int j = 0; j < S; j++)
		{
			sub_err[j] = 0;
			for(int i = 0; i < M * M; i++)
				sub_err[j] += out_err[i] * links_out[i][j];
		}
		return sub_err;
	}

	float *get_err(float *I)
	{
		set_err(I);
		return get_err();
	}

	float *get_err(float **I)
	{
		set_err(I);
		return get_err();
	}
	void update(float koof = 0.01)
	{

		const auto tr = RAND_MAX / 5 * 2;

		for(int j = 0; j < S; j++, __R__++)
		{
			for(int i = 0; i < M * M; i++)
				links_out[i][j] -= out_err[i] * koof * (rand() > tr);
		}

		Asociator<M, S> &a = *associator;
		for(int j = 0; j < S; j++)
		{
			for(int i = 0; i < M * M; i++, __R__++)
				a.links[i][j] -= sub_err[j] * koof * (rand() > tr);
		}
	}


};



int main(int argc, char *argv[])
{
	const unsigned char M = 7, S = 5;
	cv::Mat1b imj = cv::imread("test.jpg", cv::IMREAD_GRAYSCALE);

	cv::Mat1b imj2 = imj.clone();

	cv::Mat1f imj_mask = imj.clone();

	cv::Mat1f imj_temp[S + 1];

	for(int i = 0; i < S + 1; i++)
		imj_temp[i] = imj.clone();

	Asociator<M, S> ass;

	Asociator_trainer<M, S> trainer;

	trainer.associator = &ass;

	float mask[M * M];

	float err[M * M];


	for(int __I = 0; __I < 10000; __I++)
	{
		//for(int c = 0; c < imj.cols - M; c++)
		//	for(int r = 0; r < imj.rows - M; r++)
		for(int __J = 0; __J < 2000; __J++)
		{
			int c = rand() % (imj.cols - M);
			int r = rand() % (imj.rows - M);
			float med = 0;
			for(int a = 0, i = 0; a < M; a++)
				for(int b = 0; b < M; b++, i++)
					med += mask[i] = imj[r + a][c + b];

			med /= float(M * M);

			float disp = 0;

			for(int a = 0, i = 0; a < M; a++)
				for(int b = 0; b < M; b++, i++)
					disp += (mask[i] - med) * (mask[i] - med);

			disp = sqrtf(disp);

			if(disp < 100)
			{
				__J--;
				continue;
			}

			trainer.get(mask);

			for(int a = 0, i = 0; a < M; a++)
				for(int b = 0; b < M; b++, i++)
					err[i] = 10 * (trainer.out[i] - (imj[r + a][c + b] - med) * 1 + med) / (10 + (a - M / 2) * (a - M / 2) + (b - M / 2) * (b - M / 2));


			trainer.get_err(err);


			trainer.update(0.00000025);

		}

		for(int c = 0; c < imj.cols; c++)
			for(int r = 0; r < imj.rows; r++)
			{
				imj2[r][c] = imj_mask[r][c] = 0;
				for(int i = 0; i < S + 1; i++)
					imj_temp[i][r][c] = 0;
			}


		for(int c = 0; c < imj.cols - M; c++)
			for(int r = 0; r < imj.rows - M; r++)
			{
				for(int a = 0, i = 0; a < M; a++)
					for(int b = 0; b < M; b++, i++)
						mask[i] = imj[r + a][c + b];

				trainer.get(mask);

				{
					for(int i = 0; i < S; i++)
						imj_temp[i][r][c] = tanhf(trainer.associator->sub[i] * 0.004);				
				}

				for(int a = 0, i = 0; a < M; a++)
					for(int b = 0; b < M; b++, i++) if((a - M / 2) * (a - M / 2) + (b - M / 2) * (b - M / 2) < 1)
					{
						imj_temp[S][r + a][c + b] += abs(trainer.out[i]);
						imj_mask[r + a][c + b] += 1;
					}
			}

		for(int c = 0; c < imj.cols; c++)
			for(int r = 0; r < imj.rows; r++)
				imj2[r][c] = imj_temp[S][r][c] / imj_mask[r][c];
		cv::imshow("base", imj);


		//for(int i = 0; i < S; i++)
		//	cv::imshow(std::to_string(i), imj_temp[i]);


		cv::imshow("new", imj2);

		cv::waitKey(1);
	}

	return 0;
};
