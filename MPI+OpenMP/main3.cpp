#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <omp.h>
#include <cstring>

using namespace std;

double Derive(double *tmp_res, char *buf, int *buf_len, double x0, double deltaX, double *deltaY, int &i)
{
	double r = 0; int op = 0; char ch = ' ';
	while (i < buf_len[0])
	{
		ch = buf[i];
		if ((ch >= '0') && (ch <= '9'))
		{
			int temp = 0, j = 0;
			while ((buf[i] >= '0') && (buf[i] <= '9') && (i < buf_len[0]))
			{
				temp *= int(pow(double(10), j)) + buf[i] - '0'; j++; i++;
			}
			if (r == 0) r = temp;
			else
				switch (op)
			{
				case '*':
					if (buf[i] == '^') { i = i - j; r *= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i); }
					else r *= temp;
					break;
				case '/':
					if (buf[i] == '^') { i = i - j; r /= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i); }
					else r /= temp;
					break;
				case '^': r = pow(r, temp);
				case '+':
					if ((buf[i] == '*') || (buf[i] == '/') || (buf[i] == '^'))
					{
						i = i - j;
						r += Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);
					}
					else r += temp;
					break;
				case '-':
					if ((buf[i] == '*') || (buf[i] == '/') || (buf[i] == '^'))
					{
						i = i - j;
						r -= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);
					}
					else r -= temp;
					break;
			}
		}
		else if (ch == 'x')
			if (r == 0) { i++; r = x0 + deltaX; }
			else
				switch (op)
			{
				case '*':
					if (buf[i + 1] == '^') r *= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);
					else { i++; r *= (x0 + deltaX); }
					break;
				case '/':
					if (buf[i + 1] == '^') r /= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);
					else { i++; r /= (x0 + deltaX); }
					break;
				case '+':
					if ((buf[i + 1] == '*') || (buf[i + 1] == '/') || (buf[i + 1] == '^'))
						r += Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);
					else { i++; r += (x0 + deltaX); }
					break;
				case '-':
					if ((buf[i + 1] == '*') || (buf[i + 1] == '/') || (buf[i + 1] == '^'))
						r -= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);
					else { i++; r -= (x0 + deltaX); }
					break;
			}
		else if (ch == '(')
		{
			int li = i; i++;
			double temp = Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);
			if (r == 0) r = temp;
			else
				switch (op)
			{
				case '*':
					if (buf[i] == '^') { i = li; r *= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i); }
					else r *= temp;
					break;
				case '/':
					if (buf[i] == '^') { i = li; r *= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i); }
					else r /= temp;
					break;
				case '+':
					if ((buf[i] == '*') || (buf[i] == '/') || (buf[i] == '^'))
					{
						i = li;
						r += Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);
					}
					else r += temp;
					break;
				case '-':
					if ((buf[i] == '*') || (buf[i] == '/') || (buf[i] == '^'))
					{
						i = li;
						r -= Derive(tmp_res, buf, buf_len, x0, deltaX, deltaY, i);;
					}
					else r -= temp;
					break;
			}
		}
		else if (strchr(")+-/*^", ch) != NULL)
		{
			i++;
			if (ch == ')') return r;
			else if (ch == '*') op = '*';
			else if (ch == '/') op = '/';
			else if (ch == '+') op = '+';
			else if (ch == '-') op = '-';
			else if (ch == '^') op = '^';
		}
	}
	return r;
}

void RK(double *tmp_res, double x0_iter, int n, double h, char *buf, int * buf_len, int * buf_scroll)
{
	int i(0), eqlen(0);
	double * K = new double[4 * n];
	double * deltaY = new double[n];
	for (int i = 0; i < n; i++)
		deltaY[i] = 0;
	for (i = 0; i < n; i++)
	{
		eqlen = 0;
		K[i * 4] = h * Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0_iter, 0, deltaY, eqlen);
		deltaY[i] = K[i * 4] / 2;
	}
	for (i = 0; i < n; i++)
	{
		eqlen = 0;
		K[i * 4 + 1] = h * Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0_iter, h / 2, deltaY, eqlen);
		deltaY[i] = K[i * 4 + 1] / 2;
	}
	for (i = 0; i < n; i++)
	{
		eqlen = 0;
		K[i * 4 + 2] = h * Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0_iter, h / 2, deltaY, eqlen);
		deltaY[i] = K[i * 4 + 2];
	}
	for (i = 0; i < n; i++)
	{
		eqlen = 0;
		K[i * 4 + 3] = h * Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0_iter, h, deltaY, eqlen);
	}
	for (i = 0; i < n; i++)
	{
		double deltay0 = (K[i * 4] + 2 * K[i * 4 + 1] + 2 * K[i * 4 + 2] + K[i * 4 + 3]) / 6;
		tmp_res[n + i] = tmp_res[i] + deltay0;
	}
	delete[]deltaY;
	delete[]K;
}

void Adams(double *tmp_res, double *results, double x0, int n, double h, char *buf, int * buf_len, int * buf_scroll)
{
	double * deltaY = new double[n];
	int i(0), d(0);
	for (i = 0; i < n; i++)
		deltaY[i] = 0;
	for (i = 0; i < n; i++)
	{
		d = 0;
		double fk = Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0 - 3 * h, 0, deltaY, d); d = 0;
		double fk_1 = Derive(tmp_res + n, buf + buf_scroll[i], buf_len + i, x0 - 2 * h, 0, deltaY, d); d = 0;
		double fk_2 = Derive(tmp_res + 2 * n, buf + buf_scroll[i], buf_len + i, x0 - h, 0, deltaY, d); d = 0;
		double fk_3 = Derive(tmp_res + 3 * n, buf + buf_scroll[i], buf_len + i, x0, 0, deltaY, d); d = 0;
		results[i] = tmp_res[3 * n + 0 + i] + h * (55 * fk_3 - 59 * fk_2 + 37 * fk_1 - 9 * fk) / 24;
	}
	for (i = 0; i < n; i++)
	{
		d = 0;
		double fk_1 = Derive(tmp_res + n, buf + buf_scroll[i], buf_len + i, x0 - 2 * h, 0, deltaY, d); d = 0;
		double fk_2 = Derive(tmp_res + 2 * n, buf + buf_scroll[i], buf_len + i, x0 - h, 0, deltaY, d); d = 0;
		double fk_3 = Derive(tmp_res + 3 * n, buf + buf_scroll[i], buf_len + i, x0, 0, deltaY, d); d = 0;
		double fk_4 = Derive(results, buf + buf_scroll[i], buf_len + i, x0 + h, 0, deltaY, d); d = 0;
		results[i] = tmp_res[3 * n + 0 + i] + h * (9 * fk_4 + 19 * fk_3 - 5 * fk_2 + fk_1) / 24;
	}
	delete[]deltaY;
}

void Calculation(double *results, double x0, int n, int m, double h, char *buf, int * buf_len, int * buf_scroll)
{
	int tmp_size = m * n;
	double *tmp_res = new double[tmp_size]; // временные результаты
	for (int i = 0; i < n; i++)
		tmp_res[i] = results[i];
	double x0_iter = x0;
	// вычисляем значения в 3 точках методом Ранге-Кутты
	for (int i = 0; i < 3; i++)
	{
		RK(tmp_res + i * n, x0, n, h, buf, buf_len, buf_scroll);
		for (int j = 0; j < n; j++)
			results[(i + 1) * n + j] = tmp_res[(i + 1) * n + j];
		x0 += h;
	}
	// выполняем метод Адамса в оставшихся точках отрезка
	for (int i = 4; i < m; i++)
	{
		Adams(tmp_res, results + i * n, x0, n, h, buf, buf_len, buf_scroll);
		for (int j = 0; j < n; j++)
			for (int k = 0; k < 4; k++)
				if (k != 3) tmp_res[k * n + j] = tmp_res[(k + 1) * n + j];
				else tmp_res[k * n + j] = results[i * n + j];
				x0 += h;
	}
	delete[] tmp_res;
}

void RK(double * vs, double x0, int size, double step, char * local_strbuf, int * local_strbuf_len, int * local_strbuf_off, int equation_num, int procRank, int * equation_count, int * ec_offset)
{
	int i(0);
	double * K = new double[4 * equation_num];
	double * deltaY = new double[size];
	for (int i = 0; i < size; i++)
		deltaY[i] = 0;
#pragma omp parallel for shared(deltaY, K) private(i) default(shared)
	for (i = 0; i < equation_num; i++)
	{
		int eqlen(0);
		K[i * 4] = step * Derive(vs, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0, 0, deltaY, eqlen);
		deltaY[ec_offset[procRank] + i] = K[i * 4] / 2;
	}
	MPI_Allgatherv(deltaY + ec_offset[procRank], equation_count[procRank], MPI_DOUBLE,
		deltaY, equation_count, ec_offset, MPI_DOUBLE, MPI_COMM_WORLD);
#pragma omp parallel for shared(deltaY, K) private(i) default(shared)
	for (i = 0; i < equation_num; i++)
	{
		int eqlen(0);
		K[i * 4 + 1] = step * Derive(vs, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0, step / 2, deltaY, eqlen);
		deltaY[ec_offset[procRank] + i] = K[i * 4 + 1] / 2;
	}
	MPI_Allgatherv(deltaY + ec_offset[procRank], equation_count[procRank], MPI_DOUBLE,
		deltaY, equation_count, ec_offset, MPI_DOUBLE, MPI_COMM_WORLD);
#pragma omp parallel for shared(deltaY, K) private(i) default(shared)
	for (i = 0; i < equation_num; i++)
	{
		int eqlen(0);
		K[i * 4 + 2] = step * Derive(vs, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0, step / 2, deltaY, eqlen);
		deltaY[ec_offset[procRank] + i] = K[i * 4 + 2];
	}
	MPI_Allgatherv(deltaY + ec_offset[procRank], equation_count[procRank], MPI_DOUBLE,
		deltaY, equation_count, ec_offset, MPI_DOUBLE, MPI_COMM_WORLD);
#pragma omp parallel for shared(deltaY, K) private(i) default(shared)
	for (i = 0; i < equation_num; i++)
	{
		int eqlen(0); 
		K[i * 4 + 3] = step * Derive(vs, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0, step, deltaY, eqlen);
	}
#pragma omp parallel for shared(vs, K) private(i) default(shared)
	for (i = 0; i < equation_num; i++)
	{
		double deltay0 = (K[i * 4] + 2 * K[i * 4 + 1] + 2 * K[i * 4 + 2] + K[i * 4 + 3]) / 6;
		vs[size + ec_offset[procRank] + i] = vs[ec_offset[procRank] + i] + deltay0;
	}
	delete[]deltaY;
	delete[]K;
}

void Adams(double * vs, double * results, double x0, int size, double step, char * local_strbuf, int * local_strbuf_len, int * local_strbuf_off, int equation_num, int procRank, int * equation_count, int * ec_offset)
{
	int i(0),d(0);
	double * deltaY = new double[size];
	for (int i = 0; i < size; i++)
		deltaY[i] = 0;
#pragma omp parallel for shared(results,vs) private(i,d)
	for (i = 0; i < equation_num; i++)
	{
		d = 0;
		double fk = Derive(vs, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0 - 3 * step, 0, deltaY, d); d = 0;
		double fk_1 = Derive(vs + size, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0 - 2 * step, 0, deltaY, d); d = 0;
		double fk_2 = Derive(vs + 2 * size, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0 - step, 0, deltaY, d); d = 0;
		double fk_3 = Derive(vs + 3 * size, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0, 0, deltaY, d); d = 0;
		results[ec_offset[procRank] + i] = vs[3 * size + ec_offset[procRank] + i] + step * (55 * fk_3 - 59 * fk_2 + 37 * fk_1 - 9 * fk) / 24;
		d = 0;
	}
	MPI_Allgatherv(results + ec_offset[procRank], equation_count[procRank], MPI_DOUBLE,
		results, equation_count, ec_offset, MPI_DOUBLE, MPI_COMM_WORLD);
#pragma omp parallel for shared(results,vs) private(i,d)
	for (i = 0; i < equation_num; i++)
	{
		d = 0;
		double fk_1 = Derive(vs + size, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0 - 2 * step, 0, deltaY, d); d = 0;
		double fk_2 = Derive(vs + 2 * size, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0 - step, 0, deltaY, d); d = 0;
		double fk_3 = Derive(vs + 3 * size, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0, 0, deltaY, d); d = 0;
		double fk_4 = Derive(results, local_strbuf + local_strbuf_off[i], local_strbuf_len + i, x0 + step, 0, deltaY, d); d = 0;
		results[ec_offset[procRank] + i] = vs[3 * size + ec_offset[procRank] + i] + step * (9 * fk_4 + 19 * fk_3 - 5 * fk_2 + fk_1) / 24;
		d = 0;
	}
	delete[]deltaY;
}

int main(int argc, char * argv[])
{
	int procRank, procSize;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	const int very_big_number(10000000);
	omp_set_dynamic(0);
	omp_set_num_threads(2);
	int size(1);//number of equations
	if (argc > 1)
		size = atoi(argv[1]);
	char *string_buffer = 0, // global
		*local_strbuf; // local string of equation
	int *local_strbuf_len, // local string of equations length
		*length_of_equation, //length of each equation in string_buffer
		*local_strbuf_off, // local string of equations offset
		*offset_of_equation; //start of each equation in string_buffer
	double *results, *results2, *params, //buffer for size, x0, xn, step;
		x0(0), xn(0), step(0),
		time_start(0), time_end(0), paralleltime(0), serialtime(0);
	if (procRank == 0)
	{
		srand((unsigned)time(NULL));
		x0 = 1.0;
		xn = 2.0;
		step = 0.1;
		if (argc > 4)
		{
			x0 = atof(argv[2]);
			xn = atof(argv[3]);
			step = atof(argv[4]);
		}
		string_buffer = new char[very_big_number];
		for (int i = 0; i < very_big_number; i++)
			string_buffer[i] = NULL;
		length_of_equation = new int[size];
		offset_of_equation = new int[size];
		results = new double[int((xn - x0) / step + 1) * size];
		results2 = new double[int((xn - x0) / step + 1) * size];
		for (int i = 0; i < size; i++)
		{
			length_of_equation[i] = 0;
			offset_of_equation[i] = 0;
			results[i] = 1;
			results2[i] = 1;
		}

		int string_offset = 0; //global offset
		for (int i = 0; i < size; i++)
		{ //generation of equations, each of them contains size elements
			offset_of_equation[i] = string_offset;
			for (int j = 0; j < size; j++)
			{
				int temp_number = rand() % 99 + 1; // [1...99]
				int num_buffer_index = 0;
				int *num_buffer = new int[2];
				for (; temp_number / 10; temp_number /= 10, num_buffer_index++)
				{ //132
					num_buffer[num_buffer_index] = temp_number % 10 + '0'; // Number to ASCII
				}
				num_buffer[num_buffer_index] = temp_number + '0'; //num_buffer = [2 3 1]
				for (int length = string_offset + num_buffer_index + 1, i = string_offset; i < length; i++, string_offset++)
					string_buffer[i] = num_buffer[length - i - 1];
				int sign = rand() % 3;
				switch (sign){
				case 0:
					string_buffer[string_offset++] = '+';
					break;
				case 1:
					string_buffer[string_offset++] = '-';
					break;
				case 2:
					string_buffer[string_offset++] = '*';
					break;
				case 3:
					string_buffer[string_offset++] = '/';
					break;
				}
				temp_number = rand() % 6; // x_power
				for (int i = 0; i < temp_number - 1; i++){
					string_buffer[string_offset++] = 'x';
					string_buffer[string_offset++] = '*';
				}
				string_buffer[string_offset++] = 'x';

				if (j != size - 1)
				{
					sign = rand() % 3;
					switch (sign){
					case 0:
						string_buffer[string_offset++] = '+';
						break;
					case 1:
						string_buffer[string_offset++] = '-';
						break;
					case 2:
						string_buffer[string_offset++] = '*';
						break;
					case 3:
						string_buffer[string_offset++] = '/';
						break;
					}
				}
				delete[]num_buffer;
			}
			length_of_equation[i] = string_offset - offset_of_equation[i];
		}

		//////показательный пример для проверки
		//string_offset = 1;
		//string_buffer = "x";
		//length_of_equation[0] = 1;
		//offset_of_equation[0] = 0;


		time_start = MPI_Wtime();
		Calculation(results2, x0, size, int((xn - x0) / step + 1), step, string_buffer, length_of_equation, offset_of_equation);
		time_end = MPI_Wtime();
		serialtime = time_end - time_start;
		printf("SerialTime = %f\n", serialtime);

		params = new double[4];
		params[0] = (double)size;
		params[1] = x0;
		params[2] = xn;
		params[3] = step;
		time_start = MPI_Wtime();
	}
	else params = new double[4];
	MPI_Bcast(params, 4, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (procRank != 0)
	{
		size = (int)params[0];
		x0 = params[1];
		xn = params[2];
		step = params[3];
		results = new double[int((xn - x0) / step + 1) * size];
		length_of_equation = new int[size];
		offset_of_equation = new int[size];
	}
	delete[]params;
	MPI_Bcast(results, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int *equation_count = new int[procSize];
	int *ec_offset = new int[procSize];
	int *sb_offset = new int[procSize];
	int *sb_count = new int[procSize];

	MPI_Bcast(length_of_equation, size, MPI_INT, 0, MPI_COMM_WORLD);

	int temp_offset = 0, temp_count = 0;
	for (int i = 0; i < procSize; i++)
	{ 
		equation_count[i] = size / procSize;
		if (size % procSize > i) 
			equation_count[i]++;
		sb_offset[i] = temp_offset;
		ec_offset[i] = temp_count;
		int dtemp_count = temp_count + equation_count[i];
		for (temp_count; temp_count < dtemp_count; temp_count++)
			temp_offset += length_of_equation[temp_count];
		if (procRank == i)
			local_strbuf = new char[temp_offset - sb_offset[i]];
		sb_count[i] = temp_offset - sb_offset[i]; // 5
	}
	if ((procRank == 0) && (procSize > 1))
	{
		//for (size_t i = 0; i < size; i++)
		//{
		//    printf("off=%d\n", offset_of_equation[i]);
		//}
		//printf("\n");
		int temp = offset_of_equation[ec_offset[1]];
		for (int j = 0, i = equation_count[j]; i < size; i++)
		{
			//printf("i = %d sum=%d\n", i, equation_count[j + 1] + ec_offset[j + 1]);
			if (i >= (equation_count[j + 1] + ec_offset[j + 1])){
				j++;
				temp = offset_of_equation[ec_offset[j + 1]];
			}
			//printf("temp=%d\n", temp);
			offset_of_equation[i] -= temp;
		}
		//for (size_t i = 0; i < size; i++)
		//{
		//    printf("off=%d\n", offset_of_equation[i]);
		//}
	}
	local_strbuf_len = new int[equation_count[procRank]];
	local_strbuf_off = new int[equation_count[procRank]];
	MPI_Scatterv(length_of_equation, equation_count, ec_offset, MPI_INT,
		local_strbuf_len, equation_count[procRank], MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(string_buffer, sb_count, sb_offset, MPI_CHAR,
		local_strbuf, sb_count[procRank], MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Scatterv(offset_of_equation, equation_count, ec_offset, MPI_INT,
		local_strbuf_off, equation_count[procRank], MPI_INT, 0, MPI_COMM_WORLD);

	int _sz = (int((xn - x0) / step) + 1)*size;
	double *vs = new double[_sz];
	for (int i = 0; i < size; i++)
		vs[i] = results[i];
	double x0_iter = x0;

	MPI_Barrier(MPI_COMM_WORLD);

	for (int i = 0; i < 3; i++)
	{
		RK(vs + i * size, x0_iter, size, step, local_strbuf, local_strbuf_len, local_strbuf_off, equation_count[procRank], procRank, equation_count, ec_offset);
		for (int j = 0; j < equation_count[procRank]; j++)
			results[(i + 1) * size + ec_offset[procRank] + j] = vs[(i + 1) * size + ec_offset[procRank] + j];
		x0_iter += step;
		MPI_Allgatherv(results + (i + 1) * size + ec_offset[procRank], equation_count[procRank], MPI_DOUBLE,
			results + (i + 1) * size, equation_count, ec_offset, MPI_DOUBLE, MPI_COMM_WORLD);
		for (int j = 0; j < size; j++)
			vs[(i + 1) * size + j] = results[(i + 1) * size + j];
	}
	for (int i = 4; i < (int((xn - x0) / step) + 1); i++)
	{
		Adams(vs, results + i * size, x0_iter, size, step, local_strbuf, local_strbuf_len, local_strbuf_off, equation_count[procRank], procRank, equation_count, ec_offset);
		MPI_Allgatherv(results + i * size + ec_offset[procRank], equation_count[procRank], MPI_DOUBLE,
			results + i * size, equation_count, ec_offset, MPI_DOUBLE, MPI_COMM_WORLD);
		for (int j = 0; j < size; j++)
			for (int k = 0; k < 4; k++)
				if (k != 3) vs[k * size + j] = vs[(k + 1) * size + j];
				else vs[k * size + j] = results[i * size + j];
				x0_iter += step;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (procRank == 0)
	{
		time_end = MPI_Wtime();
		paralleltime = time_end - time_start;

		printf("ParallelTime = %f\n", paralleltime);

		FILE *file = fopen("out.txt", "w");
		if (file != NULL)
		{
			fprintf(file, "Точки: ");
			for (int i = 0; i < (int((xn - x0) / step) + 1); i++)
				fprintf(file, "%f ", x0 + step * i); // вывод xi
			fprintf(file, "\nРезультат: ");
			for (int i = 0; i < size; i++)
			{
				fprintf(file, "\n%d уравнение: ", i + 1); // вывод yi
				for (int j = 0; j < (int((xn - x0) / step) + 1) * size; j++)
					fprintf(file, "%f ", results[i + j * size]);
			}
			fclose(file);
		}

		printf("Boost = %f\n", serialtime / paralleltime);

		bool flag = true;
		for (int i = 0; i < size; i++)
		{
			if (results[i] != results2[i])
				flag = false;
		}
		if (flag) printf ("The results are equal\n");
		else printf ("The results are not equal\n");

		delete[] results2;
	}
	delete[]equation_count;
	delete[]ec_offset;
	delete[]length_of_equation;
	delete[]offset_of_equation;
	delete[]sb_count;
	delete[]sb_offset;
	delete[]local_strbuf;
	delete[]local_strbuf_len;
	delete[]results;
	MPI_Finalize();
	return 0;
}
