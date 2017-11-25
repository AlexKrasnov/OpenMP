#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cmath>
#include <conio.h>
#include <omp.h>
#include <ctime>

using namespace std;

// обработка уравнения и взятие производной
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
#pragma omp parallel for shared(deltaY,K) private(i,eqlen)
	for (i = 0; i < n; i++)
	{
		eqlen = 0;
		K[i * 4] = h * Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0_iter, 0, deltaY, eqlen);
		deltaY[i] = K[i * 4] / 2;
	}
#pragma omp parallel for shared(deltaY,K) private(i,eqlen) default(shared)
	for (i = 0; i < n; i++)
	{
		eqlen = 0;
		K[i * 4 + 1] = h * Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0_iter, h / 2, deltaY, eqlen);
		deltaY[i] = K[i * 4 + 1] / 2;
	}
#pragma omp parallel for shared(deltaY,K) private(i,eqlen)
	for (i = 0; i < n; i++)
	{
		eqlen = 0;
		K[i * 4 + 2] = h * Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0_iter, h / 2, deltaY, eqlen);
		deltaY[i] = K[i * 4 + 2];
	}
#pragma omp parallel for shared(deltaY,K) private(i,eqlen)
	for (i = 0; i < n; i++)
	{
		eqlen = 0;
		K[i * 4 + 3] = h * Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0_iter, h, deltaY, eqlen);
	}
#pragma omp parallel for shared(tmp_res,K) private(i)
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
#pragma omp parallel for shared(results,tmp_res) private(i,d)
	for (i = 0; i < n; i++)
	{
		d = 0;
		double fk = Derive(tmp_res, buf + buf_scroll[i], buf_len + i, x0 - 3 * h, 0, deltaY, d); d = 0;
		double fk_1 = Derive(tmp_res + n, buf + buf_scroll[i], buf_len + i, x0 - 2 * h, 0, deltaY, d); d = 0;
		double fk_2 = Derive(tmp_res + 2 * n, buf + buf_scroll[i], buf_len + i, x0 - h, 0, deltaY, d); d = 0;
		double fk_3 = Derive(tmp_res + 3 * n, buf + buf_scroll[i], buf_len + i, x0, 0, deltaY, d); d = 0;
		results[i] = tmp_res[3 * n + 0 + i] + h * (55 * fk_3 - 59 * fk_2 + 37 * fk_1 - 9 * fk) / 24;
	}
#pragma omp parallel for shared(results,tmp_res) private(i,d)
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

int main(int argc, char * argv[])
{
	setlocale(LC_ALL, "Rus");
	omp_set_dynamic(0);

	int n(300); // количество уравнений
	double x0(1.0), xn(2.0), h(0.1); // координаты отрезка, шаг
	double t1(0), t2(0), dt1(0), dt2(0), boost(0) ; // время отсчёта, ускорение
	const int size_buf(10000000); // размер буфера для записи уравнений
	char *buffer = new char[size_buf]; // буфер для записи уравнений
	int *buffer_length = new int[n], *buffer_scroll = new int[n]; // буфер длин и буфер сдвига
	int m = int((xn - x0) / h + 1); // число точек
	double *results1 = new double[m * n]; // результат последовательного алгортима
	double *results2 = new double[m * n]; // результат параллельного алгоритма

	srand(unsigned int(time(NULL)));
	for (int i = 0; i < size_buf; i++)
		buffer[i] = 0;
	for (int i = 0; i < n; i++)
	{
		buffer_length[i] = 0;
		buffer_scroll[i] = 0;
		results1[i] = 1;
		results2[i] = 1;
	}
	int scroll = 0;

	// генерация произвольных уравнений
#pragma region generation
	for (int i = 0; i < n; i++)
	{
		buffer_scroll[i] = scroll;
		for (int j = 0; j < n; j++)
		{
			int temp_number = rand() % 99 + 1; // [1...99]
			int num_buffer_index;
			int *num_buffer = new int[2];
			for (num_buffer_index = 0; temp_number / 10; temp_number /= 10, num_buffer_index++)
				num_buffer[num_buffer_index] = temp_number % 10 + '0';
			num_buffer[num_buffer_index] = temp_number + '0';
			for (int length = scroll + num_buffer_index + 1, i = scroll; i < length; i++, scroll++)
				buffer[i] = num_buffer[length - i - 1];
			int op = rand() % 3;
			switch (op)
			{
			case 0: buffer[scroll++] = '+'; break;
			case 1: buffer[scroll++] = '-'; break;
			case 2: buffer[scroll++] = '*'; break;
			case 3: buffer[scroll++] = '/'; break;
			}
			temp_number = rand() % 6; // x_power
			for (int i = 0; i < temp_number - 1; i++)
			{
				buffer[scroll++] = 'x';
				buffer[scroll++] = '*';
			}
			buffer[scroll++] = 'x'; //y

			if (j != n - 1)
			{
				op = rand() % 3;
				switch (op)
				{
				case 0: buffer[scroll++] = '+'; break;
				case 1: buffer[scroll++] = '-'; break;
				case 2: buffer[scroll++] = '*'; break;
				case 3: buffer[scroll++] = '/'; break;
				}
			}
			delete[]num_buffer;
		}
		buffer_length[i] = scroll - buffer_scroll[i];
	}
#pragma endregion

	//////показательный пример для проверки
	//scroll = 1;
	//buffer = "x";
	//buffer_length[0] = 1;
	//buffer_scroll[0] = 0;

	cout << "\t\t***Многошаговый метод Адамса***\n\n";
	cout << "OpenMP\n";
	if (n == 1) cout << "Уравнение: \n" << buffer;
	else if (n > 1 && n <= 5)
	{
		cout << "Система уравнений: \n";
		int j = 0;
		int sum = buffer_length[0];
		for (int i = 0; i < scroll; i++)
		{
			if (i == sum)
			{
				sum += buffer_length[++j];
				cout << '\n';
			}
			cout << buffer[i];
		}
		cout << endl;
	}

	t1 = omp_get_wtime();
	omp_set_num_threads(1);
	Calculation(results2, x0, n, m, h, buffer, buffer_length, buffer_scroll);
	t2 = omp_get_wtime();
	dt2 = t2 - t1;
	cout << "\n----------------------------------------------------------\n";
	printf("\nВремя работы последовательного алгоритма:  %f\n", dt2);

	t1 = omp_get_wtime();
	omp_set_num_threads(4);
	Calculation(results1, x0, n, m, h, buffer, buffer_length, buffer_scroll);
	t2 = omp_get_wtime();
	dt1 = t2 - t1;
	cout << "\n----------------------------------------------------------\n";
	printf("\nВремя работы параллельного алгоритма:  %f\n", dt1);

	boost = dt2 / dt1;
	cout << "\n----------------------------------------------------------\n";
	printf("\nУскорение  %f\n", boost);
	cout << "\n----------------------------------------------------------\n";

	FILE *file = fopen("out.txt", "w");
	if (file != NULL)
	{
		fprintf(file, "Точки: ");
		for (int i = 0; i < m; i++)
			fprintf(file, "%f ", x0 + h * i); // вывод xi
		fprintf(file, "\nРезультат: ");
		for (int i = 0; i < n; i++)
		{
			fprintf(file, "\n%d уравнение: ", i + 1); // вывод yi
			for (int j = 0; j < m; j++)
				fprintf(file, "%f ", results1[i + j * n]);
		}
		fclose(file);
	}
	bool flag = true;
	for (int i = 0; i < n * m; i++)
	{
		if (results1[i] != results2[i])
			flag = false;
	}

	if (flag) std::cout << "Результаты совпадают\n";
	else std::cout << "Результаты не совпадают\n";

	//delete[] buffer;
	delete[] buffer_scroll;
	delete[] buffer_length;
	delete[] results1;
	delete[] results2;
	_getch();
	return 0;
}