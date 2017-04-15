// Главный DLL-файл.

#include "stdafx.h"

#include "LibraryZeidel.h"

#include <stdio.h>
#include <iostream>
#include <omp.h>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <string>

using namespace LibraryZeidel;
using namespace std;

void inversion(long double **A, int N);
void MatrixMultiplication(long double **a, long double **b, long double **c, int N);
void inversionTriangulMatrix(long double **A, int N);

GaussSeidel::GaussSeidel(System::String^ file)
{
	using namespace System::Runtime::InteropServices;
	string inputFile = (const char*)(Marshal::StringToHGlobalAnsi(file)).ToPointer();

	srand(time(NULL));
	ifstream fin(inputFile);
	char s;
	int n = 0;
	char ch;
	while (!fin.eof())
	{
		s = fin.get();
		if (s == '\n')
			n++;
	}
	fin.close();

	esp = (double)1 / n;
	PreJacobi = new long double*[n];
	ImputMatrix = new long double*[n];
	Matrix = new long double*[n];
	RightPath = new long double[n];
	InputRigthPath = new long double[n];
	AnswerI = new long double[n];
	AnswerIp1 = new long double[n];
	ReversToDiagonal = new long double*[n];
	sizeMatrix = n;
	fin.open(inputFile);
	for (int i = 0; i < n; i++)
	{
		Matrix[i] = new long double[n];
		PreJacobi[i] = new long double[n];
		ImputMatrix[i] = new long double[n];
		ReversToDiagonal[i] = new long double[n];

		for (int j = 0; j < n; j++)
		{
			fin >> Matrix[i][j];
			ImputMatrix[i][j] = Matrix[i][j];
		}
		fin >> InputRigthPath[i];
		RightPath[i] = InputRigthPath[i];
	}
}

GaussSeidel::GaussSeidel(int n)
{
	esp = 0.000001;
	lastError = 0;
	PreJacobi = new long double*[n];
	ImputMatrix = new long double*[n];
	Matrix = new long double*[n];
	RightPath = new long double[n];
	InputRigthPath = new long double[n];
	AnswerI = new long double[n];
	AnswerIp1 = new long double[n];
	ReversToDiagonal = new long double*[n];
	sizeMatrix = n;

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < n; i++)
		{
			Matrix[i] = new long double[n];
			PreJacobi[i] = new long double[n];
			ImputMatrix[i] = new long double[n];
			ReversToDiagonal[i] = new long double[n];
			InputRigthPath[i] = RightPath[i] = (rand() % 100) + 10000000; //Просто для удобства заполняем матрицу случайными значениями (ввод из файла см выше)
			AnswerI[i] = 0;
			AnswerIp1[i] = 0;
			for (int j = 0; j < n; j++)
			{
				ImputMatrix[i][j] = Matrix[i][j] = (rand() % 100); //Тоже параметры для заполнения матрицы случайными зачениями
				if (i == j)
					Matrix[i][j] = ImputMatrix[i][j] = 30000; //Тоже параметры для заполнения матрицы случайными зачениями
				else
				{
					if (rand() % 10 != 0)
					{
						Matrix[i][j] = ImputMatrix[i][j] = 0;
					}
				}
			}
		}
	}
}

GaussSeidel::~GaussSeidel()
{
	for (int i = 0; i < sizeMatrix; i++)
	{
		delete Matrix[i];
		delete PreJacobi[i];
		delete ImputMatrix[i];
		delete ReversToDiagonal[i];
	}

	delete  PreJacobi;
	delete ImputMatrix;
	delete Matrix;
	delete RightPath;
	delete InputRigthPath;
	delete AnswerI;
	delete AnswerIp1;
	delete ReversToDiagonal;
}

int GaussSeidel::goWork(int r)
{

	for (int i = 0; i < sizeMatrix; i++)
	{
		AnswerI[i] = AnswerIp1[i] = 0;
	}

	//PrintMatrix("ИсходнаяМатрица.txt");
	int countIter = 0;

	if (r == 1)
	{
	//	cout << "Вычисление предобуславливателя...." << endl;
		PreconditionerSeidel();
		//PrintMatrix("После перкосорёживания.txt");
	}
	//cout << "Норма матрицы B=  " << normOfMatrixB(Matrix) << endl;

//	cout << "Начало параллельной области" << endl;
	do
	{
		countIter++;
#pragma omp parallel shared (countIter)
		{
#pragma omp for
			for (int i = 0; i < sizeMatrix; i++)
			{
				AnswerI[i] = AnswerIp1[i];
				AnswerIp1[i] = 0;
			}
#pragma omp barrier
#pragma omp for
			for (int i = 0; i < sizeMatrix; i++)
			{
				for (int j = 0; j < sizeMatrix; j++)
				{
					if (i != j)
						AnswerIp1[i] += AnswerI[j] * ((-Matrix[i][j]) / Matrix[i][i]);
				}
				AnswerIp1[i] += RightPath[i] / Matrix[i][i];
			}
		}
	} while (SolutionsError() > esp);

	if (lastError == -1)
	{
		return EXIT_FAILURE;
	}
	lastError = 0;


	return 0;
}

double GaussSeidel::normOfMatrixB(long double **x)
{
	double si = 0;
	double smax = 0;
	for (int i = 0; i < sizeMatrix; i++)
	{
		si = 0;
		for (int j = 0; j < sizeMatrix; j++)
		{
			if (i != j)
				si += fabs(x[i][j]);
		}
		si /= x[i][i];
		if (si > smax)
			smax = si;
	}
	return smax;
}

double GaussSeidel::SolutionsError()
{
	long double valueI = 0;
	long double errorI = 0;
	long double errorMax = 0;

#pragma omp parallel private (errorI, valueI)
	{
#pragma omp for
		for (int i = 0; i < sizeMatrix; i++)
		{
			valueI = 0;
			for (int j = 0; j < sizeMatrix; j++)
			{
				valueI += ImputMatrix[i][j] * AnswerIp1[j];
			}
			errorI = fabs(valueI - InputRigthPath[i]);
#pragma omp critical (error)
			{
				if (errorI > errorMax)
					errorMax = errorI;
			}
		}
	}
	if (lastError == errorMax)
	{
		lastError = -1;
		return -1;
	}
	if (errorMax > MAX_ERROR)
	{
		//cout << "Предел погрешности превышен" << endl;
		lastError = -1;
		return -1;
	}
	lastError = errorMax;
	//cout << errorMax << endl;
	return errorMax;

}



void GaussSeidel::PrintMatrix(System::String^ file)
{
	using namespace System::Runtime::InteropServices;
	
	std::string inputFile = (const char*)(Marshal::StringToHGlobalAnsi(file)).ToPointer();

	ofstream fout;
	fout.open(inputFile);
	//fout.precision(2);
	for (int i = 0; i < sizeMatrix; i++)
	{
		for (int j = 0; j < sizeMatrix; j++)
		{
			fout << Matrix[i][j] << " ";
		}
		fout << RightPath[i] << " ; " << endl;
	}
	fout.close();
}

void GaussSeidel::PrintMatrix(System::String^ file, long double** a)
{
	using namespace System::Runtime::InteropServices;
	string inputFile = (const char*)(Marshal::StringToHGlobalAnsi(file)).ToPointer();

	ofstream fout;
	fout.open(inputFile);
	for (int i = 0; i < sizeMatrix; i++)
	{
		for (int j = 0; j < sizeMatrix; j++)
		{
			fout << a[i][j] << " ";
		}
		fout << " ;" << RightPath[i] << endl;
	}
	fout.close();
}

void GaussSeidel::PrintMatrix(System::String^ file, double** a)
{
	using namespace System::Runtime::InteropServices;
	string inputFile = (const char*)(Marshal::StringToHGlobalAnsi(file)).ToPointer();


	ofstream fout;
	fout.open(inputFile);
	for (int i = 0; i < sizeMatrix; i++)
	{
		for (int j = 0; j < sizeMatrix; j++)
		{
			fout << a[i][j] << " ";
		}
		fout << " ;" << RightPath[i] << endl;
	}
	fout.close();
}

void GaussSeidel::printAnswer(System::String^ file)
{
	using namespace System::Runtime::InteropServices;
	string inputFile = (const char*)(Marshal::StringToHGlobalAnsi(file)).ToPointer();

	ofstream foutAnswer;
	foutAnswer.open(inputFile);

	for (int j = 0; j < sizeMatrix; j++)
	{
		foutAnswer << AnswerIp1[j] << " " << endl;
	}
	foutAnswer.close();
}

void GaussSeidel::PreconditionerSeidel()
{

	long double **Mrez = new long double*[sizeMatrix];
	long double *RigthRez = new long double[sizeMatrix];

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < sizeMatrix; i++)
		{
			RigthRez[i] = 0;
			Mrez[i] = new long double[sizeMatrix];
			for (int j = 0; j < sizeMatrix; j++)
			{

				if (i <= j)
					PreJacobi[i][j] = Matrix[i][j];
				else
					PreJacobi[i][j] = 0;
			}
		}
	}

	double time = omp_get_wtime();
	//PrintMatrix("ДоИнверсииТУДАААА.txt", PreJacobi);
	inversionTriangulMatrix(PreJacobi, sizeMatrix);
	//PrintMatrix("ПослеИнверсии.txt", PreJacobi);
	cout << "Время затраченное на обращение матрици: " << omp_get_wtime() - time << endl;
	//PrintMatrix("СамПредоуслЯкоби.txt", PreJacobi);	
	//PrintMatrix("ДоПредобусл1.txt");
	time = omp_get_wtime();
	MatrixMultiplication(PreJacobi, Matrix, Mrez, sizeMatrix);
	cout << "Время затраченное на умножение матриц: " << omp_get_wtime() - time << endl;
	long double s = 0;

#pragma omp parallel
	{
#pragma omp for private(s)
		for (int i = 0; i < sizeMatrix; i++)
		{
			s = 0;
			for (int j = 0; j < sizeMatrix; j++)
			{
				s += PreJacobi[i][j] * RightPath[j];
			}
			RightPath[i] = s;
			InputRigthPath[i] = s;
		}
	}

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < sizeMatrix; i++)
		{
			for (int j = 0; j < sizeMatrix; j++)
			{
				Matrix[i][j] = Mrez[i][j];
				ImputMatrix[i][j] = Mrez[i][j];
			}
		}
	}
	//PrintMatrix("ПослеПредобусл1.txt");

	return;
}

void inversionSomeMatrix(long double **A, int N)
{
	long double temp;

	long double **E = new long double *[N];

	for (int i = 0; i < N; i++)
		E[i] = new long double[N];

	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
	{
		E[i][j] = 0.0;

		if (i == j)
			E[i][j] = 1.0;
	}

	for (int k = 0; k < N; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < N; i++)
	for (int j = 0; j < N; j++)
		A[i][j] = E[i][j];

	for (int i = 0; i < N; i++)
		delete[] E[i];

	delete[] E;
}

void inversionTriangulMatrix(long double **A, int N)
{
	long double temp;

	long double **E = new long double *[N];

	for (int i = 0; i < N; i++)
	{
		E[i] = new long double[N];
	}

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			if (i < j)
			{
				E[i][j] = 0;
				continue;
			}
			if (i == j)
			{
				E[i][j] = 1 / A[i][j];
				continue;
			}
			if (i > j)
			{
				temp = 0;
				for (int k = 0; k < j - 1; k++)
				{
					temp += A[i][k] * A[k][j];
				}
				E[i][j] = temp / (-A[j][j]);
				continue;
			}

		}
	}



#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = E[i][j];
	}

	for (int i = 0; i < N; i++)
		delete[] E[i];

	delete[] E;
}

void MatrixMultiplication(long double **a, long double **b, long double **c, int N)
{

#pragma omp parallel
	{
#pragma omp for nowait
		for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			c[i][j] = 0;
			for (int k = 0; k < N; k++)
				c[i][j] += a[i][k] * b[k][j];
		}
	}
	return;
}