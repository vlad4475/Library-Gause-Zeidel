// LibraryZeidel.h
#include <iostream>

#pragma once

using namespace System;

namespace LibraryZeidel {

	public ref class GaussSeidel
	{
	public:
		const int MAX_ERROR = 1000000000;
		GaussSeidel(int n);
		GaussSeidel(System::String^ file);
		~GaussSeidel();
		void PrintMatrix(System::String^ file);
		void PrintMatrix(System::String^ file, double **);
		void PrintMatrix(System::String^ file, long double **);
		int goWork(int r);
		void printAnswer(System::String^ file);
		double SolutionsError();

	private:
		long double **PreJacobi;
		long double **Matrix;
		long double *RightPath;
		long double *AnswerI;
		long double *AnswerIp1;
		long double **ImputMatrix;
		long double *InputRigthPath;
		long double **ReversToDiagonal;
		long double lastError;
		double normOfMatrixB(long double **x);
		int sizeMatrix;

		void PreconditionerSeidel();
		double esp;
	};
}
