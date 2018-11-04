#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <iomanip>
#include <cmath>

template <typename T>
class Matrix {
	T **values; // значения матрицы
	int rows; // число строк
	int cols; // число столбцов

	void AllocateMemory();

public:
	Matrix(int n); // конструктор квадратной матрицы
	Matrix(int rows, int cols); // конструктор произвольной матрицы

	Matrix(const Matrix &matrix); // конструктуор копирования

	Matrix& operator=(const Matrix &matrix); // оператор присваивания

	int GetRows() const; // получение числа строк
	int GetCols() const; // получение числа столбцов

	T Det() const; // определитель
	T Track() const; // след матрицы
	int Rank() const; // ранг

	double& operator()(int i, int j); // получение элемента с возможностью изменения
	double operator()(int i, int j) const; // получение элемента без возможности изменения

	Matrix operator+(const Matrix &matrix) const; // сложение с матрицей
	Matrix operator-(const Matrix &matrix) const; // вычитание матрицы
	Matrix operator*(const Matrix &matrix) const; // умножение на матрицу

	Matrix operator*(const T& value) const; // умножение на скаляр

	Matrix Transpose() const; // получение транспонированной матрицы
	Matrix Inverse() const; // получение обратной матрицы
	Matrix Ortonormalize() const; // получение ортонормированной матрицы
	Matrix Ortogonalize() const; // получение ортогональной матрицы

	template <typename T1>
	friend std::ostream& operator<<(std::ostream &os, const Matrix<T1> matrix); // оператор вывода в поток вывода
};

template <typename T>
void Matrix<T>::AllocateMemory() {
	values = new T*[rows]; // выделяем память под строки

	for (int i = 0; i < rows; i++) {
		values[i] = new T[cols]; // выделяем память под столбцы

		for (int j = 0; j < cols; j++)
			values[i][j] = 0;
	}
}

template <typename T>
Matrix<T>::Matrix(int n) {
	if (n < 1)
		throw "Matrix::Matrix(n) - incorrect size of matrix";

	rows = n;
	cols = n;

	AllocateMemory(); // выделяем память под значения
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols) {
	if (rows < 1 || cols < 1)
		throw "Matrix::Matrix(rows, cols) - incorrect sizes of matrix";

	this->rows = rows;
	this->cols = cols;

	AllocateMemory(); // выделяем память под значения
}

template <typename T>
Matrix<T>::Matrix(const Matrix &matrix) {
	rows = matrix.rows;
	cols = matrix.cols;

	AllocateMemory(); // выделяем память под значения

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			values[i][j] = matrix.values[i][j]; // копируем значения
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix &matrix) {
	if (this == &matrix) // если присваивается сама себе
		return *this; 

	// очищаем старую память
	for (int i = 0; i < rows; i++)
		delete[] values[i];

	delete[] values;

	// копируем размеры
	rows = matrix.rows;
	cols = matrix.cols;

	AllocateMemory(); // выделяем память

	// копируем элементы
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			values[i][j] = matrix.values[i][j];

	return *this;
}

template <typename T>
int Matrix<T>::GetRows() const {
	return rows;
}

template <typename T>
int Matrix<T>::GetCols() const {
	return cols;
}

// определитель матрицы
template <typename T>
T Matrix<T>::Det() const {
	if (rows != cols) // если матрица не квадратная, то бросаем исключение
		throw "Matrix::Det() - Determinant of non squarable matrix does not exist";

	// если матрица 1 на 1
	if (rows == 1)
		return values[0][0]; // возвращаем элемент

	Matrix copy(*this); // копируем матрицу, чтобы не изменять исходную
	int sign = 1; // знак определителя
	T det = 1; // определитель

	// проходимся методом Гаусса по матрице
	for (int j = 0; j < cols; j++) {
		if (copy.values[j][j] == 0) { // если на главной диагонали 0, то ищем строку с не нулём в этом столбце
			int k = j;

			while (k < rows && !copy.values[k][j])
				k++;

			if (k == rows)
				return 0; // если не нашли такую строку, значит определитель равен нулю

			// иначе меняем строки
			for (int l = 0; l < cols; l++) {
				T tmp = copy.values[j][l];
				copy.values[j][l] = copy.values[k][l];
				copy.values[k][l] = tmp;
			}

			sign = -sign; // меняем знак
		}

		det *= copy.values[j][j]; // умножаем определитель на диагональный элемент

		for (int i = j + 1; i < rows; i++) {
			T v = -copy.values[i][j] / copy.values[j][j];

			for (int j1 = 0; j1 < cols; j1++)
				copy.values[i][j1] += copy.values[j][j1] * v;
		}
	}

	return det * sign; // возвращаем определитель с учётом знака
}

// след матрицы
template <typename T>
T Matrix<T>::Track() const {
	if (rows != cols) // если матрица не квадратная, то бросаем исключение
		throw "Matrix::Track() - track of non squarable matrix does not exist";

	T tr = 0;

	for (int i = 0; i < rows; i++)
		tr = tr + values[i][i];

	return tr;
}

// ранг матрицы
template <typename T>
int Matrix<T>::Rank() const {
	int n = rows > cols ? rows : cols;

	Matrix copy(*this); // копируем матрицу, чтобы не изменять исходную

	// приводим матрицу к треугольному виду
	for (int j = 0; j < n; j++) {
		int rowNum = j;

		while (rowNum < rows && copy.values[rowNum][j] == 0)
			rowNum++;

		if (rowNum != rows) {
			for (int k = 0; k < cols; k++) {
				T tmp = copy.values[j][k];
				copy.values[j][k] = copy.values[rowNum][k];
				copy.values[rowNum][k] = tmp;
			}

			for (int i = j + 1; i < rows; i++) {
				double elem = copy.values[i][j] / copy.values[j][j];

				for (int k = 0; k < cols; k++)
					copy.values[i][k] -= copy.values[j][k] * elem;
			}
		}
	}

	int rank = 0;

	// считаем количество ненулевых строк
	for (int i = 0; i < rows; i++) {
		int j = 0;

		while (j < cols && copy.values[i][j] == 0)
			j++;

		if (j < cols)
			rank++;
	}

	return rank; // возвращаем ранг
}

template <typename T>
double& Matrix<T>::operator()(int i, int j) {
	if (i < 0 || i >= rows || j < 0 || j >= cols)
		throw "Matrix::operator()(i, j) - index is out of bounds";

	return values[i][j];
}

template <typename T>
double Matrix<T>::operator()(int i, int j) const {
	if (i < 0 || i >= rows || j < 0 || j >= cols)
		throw "Matrix::operator()(i, j) - index is out of bounds";
	
	return values[i][j];
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix &matrix) const {
	if (rows != matrix.rows || cols != matrix.cols)
		throw "Matrix::operator+(const Matrix &) - matrixes have different sizes";

	Matrix result(rows, cols);

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result.values[i][j] = values[i][j] + matrix.values[i][j];

	return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix &matrix) const {
	if (rows != matrix.getRows() || cols != matrix.getCols())
		throw "Matrix::operator-(const Matrix &) - matrixes have different sizes";

	Matrix result(rows, cols);

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result.values[i][j] = values[i][j] - matrix.values[i][j];

	return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix &matrix) const {
	if (cols != matrix.rows)
		throw "Matrix::operator*(const Matrix &) - cols of first and rows of second matrixes does not match";

	Matrix result(rows, matrix.cols);

	for (int i = 0; i < result.rows; i++) {
		for (int j = 0; j < result.cols; j++) {
			T sum = 0;

			for (int k = 0; k < cols; k++)
				sum = sum + values[i][k] * matrix.values[k][j];

			result(i, j) = sum;
		}
	}

	return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T& value) const {
	Matrix result(rows, cols);

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result.values[i][j] = values[i][j] * value;

	return result;
}

// получение транспонированной матрицы
template <typename T>
Matrix<T> Matrix<T>::Transpose() const {
	Matrix<T> result(cols, rows);

	for (int i = 0; i < result.rows; i++) 
		for (int j = 0; j < result.cols; j++)
			result.values[i][j] = values[j][i];
	
	return result;
}

// получение обратной матрицы
template <typename T>
Matrix<T> Matrix<T>::Inverse() const {
	if (rows != cols) // если матрица не квадратная, то бросае исключение
		throw "Matrix::Inverse() - determinant of non squarable matrix does not exist";
	
	Matrix copy(*this); // копируем матрицу, чтобы не изменять исходную
	Matrix E(rows); // создаём квадратную матрицу

	// заполняем главную диагональ единицами, получая единичную матрицу
	for (int i = 0; i < rows; i++)
		E.values[i][i] = 1;

	// выполняем прямой ход метода Гаусса
	for (int j = 0; j < cols; j++) {
		if (copy.values[j][j] == 0) {
			int k = j;

			while (k < rows && copy.values[k][j] == 0)
				k++;

			if (k == rows)
				throw "Matrix::Inverse() - determinant = 0, inverse matrix does not exist";

			for (int l = 0; l < cols; l++) {
				T tmp = copy.values[j][l];
				copy.values[j][l] = copy.values[k][l];
				copy.values[k][l] = tmp;

				tmp = E.values[j][l];
				E.values[j][l] = E.values[k][l];
				E.values[k][l] = tmp;
			}
		}

		T v = copy.values[j][j];

		for (int k = 0; k < cols; k++) {
			copy.values[j][k] /= v;
			E.values[j][k] /= v;
		}

		for (int i = j + 1; i < rows; i++) {
			v = -copy.values[i][j];

			for (int j1 = 0; j1 < cols; j1++) {
				copy.values[i][j1] += copy.values[j][j1] * v;
				E.values[i][j1] += E.values[j][j1] * v;
			}
		}
	}

	// выполняем обратный ход метода Гаусса
	for (int j = cols - 1; j > 0; j--) {
		for (int i = j - 1; i >= 0; i--) {
			T v = -copy.values[i][j];

			for (int j1 = cols - 1; j1 >= 0; j1--) {
				copy.values[i][j1] += copy.values[j][j1] * v;
				E.values[i][j1] += E.values[j][j1] * v;
			}
		}
	}

	return E; // возвращаем матрицу Е, приведённую к обратной
}

// получение ортонормированной матрицы
template <typename T>
Matrix<T> Matrix<T>::Ortonormalize() const {
	Matrix Q(rows, cols); // ортонормированная мтарица
	Matrix B(rows, cols); // ортогональная матрица

	T beta_1 = 0; // вычисляем первый коэффициент

	for (int j = 0; j < rows; j++)
		beta_1 = beta_1 + values[0][j] * values[0][j]; // считаем сумму квадратов

	beta_1 = sqrt(beta_1); // извлекаем квадратный корень

	for (int j = 0; j < cols; j++) {
		B.values[0][j] = values[0][j]; // запоминаем исходый вектор
		Q.values[0][j] = values[0][j] / beta_1; // запоминаем нормированный вектор
	}

	// проводим процесс ортогонализации Грама-Шмидта по столбцам
	for (int i = 1; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			B.values[i][j] = values[i][j];

			for (int k = 0; k <= i - 1; k++) {
				T s = 0;
				T mod = 0;

				for (int m = 0; m < cols; m++) {
					s += values[i][m] * B.values[k][m];
					mod += B.values[k][m] * B.values[k][m];
				}

				B.values[i][j] -= s / mod * B.values[k][j];
			}
		}

		// аналогично находим нормирующий множитель
		T beta_i = 0;
		for (int j = 0; j < cols; j++)
			beta_i += B.values[i][j] * B.values[i][j];

		beta_i = sqrt(beta_i);

		for (int j = 0; j < cols; j++)
			Q.values[i][j] = B.values[i][j] / beta_i; // нормируем вектор
	}

	return Q; // возвращаем матрицу
}

// получение ортогональной матрицы
template <typename T>
Matrix<T> Matrix<T>::Ortogonalize() const {
	Matrix B(rows, cols); // ортогональная матрица

	for (int j = 0; j < cols; j++)
		B.values[0][j] = values[0][j]; // запоминаем исходый вектор

	// проводим процесс ортогонализации Грама-Шмидта по столбцам
	for (int i = 1; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			B.values[i][j] = values[i][j];

			for (int k = 0; k <= i - 1; k++) {
				T s = 0;
				T mod = 0;

				for (int m = 0; m < cols; m++) {
					s += values[i][m] * B.values[k][m];
					mod += B.values[k][m] * B.values[k][m];
				}

				B.values[i][j] -= s / mod * B.values[k][j];
			}
		}
	}

	return B; // возвращаем матрицу
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const Matrix<T> matrix) {
	for (int i = 0; i < matrix.rows; i++) {
		for (int j = 0; j < matrix.cols; j++)
			os << matrix.values[i][j] << " ";

		os << std::endl;
	}

	return os;
}

template <>
std::ostream& operator<<(std::ostream &os, const Matrix<double> matrix) {
	for (int i = 0; i < matrix.rows; i++) {
		for (int j = 0; j < matrix.cols; j++)
			os << std::setw(10) << (fabs(matrix.values[i][j]) < 1e-13 ? 0 : matrix.values[i][j]) << " ";

		os << std::endl;
	}

	return os;
}

#endif