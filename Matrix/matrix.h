// 基于C++11标准
#ifndef _MATRIX_
#define _MATRIX_

//#include <iostream>

#include <string>
#include <vector>
#include <tuple>
#include <memory>
#include <cmath>
#include <complex>
/**********************************  常用函数  *************************************************************************************************/
double sgn(const double& num);                                                 // 符号函数

																			   /***********************************  矩阵类  ********************************************/
class Matrix
{
private:
	int row;                                                                   // 行数
	int column;                                                                // 列数
	double* mat;                                                               // 采用一维数组存储元素

public:
	Matrix();
	Matrix(const int& m, const int& n);                                        // 常规构造: m行数，n列数
	Matrix(const std::initializer_list<double>& m);                            // 列表构造
	Matrix(const std::initializer_list<std::initializer_list<double>>& m);     // 列表构造
	Matrix(const Matrix& other);                                               // 拷贝构造
	Matrix(Matrix&& other);                                                    // 移动构造

																			   // 赋值运算
	Matrix& operator=(const Matrix& other)noexcept;
	Matrix& operator=(Matrix&& other)noexcept;                                 // 移动赋值
	Matrix& operator+=(const double& num);                                     // M += num
	Matrix& operator+=(const Matrix& m);                                       // M += M
	Matrix& operator-=(const double& num);                                     // M -= num
	Matrix& operator-=(const Matrix& m);                                       // M -= M
	Matrix& operator*=(const double& num);                                     // M *= num
	Matrix& operator*=(const Matrix& m);                                       // M *= M
	Matrix& operator/=(const double& num);

	// 下标索引
	inline double& operator()(int m, int n)noexcept;
	inline const double& operator()(int m, int n)const noexcept;
	inline double& at(int m, int n);
	inline const double& at(int m, int n)const;

	// 加法运算
	Matrix operator+()const;                                                   // +M
	friend Matrix operator+(const Matrix& m1, const Matrix& m2);               // M + M
	friend Matrix& operator+(const Matrix& m1, Matrix&& m2);                   // M + move(M)
	friend Matrix& operator+(Matrix&& m1, const Matrix& m2);                    // move(M) + M
	friend Matrix& operator+(Matrix&& m1, Matrix&& m2);                        // move(M) + move(M)
	friend Matrix operator+(const Matrix& m, const double& num);               // M + num
	friend Matrix operator+(const double& num, const Matrix& m);               // num + M
	friend Matrix& operator+(const double& num, Matrix&& m);                   // num + move(M)
	friend Matrix& operator+(Matrix&& m, const double& num);                   // move(M) + num

																			   // 减法运算
	Matrix operator-()const;                                                   // -M
	Matrix operator-(const Matrix& m)const;                                    // M - M
	Matrix& operator-(Matrix&& m)const;                                        // M - move(M)
	Matrix operator-(const double& num)const;                                  // M - num
	friend Matrix& operator-(Matrix&& m1, const Matrix& m2);                   // move(M) - M(or move(M))
	friend Matrix& operator-(Matrix&& m1, Matrix&& m2);                        // move(M) + move(M)
	friend Matrix operator-(const double& num, const Matrix& m);               // num - M
	friend Matrix& operator-(const double& num, Matrix&& m);                   // num - move(M)
	friend Matrix& operator-(Matrix&& m, const double& num);                   // move(M) - num

																			   // 乘法运算
	Matrix operator*(const Matrix& m)const;                                    // M * M
	Matrix operator*(const double& num)const;                                  // M * num
	friend Matrix& operator*(Matrix&& m, const double& num);                   // move(M) * num
	friend Matrix& operator*(Matrix&& m1, Matrix&& m2);                        // move(M) * move(M)
	friend Matrix operator*(const double& num, const Matrix& m);               // num * M
	friend Matrix& operator*(const double& num, Matrix&& m);                   // num * move(M)
	Matrix dotMult(const Matrix& m)const;                                      // M .* M
	Matrix& dotMult(Matrix&& m)const;                                          // M .* move(M)

																			   // 除法运算
	Matrix operator/(const double& num);                                       // M / num
	Matrix dotDiv(const Matrix& m)const;                                       // M ./ M
	Matrix& dotDiv(Matrix&& m)const;                                           // M ./ move(M)

																			   // 逻辑运算符
	bool operator==(const Matrix& m)const noexcept;                            // 等于
	bool operator!=(const Matrix& m)const noexcept;                            // 不等于

	const std::vector<int> size()const;                                        // 获取矩阵行列数
	std::vector<double> getRow(int n)const;                                    // 获取行向量
	std::vector<double> getColumn(int n)const;                                 // 获取列向量
	std::vector<double> getDiag()const;                                        // 获取对角元素
	double normOne()const;                                                     // 1范数,列和范数
	double normInf()const;                                                     // 无穷范数,行和范数
	Matrix trans()const;                                                       // 矩阵转置

	const int rank()const;                                                     // 矩阵的秩
	const double trace()const;                                                 // 矩阵的迹
	Matrix inv()const;                                                         // 逆矩阵
	const double det()const;                                                   // 矩阵行列式
	std::vector<Matrix> QR()const;                                             // 矩阵QR分解
	std::vector<std::complex<double>> eigs(double e = 0)const;                 // 矩阵特征值
	std::vector<Matrix> SVD()const;                                            // 奇异值分解，返回S、V、D三个矩阵
	Matrix subMat(int r1, int c1, int r2, int c2)const;                        // 子阵

	static Matrix eye(const int& m, const int& n);                             // 单位矩阵
	static Matrix ones(const int& m, const int& n);                            // 元素全为1的矩阵
	static Matrix diag(const std::initializer_list<double>& nums);
	static Matrix rbind(const std::initializer_list<Matrix>& M);               // [M1; M2; ...], 需要列数相等
	static Matrix cbind(const std::initializer_list<Matrix>& M);               // [M1, M2, ...], 需要行数相等

																			   // 全选主元高斯消去法
	static std::tuple<Matrix, std::unique_ptr<int[]>, std::unique_ptr<int[]>, int, Matrix>
		gaussElimination(const Matrix& A, const Matrix& b = Matrix());

	static std::tuple<Matrix, int> solve(const Matrix& A, const Matrix& b);    // 求解线性方程组 Ax = b

	virtual ~Matrix();                                                         // 析构

private:
	Matrix hess2()const;                                                       // 矩阵拟上三角分解(A = P*B*P'),P为正交矩阵,B为拟上三角阵
	std::vector<std::complex<double>>
		eigVal22(double a, double b, double c, double d)const;                 // 计算二阶矩阵特征值(按行输入)
	Matrix& iterM(Matrix& A)const;                                             // QR方法计算特征值中的矩阵迭代
																			   //void print()const;
};

/**************************************************************************************************/
// 常规构造: m行数，n列数
Matrix::Matrix() : row(0), column(0), mat(nullptr) {}
Matrix::Matrix(const int& m, const int& n) : row(0), column(0), mat(nullptr)
{
	if (m >= 0 && n >= 0)
	{
		if (m*n == 0)
		{
			mat = nullptr;
		}
		else
		{
			mat = new double[m*n]();
		}
		row = m;
		column = n;
	}
	else
	{
		throw std::length_error("Error of matrix dimension.");
	}
}
Matrix::Matrix(const std::initializer_list<double>& m) : row(0), column(0), mat(nullptr)
{
	int length = m.size();
	if (length > 0)
	{
		mat = new double[length];
		this->row = 1;
		this->column = length;

		for (int i = 0; i < length; i++)
		{
			mat[i] = *(m.begin() + i);
		}
	}
}
Matrix::Matrix(const std::initializer_list<std::initializer_list<double>>& m) :row(0), column(0), mat(nullptr)
{
	int row = m.size();
	if (row > 0)
	{
		bool dimMatch = true;
		int column = (*(m.begin())).size();

		// 判断各行元素数量是否相等，若相不等 dimMatch = false
		for (int i = 0; i < row; i++)
		{
			if (column != (*(m.begin() + i)).size())
				dimMatch = false;
		}

		if (dimMatch)
		{
			mat = new double[row * column];
			this->row = row;
			this->column = column;
			for (int i = 0; i < row; i++)
			{
				for (int j = 0; j < column; j++)
				{
					mat[i*column + j] = *((*(m.begin() + i)).begin() + j);
				}
			}
		}
		else
		{
			throw std::length_error("Dimensions do not match.");
		}
	}
}

// 复制构造
Matrix::Matrix(const Matrix& other) : row(0), column(0), mat(nullptr)
{
	mat = new double[other.row * other.column];
	row = other.row;
	column = other.column;

	int size = row * column;
	for (int i = 0; i < size; i++)
		mat[i] = other.mat[i];
}

// 移动构造
Matrix::Matrix(Matrix&& other) : row(0), column(0), mat(nullptr)
{
	//*this = std::move(other);
	mat = other.mat;
	row = other.row;
	column = other.column;
	other.mat = nullptr;
}

// 重载运算符(),用于获取矩阵元素和赋值,索引从0开始
inline double& Matrix::operator()(int m, int n)noexcept
{
	return mat[column*m + n];
}
inline const double& Matrix::operator()(int m, int n)const noexcept
{
	return mat[column*m + n];
}
inline double& Matrix::at(int m, int n)
{
	if (m >= 0 && n >= 0 && m < row && n < column)
	{
		return mat[column*m + n];
	}
	else
	{
		throw std::out_of_range("In function at():Out of index range of the matrix.");
	}
}
inline const double& Matrix::at(int m, int n)const
{
	if (m >= 0 && n >= 0 && m < row && n < column)
	{
		return mat[column*m + n];
	}
	else
	{
		throw std::out_of_range("In function at():Out of index range of the matrix.");
	}
}

// 赋值运算
Matrix& Matrix::operator=(const Matrix& other)noexcept
{
	if (this != &other)
	{
		delete[] this->mat;
		int size = other.row * other.column;
		this->mat = new double[size];
		this->row = other.row;
		this->column = other.column;

		for (int i = 0; i < size; i++)
			this->mat[i] = other.mat[i];
	}
	return *this;
}
Matrix& Matrix::operator=(Matrix&& other)noexcept
{
	if (this != &other)
	{
		delete[] this->mat;
		this->mat = other.mat;
		this->row = other.row;
		this->column = other.column;
		other.mat = nullptr;
	}
	return *this;
}
Matrix& Matrix::operator+=(const Matrix& m)
{
	if (row == m.row && column == m.column)
	{
		int size = row * column;
		for (int i = 0; i < size; i++)
		{
			mat[i] += m.mat[i];
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return *this;
}
Matrix& Matrix::operator+=(const double& num)
{
	if (row > 0 && column > 0)
	{
		int size = row * column;
		for (int i = 0; i < size; i++)
		{
			mat[i] += num;
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return *this;
}
Matrix& Matrix::operator-=(const double& num)
{
	if (row > 0 && column > 0)
	{
		int size = row * column;
		for (int i = 0; i < size; i++)
		{
			mat[i] -= num;
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return *this;
}
Matrix& Matrix::operator-=(const Matrix& m)
{
	if (row == m.row && column == m.column)
	{
		int size = row * column;
		for (int i = 0; i < size; i++)
		{
			mat[i] -= m.mat[i];
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return *this;
}
Matrix& Matrix::operator*=(const double& num)
{
	if (row > 0)
	{
		int size = row * column;
		for (int i = 0; i < size; i++)
		{
			mat[i] *= num;
		}
	}
	else
	{
		throw std::length_error("Null Matrix.");
	}
	return *this;
}
Matrix& Matrix::operator*=(const Matrix& m)
{
	if (this->size()[1] == m.size()[0])
	{
		int row = this->size()[0];
		int column = m.size()[1];
		Matrix matrix(row, column);

		int count = m.size()[0];
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < column; j++)
			{
				double sum = 0;
				for (int k = 0; k < count; k++)
				{
					sum += (*this)(i, k) * m(k, j);
				}
				matrix(i, j) = sum;
			}
		}
		*this = std::move(matrix);
	}
	else
	{
		throw std::length_error("Dimensions does not match.");
	}
	return *this;
}
Matrix& Matrix::operator/=(const double& num)
{
	if (num == 0)
	{
		throw std::invalid_argument("Error of Inf.");
	}
	else
	{
		if (row > 0)
		{
			int size = row * column;
			for (int i = 0; i < size; i++)
			{
				mat[i] /= num;
			}
		}
		else
		{
			throw std::length_error("Null Matrix.");
		}
	}
	return *this;
}

// 加法运算
Matrix operator+(const Matrix& m1, const Matrix& m2)
{
	if (m1.row == m2.row && m1.column == m2.column)
	{
		int size = m1.row * m1.column;
		Matrix matrix(m1.row, m1.column);
		for (int i = 0; i < size; i++)
		{
			matrix.mat[i] = m1.mat[i] + m2.mat[i];
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m1;
}
Matrix& operator+(const Matrix& m1, Matrix&& m2)
{
	if (m1.row == m2.row && m1.column == m2.column)
	{
		int size = m1.row * m1.column;
		for (int i = 0; i < size; i++)
		{
			m2.mat[i] += m1.mat[i];
		}
		return m2;
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m2 = m1;
}
Matrix operator+(const Matrix& m, const double& num)
{
	if (m.row > 0)
	{
		Matrix matrix(m.row, m.column);
		int size = m.row * m.column;
		for (int i = 0; i < size; i++)
		{
			matrix.mat[i] = m.mat[i] + num;
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Null Matrix.");
	}
	return m;
}
Matrix Matrix::operator+()const { return *this; }
Matrix& operator+(Matrix&& m1, const Matrix& m2)
{
	if (m1.row == m2.row && m1.column == m2.column)
	{
		int size = m1.row * m1.column;
		for (int i = 0; i < size; i++)
		{
			m1.mat[i] += m2.mat[i];
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m1;
}
Matrix& operator+(Matrix&& m1, Matrix&& m2)
{
	return std::move(m1) + m2;
}
Matrix operator+(const double& num, const Matrix& m)
{
	if (m.row > 0)
	{
		Matrix matrix(m.row, m.column);
		int size = m.row * m.column;
		for (int i = 0; i < size; i++)
		{
			matrix.mat[i] = m.mat[i] + num;
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return Matrix();
}
Matrix& operator+(const double& num, Matrix&& m)
{
	if (m.row > 0)
	{
		int size = m.row * m.column;
		for (int i = 0; i < size; i++)
		{
			m.mat[i] += num;
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m;
}
Matrix& operator+(Matrix&& m, const double& num)
{
	if (m.row > 0)
	{
		int size = m.row * m.column;
		for (int i = 0; i < size; i++)
		{
			m.mat[i] += num;
		}
	}
	else
	{
		throw std::length_error("Null Matrix.");
	}
	return m;
}

// 减法运算
Matrix Matrix::operator-(const Matrix& m)const
{
	if (m.row == row && m.column == column)
	{
		int size = row * column;
		Matrix matrix(row, column);
		for (int i = 0; i < size; i++)
		{
			matrix.mat[i] = this->mat[i] - m.mat[i];
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return *this;
}
Matrix& Matrix::operator-(Matrix&& m)const
{
	if (m.row == row && m.column == column)
	{
		int size = row * column;
		for (int i = 0; i < size; i++)
		{
			m.mat[i] = this->mat[i] - m.mat[i];
		}
		return m;
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m = *this;
}
Matrix Matrix::operator-(const double& num)const
{
	if (this->row > 0)
	{
		Matrix matrix(row, column);
		int size = row * column;
		for (int i = 0; i < size; i++)
		{
			matrix.mat[i] = this->mat[i] - num;
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Null Matrix.");
	}
	return *this;
}
Matrix Matrix::operator-()const
{
	Matrix matrix(row, column);
	int size = row * column;
	for (int i = 0; i < size; i++)
	{
		matrix.mat[i] = -(this->mat[i]);
	}
	return matrix;
}
Matrix& operator-(Matrix&& m1, const Matrix& m2)
{
	if (m1.row == m2.row && m1.column == m2.column)
	{
		int size = m1.row * m1.column;
		for (int i = 0; i < size; i++)
		{
			m1.mat[i] -= m2.mat[i];
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m1;
}
Matrix& operator-(Matrix&& m1, Matrix&& m2)
{
	return std::move(m1) - m2;
}
Matrix operator-(const double& num, const Matrix& m)
{
	if (m.row > 0)
	{
		Matrix matrix(m.row, m.column);
		int size = m.row * m.column;
		for (int i = 0; i < size; i++)
		{
			matrix.mat[i] = num - m.mat[i];
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return Matrix();
}
Matrix& operator-(const double& num, Matrix&& m)
{
	if (m.row > 0)
	{
		int size = m.row * m.column;
		for (int i = 0; i < size; i++)
		{
			m.mat[i] = num - m.mat[i];
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m;
}
Matrix& operator-(Matrix&& m, const double& num)
{
	if (m.row > 0)
	{
		int size = m.row * m.column;
		for (int i = 0; i < size; i++)
		{
			m.mat[i] -= num;
		}
	}
	else
	{
		throw std::length_error("Null Matrix.");
	}
	return m;
}

// 乘法运算
Matrix Matrix::operator*(const Matrix& m)const
{
	if (this->size()[1] == m.size()[0])
	{
		int row = this->size()[0];
		int column = m.size()[1];
		Matrix matrix(row, column);

		int count = m.size()[0];
		for (int i = 0; i < row; i++)
		{
			for (int j = 0; j < column; j++)
			{
				double sum = 0;
				for (int k = 0; k < count; k++)
				{
					sum += (*this)(i, k) * m(k, j);
				}
				matrix(i, j) = sum;
			}
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Dimensions does not match.");
	}
	return Matrix();
}
Matrix Matrix::operator*(const double& num)const
{
	Matrix matrix(row, column);
	int size = row * column;
	for (int i = 0; i < size; i++)
	{
		matrix.mat[i] = mat[i] * num;
	}
	return matrix;
}
Matrix& operator*(Matrix&& m, const double& num)
{
	int size = m.row * m.column;
	for (int i = 0; i < size; i++)
	{
		m.mat[i] *= num;
	}
	return m;
}
Matrix& operator*(Matrix&& m1, Matrix&& m2)
{
	return m1 = std::move(m1) * m2;
}
Matrix operator*(const double& num, const Matrix& m)
{
	Matrix matrix(m.row, m.column);
	int size = m.row * m.column;
	for (int i = 0; i < size; i++)
	{
		matrix.mat[i] = num * m.mat[i];
	}
	return matrix;
}
Matrix& operator*(const double& num, Matrix&& m)
{
	int size = m.row * m.column;
	for (int i = 0; i < size; i++)
	{
		m.mat[i] *= num;
	}
	return m;
}
Matrix Matrix::dotMult(const Matrix& m)const
{
	if (this->row == m.row && this->column == m.column)
	{
		Matrix matrix(this->row, this->column);
		int size = this->row * this->column;
		for (int i = 0; i < size; i++)
		{
			matrix.mat[i] = this->mat[i] * m.mat[i];
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return *this;
}
Matrix& Matrix::dotMult(Matrix&& m)const
{
	if (this->row == m.row && this->column == m.column)
	{
		int size = this->row * this->column;
		for (int i = 0; i < size; i++)
		{
			m.mat[i] *= this->mat[i];
		}
		return m;
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m = *this;
}

// 除法运算
Matrix Matrix::operator/(const double& num)
{
	if (this->row > 0)
	{
		if (num != 0)
		{
			Matrix matrix(row, column);
			int size = row * column;
			for (int i = 0; i < size; i++)
			{
				matrix.mat[i] = this->mat[i] / num;
			}
			return matrix;
		}
		else
		{
			throw std::invalid_argument("Error of Inf.");
		}
	}
	else
	{
		throw std::length_error("Null Matrix.");
	}
	return *this;
}
Matrix Matrix::dotDiv(const Matrix& m)const
{
	if (this->row == m.row && this->column == m.column)
	{
		bool flag = true;
		int size = this->row * this->column;
		for (int i = 0; i < size; i++)
		{
			if (m.mat[i] == 0)
			{
				flag = false;
			}
		}

		if (flag)
		{
			Matrix matrix(this->row, this->column);
			for (int i = 0; i < size; i++)
			{
				matrix.mat[i] = this->mat[i] / m.mat[i];
			}
			return matrix;
		}
		else
		{
			throw std::invalid_argument("Error of Inf.");
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return *this;
}
Matrix& Matrix::dotDiv(Matrix&& m)const
{
	if (this->row == m.row && this->column == m.column)
	{
		bool flag = true;
		int size = this->row * this->column;
		for (int i = 0; i < size; i++)
		{
			if (m.mat[i] == 0)
			{
				flag = false;
			}
		}

		if (flag)
		{
			for (int i = 0; i < size; i++)
			{
				m.mat[i] = this->mat[i] / m.mat[i];
			}
			return m;
		}
		else
		{
			throw std::invalid_argument("Error of Inf.");
		}
	}
	else
	{
		throw std::length_error("Dimensions do not match.");
	}
	return m = *this;
}

// 逻辑运算
bool Matrix::operator==(const Matrix& m)const noexcept
{
	if (this->row == m.row && this->column == m.column)
	{
		bool flag = true;
		int size = this->row * this->column;
		for (int i = 0; i < size; i++)
		{
			if (this->mat[i] != m.mat[i])
			{
				flag = false;
				break;
			}
		}
		return flag;
	}
	return false;
}
bool Matrix::operator!=(const Matrix& m)const noexcept
{
	return !(*this == m);
}

// 获取矩阵行列数
const std::vector<int> Matrix::size()const
{
	return std::vector<int>{row, column};
}

// 获取行向量
std::vector<double> Matrix::getRow(int n)const
{
	if (n >= 0 && n < this->row)
	{
		if (this->row == 0)
		{
			return std::vector<double>();
		}
		else
		{
			std::vector<double> vec;
			for (int i = 0; i < this->column; i++)
			{
				vec.emplace_back(mat[n*column + i]);
			}
			return vec;
		}
	}
	else
	{
		throw std::out_of_range("Out of index range.");
	}
	return std::vector<double>();
}

// 获取列向量
std::vector<double> Matrix::getColumn(int n)const
{
	if (n >= 0 && n < this->column)
	{
		if (this->column == 0)
		{
			return std::vector<double>();
		}
		else
		{
			std::vector<double> vec;
			for (int i = 0; i < this->row; i++)
			{
				vec.emplace_back(mat[i*column + n]);
			}
			return vec;
		}
	}
	else
	{
		throw std::out_of_range("Out of index range.");
	}
	return std::vector<double>();
}

// 获取对角元素
std::vector<double> Matrix::getDiag()const
{
	if (this->row == 0)
		return std::vector<double>();

	int size;
	this->row < this->column ? size = this->row : size = this->column;
	std::vector<double> vec;
	for (int i = 0; i < size; i++)
	{
		vec.emplace_back(mat[i*column + i]);
	}
	return vec;
}

// 1范数(列和范数)
double Matrix::normOne()const
{
	double norm = 0;
	if (this->row > 0)
	{
		for (int i = 0; i < this->column; i++)
		{
			double temp = 0;
			for (int j = 0; j < this->row; j++)
			{
				temp += abs(this->mat[j*column + i]);
			}
			if (temp > norm)
			{
				norm = temp;
			}
		}
	}
	return norm;
}

// 无穷范数(行和范数)
double Matrix::normInf()const
{
	double norm = 0;
	if (this->row > 0)
	{
		for (int i = 0; i < this->row; i++)
		{
			double temp = 0;
			for (int j = 0; j < this->column; j++)
			{
				temp += abs(this->mat[i*column + j]);
			}
			if (temp > norm)
			{
				norm = temp;
			}
		}
	}
	return norm;
}

// 矩阵转置
Matrix Matrix::trans()const
{
	if (row > 0)
	{
		Matrix matrix(this->column, this->row);
		for (int i = 0; i < column; i++)
		{
			for (int j = 0; j < row; j++)
			{
				matrix.mat[i*matrix.column + j] = this->mat[j*column + i];
			}
		}
		return matrix;
	}
	return Matrix();
}

// 矩阵的秩
const int Matrix::rank()const { return std::get<3>(gaussElimination(*this)); }

// 矩阵的迹
const double Matrix::trace()const
{
	if (this->row == this->column)
	{
		if (this->row == 0)
		{
			return 0; // 空矩阵的迹为0
		}
		else
		{
			double sum = 0;
			for (int i = 0; i < this->row; i++)
			{
				sum += this->mat[i*column + i];
			}
			return sum;
		}
	}
	else
	{
		throw std::length_error("The matrix must be a square matrix.");
	}
	return 0;
}

// 行列式
const double Matrix::det()const
{
	// 将矩阵化为上三角矩阵，对角元素乘积即为行列式的值
	// 正负号由行交换次数确定
	if (this->row == this->column)
	{
		if (this->row == 0)       // 空矩阵行列式为 1
			return 1;
		if (this->row == 1)
			return this->mat[0];

		int flag = 1;
		double temp = 0;
		int col = 0;
		Matrix matrix(*this);
		for (int i = 0; i < matrix.row - 1; i++)
		{
			// 选列主元
			temp = abs(matrix.mat[i*column + i]);
			col = i;
			for (int j = i + 1; j < matrix.row; j++)
			{
				if (abs(matrix.mat[j*column + i]) > temp)
				{
					temp = abs(matrix.mat[j*column + i]);
					col = j; // 记录主元所在行
				}
			}
			// 行交换
			if (col != i)
			{
				for (int k = 0; k < matrix.column; k++)
				{
					temp = matrix.mat[i*column + k];
					matrix.mat[i*column + k] = matrix.mat[col*column + k];
					matrix.mat[col*column + k] = temp;
				}
				flag *= -1;
				// 主元为0，说明方阵不满秩，行列式为0
				if (matrix.mat[i*column + i] == 0)
					return 0;
			}
			// 高斯消元
			for (int j = i + 1; j < matrix.row; j++)
			{
				temp = matrix.mat[j*column + i] / matrix.mat[i*column + i];
				for (int k = i; k < matrix.column; k++)
				{
					matrix.mat[j*column + k] -= temp * matrix.mat[i*column + k];
				}
			}
		}
		// 对角元素的乘积
		temp = 1;
		for (int i = 0; i < matrix.row; i++)
		{
			temp *= matrix.mat[i*column + i];
		}
		return temp*flag;
	}
	else
	{
		throw std::length_error("The matrix must be a square matrix.");
	}
	return 0;
}

/* 矩阵QR分解(A = QR),Q为正交矩阵,R为上三角矩阵(且主对角线元素>=0时,
分解具有唯一性)。若R满秩，主对角线元素 > 0 */
std::vector<Matrix> Matrix::QR()const
{
	if (this->row > 0 && this->row == this->column)
	{
		Matrix A(*this);
		Matrix Q = Matrix::eye(this->row, this->row);

		int n = A.row;   // 矩阵维数
		double d_r = 0;
		double c_r = 0;
		double h_r = 0;
		double* u_r = new double[n];
		Matrix U_R;
		Matrix P_R;

		for (int r = 0; r < n - 1; r++)
		{
			// 判断A(r+1,r) ~ A(n,r)是否全为0
			bool zero = true;
			for (int i = r + 1; i < n; i++)
			{
				if (A(i, r) != 0)
				{
					zero = false;
					break;
				}
			}

			if (!zero) // 如果A(r+1,r) ~ A(n,r)不全为0
			{
				// 计算d_r
				d_r = 0;
				for (int i = r; i < n; i++)
				{
					d_r += A(i, r) * A(i, r);
				}
				d_r = sqrt(d_r);

				// 计算c_r
				c_r = -sgn(A(r, r)) * d_r;
				// 计算h_r
				h_r = c_r * (c_r - A(r, r));
				// 计算向量u_r
				for (int i = 0; i < n; i++)
				{
					if (i < r)
						u_r[i] = 0;
					else
						u_r[i] = A(i, r);
				}
				u_r[r] -= c_r;

				// 利用数组u_r构造列向量U_R
				U_R.row = A.row;
				U_R.column = 1;
				U_R.mat = u_r;

				// 进行一次Q、A的迭代
				Q -= Q * U_R * U_R.trans() / h_r;
				P_R = A.trans() * U_R / h_r;
				A -= U_R * P_R.trans();
				U_R.mat = nullptr;

				// 消除上三角阵R主对角线以下元素的截断误差
				for (int i = r + 1; i < n; i++)
				{
					A(i, r) = 0;
				}
			}
		}
		delete[] u_r;

		// 迭代完成后 Q = Q_r, R = A_r
		return std::vector<Matrix>({ Q, A });
	}
	else
	{
		throw std::length_error("The matrix must be a square matrix.");
	}
	return std::vector<Matrix>();
}

// 求矩阵全部特征值(带双步位移的QR方法), e为精度水平
std::vector<std::complex<double>> Matrix::eigs(double e)const
{
	// 只有非空方阵才能求特征值
	if (this->row != this->column || this->row == 0)
	{
		if (this->row != this->column)
		{
			throw std::length_error("In Matrix::eigs: The matrix must be a square matrix.");
		}
		return std::vector<std::complex<double>>();
	}

	std::vector<std::complex<double>> vec;                     // 存储特征值的容器
	Matrix matrix(this->hess2());                              // 矩阵的拟上三角化
	e == 0 ? e = this->normOne() / 1.0e10 : e = abs(e);        // 确定精度水平,默认为矩阵1范数的1/(1e10)
	int m = this->row;

	while (true)
	{
		if (m == 1)
		{
			vec.emplace_back(matrix(0, 0), 0);
			break;
		}

		if (abs(matrix(m - 1, m - 2)) < e)
		{
			vec.emplace_back(matrix(m - 1, m - 1), 0);
			m = m - 1;
		}
		else
		{
			if (m == 2 || abs(matrix(m - 2, m - 3)) < e)
			{
				auto val = eigVal22(matrix(m - 2, m - 2), matrix(m - 2, m - 1), matrix(m - 1, m - 2), matrix(m - 1, m - 1));
				vec.push_back(val[0]);
				vec.push_back(val[1]);
				if (m == 2)
				{
					break;
				}
				m = m - 2;
			}
			else
			{
				// 计算A_k+1
				// iterM：QR方法计算特征值时的矩阵迭代
				matrix = iterM(matrix.subMat(0, 0, m - 1, m - 1));
			}
		}
	}
	return vec;
}

// 矩阵求逆
Matrix Matrix::inv()const
{
	if (this->row > 0 && this->row == this->column)
	{
		Matrix matrix(*this);
		int n = matrix.row;
		std::unique_ptr<int[]> rFlag(new int[n]);
		std::unique_ptr<int[]> cFlag(new int[n]);

		double temp;
		for (int k = 0; k < n; k++)
		{
			temp = abs(matrix.mat[k * column + k]);
			rFlag[k] = cFlag[k] = k;

			// 搜寻子阵的最大主元并记录其位置
			for (int i = k; i < n; i++)
			{
				for (int j = k; j < n; j++)
				{
					if (abs(matrix.mat[i * column + j]) > temp)
					{
						temp = abs(matrix.mat[i * column + j]);
						rFlag[k] = i;
						cFlag[k] = j;
					}
				}
			}

			// 全选主元交换行和列
			// 列交换
			if (cFlag[k] != k)
			{
				for (int i = 0; i < n; i++)
				{
					temp = matrix.mat[i * column + k];
					matrix.mat[i * column + k] = matrix.mat[i * column + cFlag[k]];
					matrix.mat[i * column + cFlag[k]] = temp;
				}
			}

			// 行交换
			if (rFlag[k] != k)
			{
				for (int i = 0; i < n; i++)
				{
					temp = matrix.mat[k * column + i];
					matrix.mat[k * column + i] = matrix.mat[rFlag[k] * column + i];
					matrix.mat[rFlag[k] * column + i] = temp;
				}
			}

			// 主元为0, 矩阵不满秩, 没有逆矩阵
			if (matrix.mat[k*column + k] == 0)
			{
				throw std::logic_error("No inverse matrix.");
				return Matrix();
			}
			else
			{
				matrix.mat[k*column + k] = 1 / matrix.mat[k*column + k];
			}

			for (int j = 0; j < n; j++)
				if (j != k)
					matrix.mat[k * column + j] *= matrix.mat[k * column + k];

			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					if (i != k && j != k)
						matrix.mat[i * column + j] -= matrix.mat[i * column + k] * matrix.mat[k * column + j];

			for (int i = 0; i < n; i++)
				if (i != k)
					matrix.mat[i * column + k] = -matrix.mat[i * column + k] * matrix.mat[k * column + k];
		}

		// 还原矩阵: (1) 先交换的行(列)后进行恢复; (2) 原来的行(列)交换用列(行)来恢复
		for (int k = n - 1; k >= 0; k--)
		{
			// 列交换
			if (rFlag[k] != k)
			{
				for (int i = 0; i < n; i++)
				{
					temp = matrix.mat[i * column + k];
					matrix.mat[i * column + k] = matrix.mat[i * column + rFlag[k]];
					matrix.mat[i * column + rFlag[k]] = temp;
				}
			}

			// 行交换
			if (cFlag[k] != k)
			{
				for (int i = 0; i < n; i++)
				{
					temp = matrix.mat[k * column + i];
					matrix.mat[k * column + i] = matrix.mat[cFlag[k] * column + i];
					matrix.mat[cFlag[k] * column + i] = temp;
				}
			}
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Number of rows must equal columns.");
	}
	return *this;
}

// 子阵,形参为左上角元素位置和右下角元素位置
Matrix Matrix::subMat(int r1, int c1, int r2, int c2)const
{
	if (r1 >= 0 && c1 >= 0 && r2 < this->row && c2 < this->column && r1 <= r2 && c1 <= c2)
	{
		Matrix matrix(r2 - r1 + 1, c2 - c1 + 1);
		for (int i = r1; i <= r2; i++)
		{
			for (int j = c1; j <= c2; j++)
			{
				matrix.mat[(i - r1)*matrix.column + j - c1] = this->mat[i*column + j];
			}
		}
		return matrix;
	}
	else
	{
		throw std::out_of_range("Invalid index of the matrix.");
	}
	return *this;
}

// 静态函数
Matrix Matrix::eye(const int& m, const int& n)
{
	if (m > 0 && n > 0)
	{
		Matrix matrix(m, n);
		int size = 0;
		m < n ? size = m : size = n;
		for (int i = 0; i < size; i++)
		{
			matrix.mat[i*n + i] = 1;
		}
		return matrix;
	}
	else
	{
		throw std::length_error("Error of dimension size.");
	}
	return Matrix();
}
Matrix Matrix::ones(const int& m, const int& n) { return Matrix(m, n) + 1; }
Matrix Matrix::diag(const std::initializer_list<double>& nums)
{
	int size = nums.size();
	Matrix matrix(size, size);
	for (int i = 0; i < size; i++)
	{
		matrix.mat[i*matrix.column + i] = *(nums.begin() + i);
	}
	return matrix;
}
Matrix Matrix::rbind(const std::initializer_list<Matrix>& M)
{
	int count = M.size();
	if (count > 0)
	{
		int column = (*(M.begin())).column;
		int row = (*(M.begin())).row;

		// 判断矩阵是否可拼接,并计算拼接后矩阵的行数
		for (int i = 1; i < count; i++)
		{
			row += (*(M.begin() + i)).row;

			if ((*(M.begin() + i)).column != column)
			{
				throw std::length_error("Dimensions do not match.");
				return Matrix();
			}
		}

		// 拼接矩阵 [M1; M2; ...]
		Matrix matrix(row, column);
		int p = 0;
		for (int i = 0; i < count; i++)
		{
			int size = (*(M.begin() + i)).row * (*(M.begin() + i)).column;
			for (int j = 0; j < size; j++)
			{
				matrix.mat[p + j] = (*(M.begin() + i)).mat[j];
			}
			p += size;
		}

		return matrix;
	}

	return Matrix();
}
Matrix Matrix::cbind(const std::initializer_list<Matrix>& M)
{
	int count = M.size();
	if (count > 0)
	{
		int row = (*(M.begin())).row;
		int column = (*(M.begin())).column;

		// 判断矩阵是否可拼接,并计算拼接后矩阵的行数
		for (int i = 1; i < count; i++)
		{
			column += (*(M.begin() + i)).column;

			if ((*(M.begin() + i)).row != row)
			{
				throw std::length_error("Dimensions do not match.");
				return Matrix();
			}
		}

		// 拼接矩阵 [M1, M2, ...]
		Matrix matrix(row, column);
		int p = 0;
		for (int k = 0; k < row; k++)
		{
			for (int i = 0; i < count; i++)
			{
				int col = (*(M.begin() + i)).column;
				for (int j = 0; j < col; j++)
				{
					matrix.mat[p + j] = (*(M.begin() + i)).mat[k*col + j];
				}
				p += col;
			}
		}

		return matrix;
	}

	return Matrix();
}
// 全选主元高斯消去
// 返回值分别为：{消元后的矩阵，行交换记录，列交换记录，矩阵的秩}
std::tuple<Matrix, std::unique_ptr<int[]>, std::unique_ptr<int[]>, int, Matrix>
Matrix::gaussElimination(const Matrix& A, const Matrix& b)
{
	/***************************************************************/
	// 形参vec[可选]，默认为空向量
	// 若vec不为空，则需传入线性方程组 Ax = b中的向量b，
	// 函数返回值中的第四个参数为向量b通过高斯消元过程中的行变换得到的新向量
	/***************************************************************/

	// 空矩阵秩为0
	if (A.row == 0)
		return std::make_tuple(Matrix(), nullptr, nullptr, 0, Matrix());

	// 形参vec为一向量，与消元矩阵具有相同的行数
	if (b != Matrix())
	{
		if (b.column != 1 || b.row != A.row)
		{
			throw std::length_error("The length of passed in vector should equal the row of the matrix.");
		}
	}

	Matrix matrix(A);
	Matrix vec(b);
	int size = matrix.row < matrix.column ? matrix.row : matrix.column;
	int row = matrix.row;
	int column = matrix.column;

	// 记录主元位置
	std::unique_ptr<int[]> R(new int[size]());          // 行交换信息
	std::unique_ptr<int[]> C(new int[size]());          // 列交换信息

	double temp;
	// 全选主元消元
	for (int k = 0; k < size; k++)
	{
		temp = abs(matrix.mat[k * column + k]);
		R[k] = C[k] = k;

		// 搜寻子阵的最大主元并记录其位置
		for (int i = k; i < matrix.row; i++)
		{
			for (int j = k; j < matrix.column; j++)
			{
				if (abs(matrix.mat[i * column + j]) > temp)
				{
					temp = abs(matrix.mat[i * column + j]);
					R[k] = i;
					C[k] = j;
				}
			}
		}

		// 全选主元交换行和列
		// 列交换
		if (C[k] != k)
		{
			for (int i = 0; i < matrix.row; i++)
			{
				temp = matrix.mat[i * column + k];
				matrix.mat[i * column + k] = matrix.mat[i * column + C[k]];
				matrix.mat[i * column + C[k]] = temp;
			}
		}

		// 行交换
		if (R[k] != k)
		{
			// 矩阵matrix进行行交换
			for (int i = 0; i < matrix.column; i++)
			{
				temp = matrix.mat[k * column + i];
				matrix.mat[k * column + i] = matrix.mat[R[k] * column + i];
				matrix.mat[R[k] * column + i] = temp;
			}

			// 向量vec进行行交换
			if (vec != Matrix())
			{
				temp = vec.mat[k];
				vec.mat[k] = vec.mat[R[k]];
				vec.mat[R[k]] = temp;
			}
		}

		// 主元为0则主元右下角元素均为0, 高斯消元完成
		if (matrix.mat[k*column + k] == 0)
			break;

		// 消元
		for (int i = k + 1; i < matrix.row; i++)
		{
			// 矩阵matrix消元
			temp = matrix.mat[i*column + k] / matrix.mat[k*column + k];
			matrix.mat[i*column + k] = 0;
			for (int j = k + 1; j < matrix.column; j++)
			{
				matrix.mat[i*column + j] -= temp * matrix.mat[k*column + j];
			}

			// 向量vec消元
			if (vec != Matrix())
			{
				vec.mat[i] -= temp * vec.mat[k];
			}
		}
	}

	// 主对角线非零元素个数为矩阵的秩
	int count = 0;
	for (int i = 0; i < size; i++)
	{
		if (matrix.mat[i*column + i] != 0)
			count++;
	}

	return std::make_tuple(matrix, std::move(R), std::move(C), count, vec);
}

// 求解线性方程组: Ax = b
std::tuple<Matrix, int> Matrix::solve(const Matrix& A, const Matrix& b)
{
	/*******************************************************/
	// 求解线性方程组: Ax = b
	// A为系数矩阵，b为目标向量
	// 采用全选主元的高斯消元法求解方程组

	// (1)若方程有唯一解：   返回 tuple(解向量x,  标志量0)
	// (2)若方程有无穷多解： 返回 tuple(一个特解x, 标志量1)
	// (3)若方程无精确解：   返回 tuple(近似解x*,  标志量2)
	//    近似解x*使得|Ax* - b|最小
	/********************************************************/

	// 判断系数矩阵的行数宇目标向量的行数是否相等
	if (A.row != b.row)
	{
		throw std::length_error("Size of coefficient matrix and target vector does not match.");
	}

	auto temp = gaussElimination(A, b);                // 对系数矩阵进行全选主元的高斯消元
	Matrix& matrix = std::get<0>(temp);                // 全选主元高斯消元后的矩阵
	std::unique_ptr<int[]>& R = std::get<1>(temp);     // 行交换信息
	std::unique_ptr<int[]>& C = std::get<2>(temp);     // 列交换信息
	int rankA = std::get<3>(temp);                     // 系数矩阵的秩rank(A)
	int rankAb = cbind({ A,b }).rank();                // rank(A,b)
	int n = A.column;                                  // 未知数个数
	Matrix& vec = std::get<4>(temp);                   // 高斯消元后方程组的目标向量

													   // rank(A) == rank(A,b) && rank(A) == n: 方程有唯一解
													   // rank(A) == rank(A,b) && rank(A) < n:  方程有无穷多解
													   // others: 方程无解( rank(A)不可能大于n )

													   // (1) 方程有解
	if (rankA == rankAb)
	{
		Matrix x(n, 1);   // 方程的解

						  // 求解高斯消元后的上三角方程组
		for (int i = rankA - 1; i >= 0; i--)
		{
			for (int j = 0; j < rankA - i - 1; j++)
			{
				vec(i, 0) -= matrix(i, rankA - 1 - j) * x(rankA - 1 - j, 0);
			}

			x(i, 0) = 1.0 * vec(i, 0) / matrix(i, i);
		}

		// 方程有无穷多个解时，另解向量x部分元素为0，得到一个特解
		if (rankA < n)
		{
			for (int i = rankA; i < n; i++)
			{
				x(i, 0) = 0;
			}
		}

		// 依据列交换信息,恢复向量x中元素的顺序
		auto temp = x(rankA - 1, 0);
		for (int i = rankA - 1; i >= 0; i--)
		{
			if (C[i] != i)
			{
				temp = x(i, 0);
				x(i, 0) = x(C[i], 0);
				x(C[i], 0) = temp;
			}
		}

		// (1.1) 方程有唯一解
		if (rankA == n)
		{
			return std::make_tuple(x, 0);
		}

		// (1.2) 方程有无穷多个解
		if (rankA < n)
		{
			return std::make_tuple(x, 1);
		}
	}

	// (2)方程无解
	if (rankA != rankAb)
	{
		Matrix AT(A.trans());
		return std::make_tuple(std::get<0>(solve(AT*A, AT*b)), 2);
	}

	return std::make_tuple(Matrix(), 2);
}

// 析构
Matrix::~Matrix() { delete[] mat; }

// 矩阵拟上三角化
Matrix Matrix::hess2()const
{
	if (this->row > 0 && this->row == this->column)
	{
		Matrix A(*this);

		int n = A.row;   // 矩阵维数
		double d_r = 0;
		double c_r = 0;
		double h_r = 0;
		double* u_r = new double[n];
		Matrix U_R;
		Matrix P_R;
		Matrix Q_R;
		double t_r;

		for (int r = 0; r < n - 2; r++)
		{
			// 判断A(r+2,r) ~ A(n,r)是否全为0
			bool zero = true;
			for (int i = r + 2; i < n; i++)
			{
				if (A(i, r) != 0)
				{
					zero = false;
					break;
				}
			}

			if (!zero) // 如果A(r+1,r) ~ A(n,r)不全为0
			{
				// 计算d_r
				d_r = 0;
				for (int i = r + 1; i < n; i++)
				{
					d_r += A(i, r) * A(i, r);
				}
				d_r = sqrt(d_r);

				// 计算c_r
				c_r = -sgn(A(r + 1, r)) * d_r;
				// 计算h_r
				h_r = c_r * (c_r - A(r + 1, r));
				// 计算向量u_r
				for (int i = 0; i < n; i++)
				{
					if (i < r + 1)
						u_r[i] = 0;
					else
						u_r[i] = A(i, r);
				}
				u_r[r + 1] -= c_r;

				// 利用数组u_r构造列向量U_R
				U_R.row = A.row;
				U_R.column = 1;
				U_R.mat = u_r;

				// 进行一次Q、A的迭代
				//Q -= Q * U_R * U_R.trans() / h_r;
				P_R = A.trans() * U_R / h_r;
				Q_R = A*U_R / h_r;
				t_r = (P_R.trans()*U_R / h_r)(0, 0);
				A -= (Q_R - t_r*U_R)*U_R.trans() + U_R*P_R.trans();
				U_R.mat = nullptr;

				// 消除上三角阵R主对角线以下元素的截断误差
				for (int i = r + 2; i < n; i++)
				{
					A(i, r) = 0;
				}
			}
		}
		delete[] u_r;
		return A;
	}
	else
	{
		throw std::length_error("The matrix must be a square matrix.");
	}
	return Matrix();
}

// 计算二阶矩阵特征值,形参按行输入
std::vector<std::complex<double>> Matrix::eigVal22(double a, double b, double c, double d)const
{
	double num = (a + d)*(a + d) - 4 * (a*d - b*c);
	std::vector<std::complex<double>> vec;
	if (num < 0)
	{
		double real = (a + d) / 2;
		double imag = sqrt(-num) / 2;
		vec.emplace_back(real, imag);
		vec.emplace_back(real, -imag);
	}
	else
	{
		double realA = (a + d) / 2;
		double realB = sqrt(num) / 2;
		vec.emplace_back(realA + realB, 0);
		vec.emplace_back(realA - realB, 0);
	}
	return vec;
}

// QR方法计算特征值时的矩阵迭代
Matrix& Matrix::iterM(Matrix& A)const
{
	int m = A.row;
	if (m < 2)
	{
		throw std::length_error("The matrix dimension should be larger than 2.");
		return A;
	}
	double s = A(m - 2, m - 2) + A(m - 1, m - 1);
	double t = A(m - 2, m - 2) * A(m - 1, m - 1) - A(m - 2, m - 1) * A(m - 1, m - 2);
	Matrix M = A*A - s*A + t*(Matrix::eye(A.row, A.row));
	Matrix& B = M;

	double d_r = 0;
	double c_r = 0;
	double h_r = 0;
	double t_r = 0;
	double* u_r = new double[m];

	// 迭代所用向量
	Matrix U_R;
	Matrix P_R;
	Matrix V_R;
	Matrix Q_R;

	for (int r = 0; r < m - 1; r++)
	{
		// 判断A(r+1,r) ~ A(m,r)是否全为0
		bool zero = true;
		for (int i = r + 1; i < m; i++)
		{
			if (B(i, r) != 0)
			{
				zero = false;
				break;
			}
		}

		if (!zero) // 如果B(r+1,r) ~ A(m,r)不全为0
		{
			// 计算d_r
			d_r = 0;
			for (int i = r; i < m; i++)
			{
				d_r += B(i, r) * B(i, r);
			}
			d_r = sqrt(d_r);

			// 计算c_r
			c_r = -sgn(B(r, r)) * d_r;
			// 计算h_r
			h_r = c_r * (c_r - B(r, r));
			// 计算向量u_r
			for (int i = 0; i < m; i++)
			{
				if (i < r)
					u_r[i] = 0;
				else
					u_r[i] = B(i, r);
			}
			u_r[r] -= c_r;

			// 利用数组u_r构造列向量U_R
			U_R.row = m;
			U_R.column = 1;
			U_R.mat = u_r;

			// 进行一次Q、A的迭代
			V_R = B.trans() * U_R / h_r;
			B -= U_R * V_R.trans();
			P_R = A.trans() * U_R / h_r;
			Q_R = A * U_R / h_r;
			t_r = (P_R.trans() * U_R / h_r).mat[0];
			A -= (Q_R - t_r*U_R)*U_R.trans() + U_R * P_R.trans();
			U_R.mat = nullptr;
		}
	}
	delete[] u_r;
	return A;
}

//void Matrix::print()const
//{
//	int r = this->row;
//	int c = this->column;
//	for (int i = 0; i < r; i++)
//	{
//		for (int j = 0; j < c; j++)
//		{
//			std::cout << mat[i*c + j] << "\t";
//		}
//		std::cout << std::endl;
//	}
//	std::cout << std::endl;
//}
/***********************************  常用函数  ****************************************/
double sgn(const double& num)
{
	return num < 0 ? -1 : 1;
}

#endif