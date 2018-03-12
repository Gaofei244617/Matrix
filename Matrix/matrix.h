// ����C++11��׼
#ifndef _MATRIX_
#define _MATRIX_

//#include <iostream>

#include <string>
#include <vector>
#include <tuple>
#include <memory>
#include <cmath>
#include <complex>
/**********************************  ���ú���  *************************************************************************************************/
double sgn(const double& num);                                                 // ���ź���

																			   /***********************************  ������  ********************************************/
class Matrix
{
private:
	int row;                                                                   // ����
	int column;                                                                // ����
	double* mat;                                                               // ����һά����洢Ԫ��

public:
	Matrix();
	Matrix(const int& m, const int& n);                                        // ���湹��: m������n����
	Matrix(const std::initializer_list<double>& m);                            // �б���
	Matrix(const std::initializer_list<std::initializer_list<double>>& m);     // �б���
	Matrix(const Matrix& other);                                               // ��������
	Matrix(Matrix&& other);                                                    // �ƶ�����

																			   // ��ֵ����
	Matrix& operator=(const Matrix& other)noexcept;
	Matrix& operator=(Matrix&& other)noexcept;                                 // �ƶ���ֵ
	Matrix& operator+=(const double& num);                                     // M += num
	Matrix& operator+=(const Matrix& m);                                       // M += M
	Matrix& operator-=(const double& num);                                     // M -= num
	Matrix& operator-=(const Matrix& m);                                       // M -= M
	Matrix& operator*=(const double& num);                                     // M *= num
	Matrix& operator*=(const Matrix& m);                                       // M *= M
	Matrix& operator/=(const double& num);

	// �±�����
	inline double& operator()(int m, int n)noexcept;
	inline const double& operator()(int m, int n)const noexcept;
	inline double& at(int m, int n);
	inline const double& at(int m, int n)const;

	// �ӷ�����
	Matrix operator+()const;                                                   // +M
	friend Matrix operator+(const Matrix& m1, const Matrix& m2);               // M + M
	friend Matrix& operator+(const Matrix& m1, Matrix&& m2);                   // M + move(M)
	friend Matrix& operator+(Matrix&& m1, const Matrix& m2);                    // move(M) + M
	friend Matrix& operator+(Matrix&& m1, Matrix&& m2);                        // move(M) + move(M)
	friend Matrix operator+(const Matrix& m, const double& num);               // M + num
	friend Matrix operator+(const double& num, const Matrix& m);               // num + M
	friend Matrix& operator+(const double& num, Matrix&& m);                   // num + move(M)
	friend Matrix& operator+(Matrix&& m, const double& num);                   // move(M) + num

																			   // ��������
	Matrix operator-()const;                                                   // -M
	Matrix operator-(const Matrix& m)const;                                    // M - M
	Matrix& operator-(Matrix&& m)const;                                        // M - move(M)
	Matrix operator-(const double& num)const;                                  // M - num
	friend Matrix& operator-(Matrix&& m1, const Matrix& m2);                   // move(M) - M(or move(M))
	friend Matrix& operator-(Matrix&& m1, Matrix&& m2);                        // move(M) + move(M)
	friend Matrix operator-(const double& num, const Matrix& m);               // num - M
	friend Matrix& operator-(const double& num, Matrix&& m);                   // num - move(M)
	friend Matrix& operator-(Matrix&& m, const double& num);                   // move(M) - num

																			   // �˷�����
	Matrix operator*(const Matrix& m)const;                                    // M * M
	Matrix operator*(const double& num)const;                                  // M * num
	friend Matrix& operator*(Matrix&& m, const double& num);                   // move(M) * num
	friend Matrix& operator*(Matrix&& m1, Matrix&& m2);                        // move(M) * move(M)
	friend Matrix operator*(const double& num, const Matrix& m);               // num * M
	friend Matrix& operator*(const double& num, Matrix&& m);                   // num * move(M)
	Matrix dotMult(const Matrix& m)const;                                      // M .* M
	Matrix& dotMult(Matrix&& m)const;                                          // M .* move(M)

																			   // ��������
	Matrix operator/(const double& num);                                       // M / num
	Matrix dotDiv(const Matrix& m)const;                                       // M ./ M
	Matrix& dotDiv(Matrix&& m)const;                                           // M ./ move(M)

																			   // �߼������
	bool operator==(const Matrix& m)const noexcept;                            // ����
	bool operator!=(const Matrix& m)const noexcept;                            // ������

	const std::vector<int> size()const;                                        // ��ȡ����������
	std::vector<double> getRow(int n)const;                                    // ��ȡ������
	std::vector<double> getColumn(int n)const;                                 // ��ȡ������
	std::vector<double> getDiag()const;                                        // ��ȡ�Խ�Ԫ��
	double normOne()const;                                                     // 1����,�кͷ���
	double normInf()const;                                                     // �����,�кͷ���
	Matrix trans()const;                                                       // ����ת��

	const int rank()const;                                                     // �������
	const double trace()const;                                                 // ����ļ�
	Matrix inv()const;                                                         // �����
	const double det()const;                                                   // ��������ʽ
	std::vector<Matrix> QR()const;                                             // ����QR�ֽ�
	std::vector<std::complex<double>> eigs(double e = 0)const;                 // ��������ֵ
	std::vector<Matrix> SVD()const;                                            // ����ֵ�ֽ⣬����S��V��D��������
	Matrix subMat(int r1, int c1, int r2, int c2)const;                        // ����

	static Matrix eye(const int& m, const int& n);                             // ��λ����
	static Matrix ones(const int& m, const int& n);                            // Ԫ��ȫΪ1�ľ���
	static Matrix diag(const std::initializer_list<double>& nums);
	static Matrix rbind(const std::initializer_list<Matrix>& M);               // [M1; M2; ...], ��Ҫ�������
	static Matrix cbind(const std::initializer_list<Matrix>& M);               // [M1, M2, ...], ��Ҫ�������

																			   // ȫѡ��Ԫ��˹��ȥ��
	static std::tuple<Matrix, std::unique_ptr<int[]>, std::unique_ptr<int[]>, int, Matrix>
		gaussElimination(const Matrix& A, const Matrix& b = Matrix());

	static std::tuple<Matrix, int> solve(const Matrix& A, const Matrix& b);    // ������Է����� Ax = b

	virtual ~Matrix();                                                         // ����

private:
	Matrix hess2()const;                                                       // �����������Ƿֽ�(A = P*B*P'),PΪ��������,BΪ����������
	std::vector<std::complex<double>>
		eigVal22(double a, double b, double c, double d)const;                 // ������׾�������ֵ(��������)
	Matrix& iterM(Matrix& A)const;                                             // QR������������ֵ�еľ������
																			   //void print()const;
};

/**************************************************************************************************/
// ���湹��: m������n����
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

		// �жϸ���Ԫ�������Ƿ���ȣ����಻�� dimMatch = false
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

// ���ƹ���
Matrix::Matrix(const Matrix& other) : row(0), column(0), mat(nullptr)
{
	mat = new double[other.row * other.column];
	row = other.row;
	column = other.column;

	int size = row * column;
	for (int i = 0; i < size; i++)
		mat[i] = other.mat[i];
}

// �ƶ�����
Matrix::Matrix(Matrix&& other) : row(0), column(0), mat(nullptr)
{
	//*this = std::move(other);
	mat = other.mat;
	row = other.row;
	column = other.column;
	other.mat = nullptr;
}

// ���������(),���ڻ�ȡ����Ԫ�غ͸�ֵ,������0��ʼ
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

// ��ֵ����
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

// �ӷ�����
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

// ��������
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

// �˷�����
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

// ��������
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

// �߼�����
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

// ��ȡ����������
const std::vector<int> Matrix::size()const
{
	return std::vector<int>{row, column};
}

// ��ȡ������
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

// ��ȡ������
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

// ��ȡ�Խ�Ԫ��
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

// 1����(�кͷ���)
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

// �����(�кͷ���)
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

// ����ת��
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

// �������
const int Matrix::rank()const { return std::get<3>(gaussElimination(*this)); }

// ����ļ�
const double Matrix::trace()const
{
	if (this->row == this->column)
	{
		if (this->row == 0)
		{
			return 0; // �վ���ļ�Ϊ0
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

// ����ʽ
const double Matrix::det()const
{
	// ������Ϊ�����Ǿ��󣬶Խ�Ԫ�س˻���Ϊ����ʽ��ֵ
	// ���������н�������ȷ��
	if (this->row == this->column)
	{
		if (this->row == 0)       // �վ�������ʽΪ 1
			return 1;
		if (this->row == 1)
			return this->mat[0];

		int flag = 1;
		double temp = 0;
		int col = 0;
		Matrix matrix(*this);
		for (int i = 0; i < matrix.row - 1; i++)
		{
			// ѡ����Ԫ
			temp = abs(matrix.mat[i*column + i]);
			col = i;
			for (int j = i + 1; j < matrix.row; j++)
			{
				if (abs(matrix.mat[j*column + i]) > temp)
				{
					temp = abs(matrix.mat[j*column + i]);
					col = j; // ��¼��Ԫ������
				}
			}
			// �н���
			if (col != i)
			{
				for (int k = 0; k < matrix.column; k++)
				{
					temp = matrix.mat[i*column + k];
					matrix.mat[i*column + k] = matrix.mat[col*column + k];
					matrix.mat[col*column + k] = temp;
				}
				flag *= -1;
				// ��ԪΪ0��˵���������ȣ�����ʽΪ0
				if (matrix.mat[i*column + i] == 0)
					return 0;
			}
			// ��˹��Ԫ
			for (int j = i + 1; j < matrix.row; j++)
			{
				temp = matrix.mat[j*column + i] / matrix.mat[i*column + i];
				for (int k = i; k < matrix.column; k++)
				{
					matrix.mat[j*column + k] -= temp * matrix.mat[i*column + k];
				}
			}
		}
		// �Խ�Ԫ�صĳ˻�
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

/* ����QR�ֽ�(A = QR),QΪ��������,RΪ�����Ǿ���(�����Խ���Ԫ��>=0ʱ,
�ֽ����Ψһ��)����R���ȣ����Խ���Ԫ�� > 0 */
std::vector<Matrix> Matrix::QR()const
{
	if (this->row > 0 && this->row == this->column)
	{
		Matrix A(*this);
		Matrix Q = Matrix::eye(this->row, this->row);

		int n = A.row;   // ����ά��
		double d_r = 0;
		double c_r = 0;
		double h_r = 0;
		double* u_r = new double[n];
		Matrix U_R;
		Matrix P_R;

		for (int r = 0; r < n - 1; r++)
		{
			// �ж�A(r+1,r) ~ A(n,r)�Ƿ�ȫΪ0
			bool zero = true;
			for (int i = r + 1; i < n; i++)
			{
				if (A(i, r) != 0)
				{
					zero = false;
					break;
				}
			}

			if (!zero) // ���A(r+1,r) ~ A(n,r)��ȫΪ0
			{
				// ����d_r
				d_r = 0;
				for (int i = r; i < n; i++)
				{
					d_r += A(i, r) * A(i, r);
				}
				d_r = sqrt(d_r);

				// ����c_r
				c_r = -sgn(A(r, r)) * d_r;
				// ����h_r
				h_r = c_r * (c_r - A(r, r));
				// ��������u_r
				for (int i = 0; i < n; i++)
				{
					if (i < r)
						u_r[i] = 0;
					else
						u_r[i] = A(i, r);
				}
				u_r[r] -= c_r;

				// ��������u_r����������U_R
				U_R.row = A.row;
				U_R.column = 1;
				U_R.mat = u_r;

				// ����һ��Q��A�ĵ���
				Q -= Q * U_R * U_R.trans() / h_r;
				P_R = A.trans() * U_R / h_r;
				A -= U_R * P_R.trans();
				U_R.mat = nullptr;

				// ������������R���Խ�������Ԫ�صĽض����
				for (int i = r + 1; i < n; i++)
				{
					A(i, r) = 0;
				}
			}
		}
		delete[] u_r;

		// ������ɺ� Q = Q_r, R = A_r
		return std::vector<Matrix>({ Q, A });
	}
	else
	{
		throw std::length_error("The matrix must be a square matrix.");
	}
	return std::vector<Matrix>();
}

// �����ȫ������ֵ(��˫��λ�Ƶ�QR����), eΪ����ˮƽ
std::vector<std::complex<double>> Matrix::eigs(double e)const
{
	// ֻ�зǿշ������������ֵ
	if (this->row != this->column || this->row == 0)
	{
		if (this->row != this->column)
		{
			throw std::length_error("In Matrix::eigs: The matrix must be a square matrix.");
		}
		return std::vector<std::complex<double>>();
	}

	std::vector<std::complex<double>> vec;                     // �洢����ֵ������
	Matrix matrix(this->hess2());                              // ������������ǻ�
	e == 0 ? e = this->normOne() / 1.0e10 : e = abs(e);        // ȷ������ˮƽ,Ĭ��Ϊ����1������1/(1e10)
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
				// ����A_k+1
				// iterM��QR������������ֵʱ�ľ������
				matrix = iterM(matrix.subMat(0, 0, m - 1, m - 1));
			}
		}
	}
	return vec;
}

// ��������
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

			// ��Ѱ����������Ԫ����¼��λ��
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

			// ȫѡ��Ԫ�����к���
			// �н���
			if (cFlag[k] != k)
			{
				for (int i = 0; i < n; i++)
				{
					temp = matrix.mat[i * column + k];
					matrix.mat[i * column + k] = matrix.mat[i * column + cFlag[k]];
					matrix.mat[i * column + cFlag[k]] = temp;
				}
			}

			// �н���
			if (rFlag[k] != k)
			{
				for (int i = 0; i < n; i++)
				{
					temp = matrix.mat[k * column + i];
					matrix.mat[k * column + i] = matrix.mat[rFlag[k] * column + i];
					matrix.mat[rFlag[k] * column + i] = temp;
				}
			}

			// ��ԪΪ0, ��������, û�������
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

		// ��ԭ����: (1) �Ƚ�������(��)����лָ�; (2) ԭ������(��)��������(��)���ָ�
		for (int k = n - 1; k >= 0; k--)
		{
			// �н���
			if (rFlag[k] != k)
			{
				for (int i = 0; i < n; i++)
				{
					temp = matrix.mat[i * column + k];
					matrix.mat[i * column + k] = matrix.mat[i * column + rFlag[k]];
					matrix.mat[i * column + rFlag[k]] = temp;
				}
			}

			// �н���
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

// ����,�β�Ϊ���Ͻ�Ԫ��λ�ú����½�Ԫ��λ��
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

// ��̬����
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

		// �жϾ����Ƿ��ƴ��,������ƴ�Ӻ���������
		for (int i = 1; i < count; i++)
		{
			row += (*(M.begin() + i)).row;

			if ((*(M.begin() + i)).column != column)
			{
				throw std::length_error("Dimensions do not match.");
				return Matrix();
			}
		}

		// ƴ�Ӿ��� [M1; M2; ...]
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

		// �жϾ����Ƿ��ƴ��,������ƴ�Ӻ���������
		for (int i = 1; i < count; i++)
		{
			column += (*(M.begin() + i)).column;

			if ((*(M.begin() + i)).row != row)
			{
				throw std::length_error("Dimensions do not match.");
				return Matrix();
			}
		}

		// ƴ�Ӿ��� [M1, M2, ...]
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
// ȫѡ��Ԫ��˹��ȥ
// ����ֵ�ֱ�Ϊ��{��Ԫ��ľ����н�����¼���н�����¼���������}
std::tuple<Matrix, std::unique_ptr<int[]>, std::unique_ptr<int[]>, int, Matrix>
Matrix::gaussElimination(const Matrix& A, const Matrix& b)
{
	/***************************************************************/
	// �β�vec[��ѡ]��Ĭ��Ϊ������
	// ��vec��Ϊ�գ����贫�����Է����� Ax = b�е�����b��
	// ��������ֵ�еĵ��ĸ�����Ϊ����bͨ����˹��Ԫ�����е��б任�õ���������
	/***************************************************************/

	// �վ�����Ϊ0
	if (A.row == 0)
		return std::make_tuple(Matrix(), nullptr, nullptr, 0, Matrix());

	// �β�vecΪһ����������Ԫ���������ͬ������
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

	// ��¼��Ԫλ��
	std::unique_ptr<int[]> R(new int[size]());          // �н�����Ϣ
	std::unique_ptr<int[]> C(new int[size]());          // �н�����Ϣ

	double temp;
	// ȫѡ��Ԫ��Ԫ
	for (int k = 0; k < size; k++)
	{
		temp = abs(matrix.mat[k * column + k]);
		R[k] = C[k] = k;

		// ��Ѱ����������Ԫ����¼��λ��
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

		// ȫѡ��Ԫ�����к���
		// �н���
		if (C[k] != k)
		{
			for (int i = 0; i < matrix.row; i++)
			{
				temp = matrix.mat[i * column + k];
				matrix.mat[i * column + k] = matrix.mat[i * column + C[k]];
				matrix.mat[i * column + C[k]] = temp;
			}
		}

		// �н���
		if (R[k] != k)
		{
			// ����matrix�����н���
			for (int i = 0; i < matrix.column; i++)
			{
				temp = matrix.mat[k * column + i];
				matrix.mat[k * column + i] = matrix.mat[R[k] * column + i];
				matrix.mat[R[k] * column + i] = temp;
			}

			// ����vec�����н���
			if (vec != Matrix())
			{
				temp = vec.mat[k];
				vec.mat[k] = vec.mat[R[k]];
				vec.mat[R[k]] = temp;
			}
		}

		// ��ԪΪ0����Ԫ���½�Ԫ�ؾ�Ϊ0, ��˹��Ԫ���
		if (matrix.mat[k*column + k] == 0)
			break;

		// ��Ԫ
		for (int i = k + 1; i < matrix.row; i++)
		{
			// ����matrix��Ԫ
			temp = matrix.mat[i*column + k] / matrix.mat[k*column + k];
			matrix.mat[i*column + k] = 0;
			for (int j = k + 1; j < matrix.column; j++)
			{
				matrix.mat[i*column + j] -= temp * matrix.mat[k*column + j];
			}

			// ����vec��Ԫ
			if (vec != Matrix())
			{
				vec.mat[i] -= temp * vec.mat[k];
			}
		}
	}

	// ���Խ��߷���Ԫ�ظ���Ϊ�������
	int count = 0;
	for (int i = 0; i < size; i++)
	{
		if (matrix.mat[i*column + i] != 0)
			count++;
	}

	return std::make_tuple(matrix, std::move(R), std::move(C), count, vec);
}

// ������Է�����: Ax = b
std::tuple<Matrix, int> Matrix::solve(const Matrix& A, const Matrix& b)
{
	/*******************************************************/
	// ������Է�����: Ax = b
	// AΪϵ������bΪĿ������
	// ����ȫѡ��Ԫ�ĸ�˹��Ԫ����ⷽ����

	// (1)��������Ψһ�⣺   ���� tuple(������x,  ��־��0)
	// (2)�������������⣺ ���� tuple(һ���ؽ�x, ��־��1)
	// (3)�������޾�ȷ�⣺   ���� tuple(���ƽ�x*,  ��־��2)
	//    ���ƽ�x*ʹ��|Ax* - b|��С
	/********************************************************/

	// �ж�ϵ�������������Ŀ�������������Ƿ����
	if (A.row != b.row)
	{
		throw std::length_error("Size of coefficient matrix and target vector does not match.");
	}

	auto temp = gaussElimination(A, b);                // ��ϵ���������ȫѡ��Ԫ�ĸ�˹��Ԫ
	Matrix& matrix = std::get<0>(temp);                // ȫѡ��Ԫ��˹��Ԫ��ľ���
	std::unique_ptr<int[]>& R = std::get<1>(temp);     // �н�����Ϣ
	std::unique_ptr<int[]>& C = std::get<2>(temp);     // �н�����Ϣ
	int rankA = std::get<3>(temp);                     // ϵ���������rank(A)
	int rankAb = cbind({ A,b }).rank();                // rank(A,b)
	int n = A.column;                                  // δ֪������
	Matrix& vec = std::get<4>(temp);                   // ��˹��Ԫ�󷽳����Ŀ������

													   // rank(A) == rank(A,b) && rank(A) == n: ������Ψһ��
													   // rank(A) == rank(A,b) && rank(A) < n:  ������������
													   // others: �����޽�( rank(A)�����ܴ���n )

													   // (1) �����н�
	if (rankA == rankAb)
	{
		Matrix x(n, 1);   // ���̵Ľ�

						  // ����˹��Ԫ��������Ƿ�����
		for (int i = rankA - 1; i >= 0; i--)
		{
			for (int j = 0; j < rankA - i - 1; j++)
			{
				vec(i, 0) -= matrix(i, rankA - 1 - j) * x(rankA - 1 - j, 0);
			}

			x(i, 0) = 1.0 * vec(i, 0) / matrix(i, i);
		}

		// ��������������ʱ���������x����Ԫ��Ϊ0���õ�һ���ؽ�
		if (rankA < n)
		{
			for (int i = rankA; i < n; i++)
			{
				x(i, 0) = 0;
			}
		}

		// �����н�����Ϣ,�ָ�����x��Ԫ�ص�˳��
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

		// (1.1) ������Ψһ��
		if (rankA == n)
		{
			return std::make_tuple(x, 0);
		}

		// (1.2) ��������������
		if (rankA < n)
		{
			return std::make_tuple(x, 1);
		}
	}

	// (2)�����޽�
	if (rankA != rankAb)
	{
		Matrix AT(A.trans());
		return std::make_tuple(std::get<0>(solve(AT*A, AT*b)), 2);
	}

	return std::make_tuple(Matrix(), 2);
}

// ����
Matrix::~Matrix() { delete[] mat; }

// �����������ǻ�
Matrix Matrix::hess2()const
{
	if (this->row > 0 && this->row == this->column)
	{
		Matrix A(*this);

		int n = A.row;   // ����ά��
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
			// �ж�A(r+2,r) ~ A(n,r)�Ƿ�ȫΪ0
			bool zero = true;
			for (int i = r + 2; i < n; i++)
			{
				if (A(i, r) != 0)
				{
					zero = false;
					break;
				}
			}

			if (!zero) // ���A(r+1,r) ~ A(n,r)��ȫΪ0
			{
				// ����d_r
				d_r = 0;
				for (int i = r + 1; i < n; i++)
				{
					d_r += A(i, r) * A(i, r);
				}
				d_r = sqrt(d_r);

				// ����c_r
				c_r = -sgn(A(r + 1, r)) * d_r;
				// ����h_r
				h_r = c_r * (c_r - A(r + 1, r));
				// ��������u_r
				for (int i = 0; i < n; i++)
				{
					if (i < r + 1)
						u_r[i] = 0;
					else
						u_r[i] = A(i, r);
				}
				u_r[r + 1] -= c_r;

				// ��������u_r����������U_R
				U_R.row = A.row;
				U_R.column = 1;
				U_R.mat = u_r;

				// ����һ��Q��A�ĵ���
				//Q -= Q * U_R * U_R.trans() / h_r;
				P_R = A.trans() * U_R / h_r;
				Q_R = A*U_R / h_r;
				t_r = (P_R.trans()*U_R / h_r)(0, 0);
				A -= (Q_R - t_r*U_R)*U_R.trans() + U_R*P_R.trans();
				U_R.mat = nullptr;

				// ������������R���Խ�������Ԫ�صĽض����
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

// ������׾�������ֵ,�βΰ�������
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

// QR������������ֵʱ�ľ������
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

	// ������������
	Matrix U_R;
	Matrix P_R;
	Matrix V_R;
	Matrix Q_R;

	for (int r = 0; r < m - 1; r++)
	{
		// �ж�A(r+1,r) ~ A(m,r)�Ƿ�ȫΪ0
		bool zero = true;
		for (int i = r + 1; i < m; i++)
		{
			if (B(i, r) != 0)
			{
				zero = false;
				break;
			}
		}

		if (!zero) // ���B(r+1,r) ~ A(m,r)��ȫΪ0
		{
			// ����d_r
			d_r = 0;
			for (int i = r; i < m; i++)
			{
				d_r += B(i, r) * B(i, r);
			}
			d_r = sqrt(d_r);

			// ����c_r
			c_r = -sgn(B(r, r)) * d_r;
			// ����h_r
			h_r = c_r * (c_r - B(r, r));
			// ��������u_r
			for (int i = 0; i < m; i++)
			{
				if (i < r)
					u_r[i] = 0;
				else
					u_r[i] = B(i, r);
			}
			u_r[r] -= c_r;

			// ��������u_r����������U_R
			U_R.row = m;
			U_R.column = 1;
			U_R.mat = u_r;

			// ����һ��Q��A�ĵ���
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
/***********************************  ���ú���  ****************************************/
double sgn(const double& num)
{
	return num < 0 ? -1 : 1;
}

#endif