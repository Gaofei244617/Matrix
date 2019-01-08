#include "mat_def.h"
#include "matrix.h"

namespace mat
{
    // 常规构造: m行数，n列数
    Matrix::Matrix() : row(0), column(0), mat_data(nullptr) {}
    Matrix::Matrix(const std::size_t& m, const std::size_t& n) : row(0), column(0), mat_data(nullptr)
    {
        if (m > 0 && n > 0)
        {
            mat_data = new double[m*n]();
            row = m;
            column = n;
        }
    }
    Matrix::Matrix(const std::initializer_list<double>& m) : row(0), column(0), mat_data(nullptr)
    {
        std::size_t length = m.size();
        if (length > 0)
        {
            mat_data = new double[length];
            this->row = 1;
            this->column = length;

            for (size_t i = 0; i < length; i++)
            {
                mat_data[i] = *(m.begin() + i);
            }
        }
    }
    Matrix::Matrix(const std::initializer_list<std::initializer_list<double>>& m) :row(0), column(0), mat_data(nullptr)
    {
        std::size_t row = m.size();
        if (row > 0)
        {
            bool dimMatch = true;
            std::size_t column = (*(m.begin())).size();

            // 判断各行元素数量是否相等，若相不等 dimMatch = false
            for (std::size_t i = 0; i < row; i++)
            {
                if (column != (*(m.begin() + i)).size())
                {
                    dimMatch = false;
                }
            }

            if (dimMatch)
            {
                mat_data = new double[row * column];
                this->row = row;
                this->column = column;
                for (std::size_t i = 0; i < row; i++)
                {
                    for (std::size_t j = 0; j < column; j++)
                    {
                        mat_data[i*column + j] = *((*(m.begin() + i)).begin() + j);
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
    Matrix::Matrix(const Matrix& other) : row(0), column(0), mat_data(nullptr)
    {
        mat_data = new double[other.row * other.column];
        row = other.row;
        column = other.column;

        std::size_t size = row * column;
        for (std::size_t i = 0; i < size; i++)
        {
            mat_data[i] = other.mat_data[i];
        }
    }

    // 移动构造
    Matrix::Matrix(Matrix&& other) : row(0), column(0), mat_data(nullptr)
    {
        //*this = std::move(other);
        mat_data = other.mat_data;
        row = other.row;
        column = other.column;
        other.mat_data = nullptr;
    }

    // 析构函数
    Matrix::~Matrix()
    {
        delete[] mat_data;
    }

    // 元素索引
    double& Matrix::at(std::size_t m, std::size_t n)
    {
        if (m < row && n < column)
        {
            return mat_data[column*m + n];
        }
        else
        {
            throw std::out_of_range("In function at():Out of index range of the matrix.");
        }
    }
    const double& Matrix::at(std::size_t m, std::size_t n)const
    {
        if (m < row && n < column)
        {
            return mat_data[column*m + n];
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
            delete[] this->mat_data;
            std::size_t size = other.row * other.column;
            this->mat_data = new double[size];
            this->row = other.row;
            this->column = other.column;

            for (std::size_t i = 0; i < size; i++)
            {
                this->mat_data[i] = other.mat_data[i];
            }
        }
        return *this;
    }
    Matrix& Matrix::operator=(Matrix&& other)noexcept
    {
        if (this != &other)
        {
            delete[] this->mat_data;
            this->mat_data = other.mat_data;
            this->row = other.row;
            this->column = other.column;
            other.mat_data = nullptr;
        }
        return *this;
    }
    Matrix& Matrix::operator+=(const Matrix& m)
    {
        if (row == m.row && column == m.column)
        {
            std::size_t size = row * column;
            for (std::size_t i = 0; i < size; i++)
            {
                mat_data[i] += m.mat_data[i];
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
            std::size_t size = row * column;
            for (std::size_t i = 0; i < size; i++)
            {
                mat_data[i] += num;
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
            std::size_t size = row * column;
            for (std::size_t i = 0; i < size; i++)
            {
                mat_data[i] -= num;
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
            std::size_t size = row * column;
            for (std::size_t i = 0; i < size; i++)
            {
                mat_data[i] -= m.mat_data[i];
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
            std::size_t size = row * column;
            for (std::size_t i = 0; i < size; i++)
            {
                mat_data[i] *= num;
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
        if (std::get<1>(this->size()) == std::get<0>(m.size()))
        {
            std::size_t row = std::get<0>(this->size());
            std::size_t column = std::get<1>(m.size());
            Matrix matrix(row, column);

            std::size_t count = std::get<0>(m.size());
            for (std::size_t i = 0; i < row; i++)
            {
                for (std::size_t j = 0; j < column; j++)
                {
                    double sum = 0;
                    for (std::size_t k = 0; k < count; k++)
                    {
                        sum += (*this)[i][k] * m[k][j];
                    }
                    matrix[i][j] = sum;
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
                std::size_t size = row * column;
                for (std::size_t i = 0; i < size; i++)
                {
                    mat_data[i] /= num;
                }
            }
            else
            {
                throw std::length_error("Null Matrix.");
            }
        }
        return *this;
    }

    // 点乘(除)运算
    Matrix Matrix::dotMult(const Matrix& m)const
    {
        if (this->row == m.row && this->column == m.column)
        {
            Matrix matrix(this->row, this->column);
            std::size_t size = this->row * this->column;
            for (std::size_t i = 0; i < size; i++)
            {
                matrix.mat_data[i] = this->mat_data[i] * m.mat_data[i];
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
            std::size_t size = this->row * this->column;
            for (std::size_t i = 0; i < size; i++)
            {
                m.mat_data[i] *= this->mat_data[i];
            }
            return m;
        }
        else
        {
            throw std::length_error("Dimensions do not match.");
        }
        return m = *this;
    }
    Matrix Matrix::dotDiv(const Matrix& m)const
    {
        if (this->row == m.row && this->column == m.column)
        {
            bool flag = true;
            std::size_t size = this->row * this->column;
            for (std::size_t i = 0; i < size; i++)
            {
                if (m.mat_data[i] == 0)
                {
                    flag = false;
                }
            }

            if (flag)
            {
                Matrix matrix(this->row, this->column);
                for (std::size_t i = 0; i < size; i++)
                {
                    matrix.mat_data[i] = this->mat_data[i] / m.mat_data[i];
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
                if (m.mat_data[i] == 0)
                {
                    flag = false;
                }
            }

            if (flag)
            {
                for (int i = 0; i < size; i++)
                {
                    m.mat_data[i] = this->mat_data[i] / m.mat_data[i];
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

    // 获取矩阵行列数
    const std::pair<std::size_t, std::size_t> Matrix::size()const
    {
        return std::make_pair(row, column);
    }

    // 获取行向量
    Matrix Matrix::getRow(std::size_t n)const
    {
        if (this->row > 0)
        {
            Matrix vec(1, this->column);
            for (std::size_t i = 0; i < this->column; i++)
            {
                vec[0][i] = mat_data[n*column + i];
            }
            return vec;
        }
        return Matrix();
    }

    // 获取列向量
    Matrix Matrix::getColumn(std::size_t n)const
    {
        if (this->row > 0)
        {
            Matrix vec(this->row, 1);
            for (std::size_t i = 0; i < this->row; i++)
            {
                vec[1][0] = mat_data[i*column + n];
            }
            return vec;
        }
        return Matrix();
    }

    // 获取对角元素
    std::vector<double> Matrix::getDiag()const
    {
        if (this->row == 0)
        {
            return std::vector<double>();
        }

        std::size_t size;
        this->row < this->column ? size = this->row : size = this->column;
        std::vector<double> vec;
        for (std::size_t i = 0; i < size; i++)
        {
            vec.emplace_back(mat_data[i*column + i]);
        }
        return vec;
    }

    // 1范数(列和范数)
    double Matrix::normOne()const
    {
        double norm = 0;
        if (this->row > 0)
        {
            for (std::size_t i = 0; i < this->column; i++)
            {
                double temp = 0;
                for (std::size_t j = 0; j < this->row; j++)
                {
                    temp += abs(this->mat_data[j*column + i]);
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
            for (std::size_t i = 0; i < this->row; i++)
            {
                double temp = 0;
                for (std::size_t j = 0; j < this->column; j++)
                {
                    temp += abs(this->mat_data[i*column + j]);
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
            for (std::size_t i = 0; i < column; i++)
            {
                for (std::size_t j = 0; j < row; j++)
                {
                    matrix.mat_data[i*matrix.column + j] = this->mat_data[j*column + i];
                }
            }
            return matrix;
        }
        return Matrix();
    }

    // 矩阵的秩
    std::size_t Matrix::rank()const
    {
        return std::get<3>(gaussElimination(*this));
    }

    // 矩阵的迹
    double Matrix::trace()const
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
                for (std::size_t i = 0; i < this->row; i++)
                {
                    sum += this->mat_data[i*column + i];
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
    double Matrix::det()const
    {
        // 将矩阵化为上三角矩阵，对角元素乘积即为行列式的值
        // 正负号由行交换次数确定
        if (this->row == this->column)
        {
            if (this->row == 0)       // 空矩阵行列式为 1
                return 1;
            if (this->row == 1)
                return this->mat_data[0];

            int flag = 1;
            double temp = 0;
            int col = 0;
            Matrix matrix(*this);
            for (std::size_t i = 0; i < matrix.row - 1; i++)
            {
                // 选列主元
                temp = abs(matrix.mat_data[i*column + i]);
                col = i;
                for (std::size_t j = i + 1; j < matrix.row; j++)
                {
                    if (abs(matrix.mat_data[j*column + i]) > temp)
                    {
                        temp = abs(matrix.mat_data[j*column + i]);
                        col = j; // 记录主元所在行
                    }
                }
                // 行交换
                if (col != i)
                {
                    for (std::size_t k = 0; k < matrix.column; k++)
                    {
                        temp = matrix.mat_data[i*column + k];
                        matrix.mat_data[i*column + k] = matrix.mat_data[col*column + k];
                        matrix.mat_data[col*column + k] = temp;
                    }
                    flag *= -1;
                    // 主元为0，说明方阵不满秩，行列式为0
                    if (matrix.mat_data[i*column + i] == 0)
                        return 0;
                }
                // 高斯消元
                for (std::size_t j = i + 1; j < matrix.row; j++)
                {
                    temp = matrix.mat_data[j*column + i] / matrix.mat_data[i*column + i];
                    for (std::size_t k = i; k < matrix.column; k++)
                    {
                        matrix.mat_data[j*column + k] -= temp * matrix.mat_data[i*column + k];
                    }
                }
            }
            // 对角元素的乘积
            temp = 1;
            for (std::size_t i = 0; i < matrix.row; i++)
            {
                temp *= matrix.mat_data[i*column + i];
            }
            return temp * flag;
        }
        else
        {
            throw std::length_error("The matrix must be a square matrix.");
        }
        return 0;
    }

    /* 矩阵QR分解(A = QR),Q为正交矩阵,R为上三角矩阵(且主对角线元素>=0时,
    分解具有唯一性)。若R满秩，主对角线元素 > 0 */
    std::pair<Matrix, Matrix> Matrix::QR()const
    {
        if (this->row > 0 && this->row == this->column)
        {
            Matrix A(*this);
            Matrix Q = eye(this->row, this->row);

            std::size_t n = A.row;   // 矩阵维数
            double d_r = 0;
            double c_r = 0;
            double h_r = 0;
            double* u_r = new double[n];
            Matrix U_R;
            Matrix P_R;

            for (std::size_t r = 0; r < n - 1; r++)
            {
                // 判断A(r+1,r) ~ A(n,r)是否全为0
                bool zero = true;
                for (std::size_t i = r + 1; i < n; i++)
                {
                    if (A[i][r] != 0)
                    {
                        zero = false;
                        break;
                    }
                }

                if (!zero) // 如果A(r+1,r) ~ A(n,r)不全为0
                {
                    // 计算d_r
                    d_r = 0;
                    for (std::size_t i = r; i < n; i++)
                    {
                        d_r += A[i][r] * A[i][r];
                    }
                    d_r = sqrt(d_r);

                    // 计算c_r
                    c_r = -sgn(A[r][r]) * d_r;
                    // 计算h_r
                    h_r = c_r * (c_r - A[r][r]);
                    // 计算向量u_r
                    for (std::size_t i = 0; i < n; i++)
                    {
                        if (i < r)
                            u_r[i] = 0;
                        else
                            u_r[i] = A[i][r];
                    }
                    u_r[r] -= c_r;

                    // 利用数组u_r构造列向量U_R
                    U_R.row = A.row;
                    U_R.column = 1;
                    U_R.mat_data = u_r;

                    // 进行一次Q、A的迭代
                    Q -= Q * U_R * U_R.trans() / h_r;
                    P_R = A.trans() * U_R / h_r;
                    A -= U_R * P_R.trans();
                    U_R.mat_data = nullptr;

                    // 消除上三角阵R主对角线以下元素的截断误差
                    for (std::size_t i = r + 1; i < n; i++)
                    {
                        A[i][r] = 0;
                    }
                }
            }
            delete[] u_r;

            // 迭代完成后 Q = Q_r, R = A_r
            return std::make_pair(Q, A);
        }
        else
        {
            throw std::length_error("The matrix must be a square matrix.");
        }
        return std::make_pair(Matrix(), Matrix());
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
                vec.emplace_back(matrix[0][0], 0);
                break;
            }

            if (abs(matrix[m - 1][m - 2]) < e)
            {
                vec.emplace_back(matrix[m - 1][m - 1], 0);
                m = m - 1;
            }
            else
            {
                if (m == 2 || abs(matrix[m - 2][m - 3]) < e)
                {
                    auto val = eigVal22(matrix[m - 2][m - 2], matrix[m - 2][m - 1], matrix[m - 1][m - 2], matrix[m - 1][m - 1]);
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
                    Matrix temp = matrix.subMat(0, 0, m - 1, m - 1);
                    matrix = iterM(temp);
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
            std::size_t n = matrix.row;
            std::unique_ptr<std::size_t[]> rFlag(new std::size_t[n]);
            std::unique_ptr<std::size_t[]> cFlag(new std::size_t[n]);

            double temp;
            for (std::size_t k = 0; k < n; k++)
            {
                temp = abs(matrix.mat_data[k * column + k]);
                rFlag[k] = cFlag[k] = k;

                // 搜寻子阵的最大主元并记录其位置
                for (std::size_t i = k; i < n; i++)
                {
                    for (std::size_t j = k; j < n; j++)
                    {
                        if (abs(matrix.mat_data[i * column + j]) > temp)
                        {
                            temp = abs(matrix.mat_data[i * column + j]);
                            rFlag[k] = i;
                            cFlag[k] = j;
                        }
                    }
                }

                // 全选主元交换行和列
                // 列交换
                if (cFlag[k] != k)
                {
                    for (std::size_t i = 0; i < n; i++)
                    {
                        temp = matrix.mat_data[i * column + k];
                        matrix.mat_data[i * column + k] = matrix.mat_data[i * column + cFlag[k]];
                        matrix.mat_data[i * column + cFlag[k]] = temp;
                    }
                }

                // 行交换
                if (rFlag[k] != k)
                {
                    for (std::size_t i = 0; i < n; i++)
                    {
                        temp = matrix.mat_data[k * column + i];
                        matrix.mat_data[k * column + i] = matrix.mat_data[rFlag[k] * column + i];
                        matrix.mat_data[rFlag[k] * column + i] = temp;
                    }
                }

                // 主元为0, 矩阵不满秩, 没有逆矩阵
                if (matrix.mat_data[k*column + k] == 0)
                {
                    throw std::logic_error("No inverse matrix.");
                    return Matrix();
                }
                else
                {
                    matrix.mat_data[k*column + k] = 1 / matrix.mat_data[k*column + k];
                }

                for (std::size_t j = 0; j < n; j++)
                {
                    if (j != k)
                    {
                        matrix.mat_data[k * column + j] *= matrix.mat_data[k * column + k];
                    }
                }

                for (std::size_t i = 0; i < n; i++)
                {
                    for (std::size_t j = 0; j < n; j++)
                    {
                        if (i != k && j != k)
                        {
                            matrix.mat_data[i * column + j] -= matrix.mat_data[i * column + k] * matrix.mat_data[k * column + j];
                        }
                    }
                }

                for (std::size_t i = 0; i < n; i++)
                {
                    if (i != k)
                    {
                        matrix.mat_data[i * column + k] = -matrix.mat_data[i * column + k] * matrix.mat_data[k * column + k];
                    }
                }
            }

            // 还原矩阵: (1) 先交换的行(列)后进行恢复; (2) 原来的行(列)交换用列(行)来恢复
            for (int k = n - 1; k >= 0; k--)
            {
                // 列交换
                if (rFlag[k] != k)
                {
                    for (std::size_t i = 0; i < n; i++)
                    {
                        temp = matrix.mat_data[i * column + k];
                        matrix.mat_data[i * column + k] = matrix.mat_data[i * column + rFlag[k]];
                        matrix.mat_data[i * column + rFlag[k]] = temp;
                    }
                }

                // 行交换
                if (cFlag[k] != k)
                {
                    for (std::size_t i = 0; i < n; i++)
                    {
                        temp = matrix.mat_data[k * column + i];
                        matrix.mat_data[k * column + i] = matrix.mat_data[cFlag[k] * column + i];
                        matrix.mat_data[cFlag[k] * column + i] = temp;
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
    Matrix Matrix::subMat(std::size_t r1, std::size_t c1, std::size_t r2, std::size_t c2)const
    {
        if (r2 < this->row && c2 < this->column && r1 <= r2 && c1 <= c2)
        {
            Matrix matrix(r2 - r1 + 1, c2 - c1 + 1);
            for (std::size_t i = r1; i <= r2; i++)
            {
                for (std::size_t j = c1; j <= c2; j++)
                {
                    matrix.mat_data[(i - r1)*matrix.column + j - c1] = this->mat_data[i*column + j];
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

    // 全选主元高斯消去
    // 返回值分别为：{消元后的矩阵，行交换记录，列交换记录，矩阵的秩，消元后的向量}
    std::tuple<Matrix, std::unique_ptr<std::size_t[]>, std::unique_ptr<std::size_t[]>, std::size_t, Matrix>
        Matrix::gaussElimination(const Matrix& A, const Matrix& b)
    {
        /***************************************************************/
        // 形参vec[可选]，默认为空向量
        // 若vec不为空，则需传入线性方程组 Ax = b中的向量b，
        // 函数返回值中的第四个参数为向量b通过高斯消元过程中的行变换得到的新向量
        /***************************************************************/

        // 空矩阵秩为0
        if (A.row == 0)
        {
            return std::make_tuple(Matrix(), nullptr, nullptr, 0, Matrix());
        }

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
        std::size_t size = matrix.row < matrix.column ? matrix.row : matrix.column;
        std::size_t row = matrix.row;
        std::size_t column = matrix.column;

        // 记录主元位置
        std::unique_ptr<std::size_t[]> R(new std::size_t[size]());          // 行交换信息
        std::unique_ptr<std::size_t[]> C(new std::size_t[size]());          // 列交换信息

        double temp;
        // 全选主元消元
        for (std::size_t k = 0; k < size; k++)
        {
            temp = abs(matrix.mat_data[k * column + k]);
            R[k] = C[k] = k;

            // 搜寻子阵的最大主元并记录其位置
            for (std::size_t i = k; i < matrix.row; i++)
            {
                for (std::size_t j = k; j < matrix.column; j++)
                {
                    if (abs(matrix.mat_data[i * column + j]) > temp)
                    {
                        temp = abs(matrix.mat_data[i * column + j]);
                        R[k] = i;
                        C[k] = j;
                    }
                }
            }

            // 全选主元交换行和列
            // 列交换
            if (C[k] != k)
            {
                for (std::size_t i = 0; i < matrix.row; i++)
                {
                    temp = matrix.mat_data[i * column + k];
                    matrix.mat_data[i * column + k] = matrix.mat_data[i * column + C[k]];
                    matrix.mat_data[i * column + C[k]] = temp;
                }
            }

            // 行交换
            if (R[k] != k)
            {
                // 矩阵matrix进行行交换
                for (std::size_t i = 0; i < matrix.column; i++)
                {
                    temp = matrix.mat_data[k * column + i];
                    matrix.mat_data[k * column + i] = matrix.mat_data[R[k] * column + i];
                    matrix.mat_data[R[k] * column + i] = temp;
                }

                // 向量vec进行行交换
                if (vec != Matrix())
                {
                    temp = vec.mat_data[k];
                    vec.mat_data[k] = vec.mat_data[R[k]];
                    vec.mat_data[R[k]] = temp;
                }
            }

            // 主元为0则主元右下角元素均为0, 高斯消元完成
            if (matrix.mat_data[k*column + k] == 0)
            {
                break;
            }

            // 消元
            for (std::size_t i = k + 1; i < matrix.row; i++)
            {
                // 矩阵matrix消元
                temp = matrix.mat_data[i*column + k] / matrix.mat_data[k*column + k];
                matrix.mat_data[i*column + k] = 0;
                for (std::size_t j = k + 1; j < matrix.column; j++)
                {
                    matrix.mat_data[i*column + j] -= temp * matrix.mat_data[k*column + j];
                }

                // 向量vec消元
                if (vec != Matrix())
                {
                    vec.mat_data[i] -= temp * vec.mat_data[k];
                }
            }
        }

        // 主对角线非零元素个数为矩阵的秩
        std::size_t count = 0;
        for (std::size_t i = 0; i < size; i++)
        {
            if (matrix.mat_data[i*column + i] != 0)
            {
                count++;
            }
        }

        return std::make_tuple(matrix, std::move(R), std::move(C), count, vec);
    }

    // 矩阵拟上三角化
    Matrix Matrix::hess2()const
    {
        if (this->row > 0 && this->row == this->column)
        {
            Matrix A(*this);

            std::size_t n = A.row;   // 矩阵维数
            double d_r = 0;
            double c_r = 0;
            double h_r = 0;
            double* u_r = new double[n];
            Matrix U_R;
            Matrix P_R;
            Matrix Q_R;
            double t_r;

            for (std::size_t r = 0; r < n - 2; r++)
            {
                // 判断A(r+2,r) ~ A(n,r)是否全为0
                bool zero = true;
                for (std::size_t i = r + 2; i < n; i++)
                {
                    if (A[i][r] != 0)
                    {
                        zero = false;
                        break;
                    }
                }

                if (!zero) // 如果A(r+1,r) ~ A(n,r)不全为0
                {
                    // 计算d_r
                    d_r = 0;
                    for (std::size_t i = r + 1; i < n; i++)
                    {
                        d_r += A[i][r] * A[i][r];
                    }
                    d_r = sqrt(d_r);

                    // 计算c_r
                    c_r = -sgn(A[r + 1][r]) * d_r;
                    // 计算h_r
                    h_r = c_r * (c_r - A[r + 1][r]);
                    // 计算向量u_r
                    for (std::size_t i = 0; i < n; i++)
                    {
                        if (i < r + 1)
                            u_r[i] = 0;
                        else
                            u_r[i] = A[i][r];
                    }
                    u_r[r + 1] -= c_r;

                    // 利用数组u_r构造列向量U_R
                    U_R.row = A.row;
                    U_R.column = 1;
                    U_R.mat_data = u_r;

                    // 进行一次Q、A的迭代
                    //Q -= Q * U_R * U_R.trans() / h_r;
                    P_R = A.trans() * U_R / h_r;
                    Q_R = A * U_R / h_r;
                    t_r = (P_R.trans()*U_R / h_r)[0][0];
                    A -= (Q_R - t_r * U_R)*U_R.trans() + U_R * P_R.trans();
                    U_R.mat_data = nullptr;

                    // 消除上三角阵R主对角线以下元素的截断误差
                    for (std::size_t i = r + 2; i < n; i++)
                    {
                        A[i][r] = 0;
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
        double num = (a + d)*(a + d) - 4 * (a*d - b * c);
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
        std::size_t m = A.row;
        if (m < 2)
        {
            throw std::length_error("The matrix dimension should be larger than 2.");
            return A;
        }
        double s = A[m - 2][m - 2] + A[m - 1][m - 1];
        double t = A[m - 2][m - 2] * A[m - 1][m - 1] - A[m - 2][m - 1] * A[m - 1][m - 2];
        Matrix M = A * A - s * A + t * eye(A.row, A.row);
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

        for (std::size_t r = 0; r < m - 1; r++)
        {
            // 判断A(r+1,r) ~ A(m,r)是否全为0
            bool zero = true;
            for (std::size_t i = r + 1; i < m; i++)
            {
                if (B[i][r] != 0)
                {
                    zero = false;
                    break;
                }
            }

            if (!zero) // 如果B(r+1,r) ~ A(m,r)不全为0
            {
                // 计算d_r
                d_r = 0;
                for (std::size_t i = r; i < m; i++)
                {
                    d_r += B[i][r] * B[i][r];
                }
                d_r = sqrt(d_r);

                // 计算c_r
                c_r = -sgn(B[r][r]) * d_r;
                // 计算h_r
                h_r = c_r * (c_r - B[r][r]);
                // 计算向量u_r
                for (std::size_t i = 0; i < m; i++)
                {
                    if (i < r)
                        u_r[i] = 0;
                    else
                        u_r[i] = B[i][r];
                }
                u_r[r] -= c_r;

                // 利用数组u_r构造列向量U_R
                U_R.row = m;
                U_R.column = 1;
                U_R.mat_data = u_r;

                // 进行一次Q、A的迭代
                V_R = B.trans() * U_R / h_r;
                B -= U_R * V_R.trans();
                P_R = A.trans() * U_R / h_r;
                Q_R = A * U_R / h_r;
                t_r = (P_R.trans() * U_R / h_r).mat_data[0];
                A -= (Q_R - t_r * U_R)*U_R.trans() + U_R * P_R.trans();
                U_R.mat_data = nullptr;
            }
        }
        delete[] u_r;
        return A;
    }
}