#include "mat_operator.h"

namespace mat
{
    // 加法运算
    Matrix operator+(const Matrix& m)
    {
        return m;
    }
    Matrix operator+(const Matrix& m1, const Matrix& m2)
    {
        if (m1.row == m2.row && m1.column == m2.column)
        {
            usize size = m1.row * m1.column;
            Matrix matrix(m1.row, m1.column);
            for (usize i = 0; i < size; i++)
            {
                matrix.mat_data[i] = m1.mat_data[i] + m2.mat_data[i];
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
            usize size = m1.row * m1.column;
            for (usize i = 0; i < size; i++)
            {
                m2.mat_data[i] += m1.mat_data[i];
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
            usize size = m.row * m.column;
            for (usize i = 0; i < size; i++)
            {
                matrix.mat_data[i] = m.mat_data[i] + num;
            }
            return matrix;
        }
        else
        {
            throw std::length_error("Null Matrix.");
        }
        return m;
    }
    Matrix& operator+(Matrix&& m1, const Matrix& m2)
    {
        if (m1.row == m2.row && m1.column == m2.column)
        {
            usize size = m1.row * m1.column;
            for (usize i = 0; i < size; i++)
            {
                m1.mat_data[i] += m2.mat_data[i];
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
            usize size = m.row * m.column;
            for (usize i = 0; i < size; i++)
            {
                matrix.mat_data[i] = m.mat_data[i] + num;
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
            usize size = m.row * m.column;
            for (usize i = 0; i < size; i++)
            {
                m.mat_data[i] += num;
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
            usize size = m.row * m.column;
            for (usize i = 0; i < size; i++)
            {
                m.mat_data[i] += num;
            }
        }
        else
        {
            throw std::length_error("Null Matrix.");
        }
        return m;
    }

    // 减法运算
    Matrix operator-(const Matrix& m)
    {
        Matrix matrix(m.row, m.column);
        usize size = m.row * m.column;
        for (usize i = 0; i < size; i++)
        {
            matrix.mat_data[i] = -(m.mat_data[i]);
        }
        return matrix;
    }
    Matrix operator-(const Matrix& m1, const Matrix& m2)
    {
        if (m1.row == m2.row && m1.column == m2.column)
        {
            usize size = m1.row * m1.column;
            Matrix matrix(m1.row, m1.column);
            for (usize i = 0; i < size; i++)
            {
                matrix.mat_data[i] = m1.mat_data[i] - m2.mat_data[i];
            }
            return matrix;
        }
        else
        {
            throw std::length_error("Dimensions do not match.");
        }
        return m1;
    }
    Matrix& operator-(const Matrix& m1, Matrix&& m2)
    {
        if (m1.row == m2.row && m1.column == m2.column)
        {
            usize size = m1.row * m1.column;
            for (usize i = 0; i < size; i++)
            {
                m2.mat_data[i] = m1.mat_data[i] - m2.mat_data[i];
            }
            return m2;
        }
        else
        {
            throw std::length_error("Dimensions do not match.");
        }
        return m2 = m1;
    }
    Matrix operator-(const Matrix& m, const double& num)
    {
        if (m.row > 0)
        {
            Matrix matrix(m.row, m.column);
            usize size = m.row * m.column;
            for (usize i = 0; i < size; i++)
            {
                matrix.mat_data[i] = m.mat_data[i] - num;
            }
            return matrix;
        }
        else
        {
            throw std::length_error("Null Matrix.");
        }
        return m;
    }
    Matrix& operator-(Matrix&& m1, const Matrix& m2)
    {
        if (m1.row == m2.row && m1.column == m2.column)
        {
            usize size = m1.row * m1.column;
            for (usize i = 0; i < size; i++)
            {
                m1.mat_data[i] -= m2.mat_data[i];
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
            usize size = m.row * m.column;
            for (usize i = 0; i < size; i++)
            {
                matrix.mat_data[i] = num - m.mat_data[i];
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
            usize size = m.row * m.column;
            for (usize i = 0; i < size; i++)
            {
                m.mat_data[i] = num - m.mat_data[i];
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
            usize size = m.row * m.column;
            for (usize i = 0; i < size; i++)
            {
                m.mat_data[i] -= num;
            }
        }
        else
        {
            throw std::length_error("Null Matrix.");
        }
        return m;
    }

    // 乘法运算
    Matrix operator*(const Matrix& m1, const Matrix& m2)
    {
        if (std::get<1>(m1.size()) == std::get<0>(m2.size()))
        {
            usize row = std::get<0>(m1.size());
            usize column = std::get<1>(m2.size());
            Matrix matrix(row, column);

            usize count = std::get<0>(m2.size());
            for (usize i = 0; i < row; i++)
            {
                for (usize j = 0; j < column; j++)
                {
                    double sum = 0;
                    for (usize k = 0; k < count; k++)
                    {
                        sum += m1[i][k] * m2[k][j];
                    }
                    matrix[i][j] = sum;
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
    Matrix operator*(const Matrix& m, const double& num)
    {
        Matrix matrix(m.row, m.column);
        usize size = m.row * m.column;
        for (usize i = 0; i < size; i++)
        {
            matrix.mat_data[i] = m.mat_data[i] * num;
        }
        return matrix;
    }
    Matrix& operator*(Matrix&& m, const double& num)
    {
        usize size = m.row * m.column;
        for (usize i = 0; i < size; i++)
        {
            m.mat_data[i] *= num;
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
        usize size = m.row * m.column;
        for (usize i = 0; i < size; i++)
        {
            matrix.mat_data[i] = num * m.mat_data[i];
        }
        return matrix;
    }
    Matrix& operator*(const double& num, Matrix&& m)
    {
        usize size = m.row * m.column;
        for (usize i = 0; i < size; i++)
        {
            m.mat_data[i] *= num;
        }
        return m;
    }

    // 除法运算
    Matrix operator/(const Matrix& m, const double& num)
    {
        if (m.row > 0)
        {
            if (num != 0)
            {
                Matrix matrix(m.row, m.column);
                usize size = m.row * m.column;
                for (usize i = 0; i < size; i++)
                {
                    matrix.mat_data[i] = m.mat_data[i] / num;
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
        return m;
    }
    Matrix& operator/(Matrix&& m, const double& num)
    {
        if (m.row > 0)
        {
            if (num != 0)
            {
                usize size = m.row * m.column;
                for (usize i = 0; i < size; i++)
                {
                    m.mat_data[i] /= num;
                }
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
        return m;
    }

    // 逻辑运算
    bool operator==(const Matrix& m1, const Matrix& m2)noexcept
    {
        if (m1.row == m2.row && m1.column == m2.column)
        {
            bool flag = true;
            usize size = m1.row * m1.column;
            for (usize i = 0; i < size; i++)
            {
                if (m1.mat_data[i] != m2.mat_data[i])
                {
                    flag = false;
                    break;
                }
            }
            return flag;
        }
        return false;
    }
    bool operator!=(const Matrix& m1, const Matrix& m2)noexcept
    {
        return !(m1 == m2);
    }
}