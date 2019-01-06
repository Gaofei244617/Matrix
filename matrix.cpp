#include "matrix.h"

namespace mat
{
    Matrix eye(const std::size_t& m, const std::size_t& n)
    {
        Matrix matrix(m, n);
        std::size_t size = 0;
        m < n ? size = m : size = n;
        for (std::size_t i = 0; i < size; i++)
        {
            matrix(i, i) = 1;
        }
        return matrix;
    }

    Matrix ones(const std::size_t& m, const std::size_t& n)
    {
        return Matrix(m, n) + 1;
    }

    Matrix diag(const std::initializer_list<double>& nums)
    {
        std::size_t size = nums.size();
        Matrix matrix(size, size);
        for (std::size_t i = 0; i < size; i++)
        {
            matrix(i, i) = *(nums.begin() + i);
        }
        return matrix;
    }

    Matrix rbind(const std::initializer_list<Matrix>& M)
    {
        std::size_t count = M.size();
        if (count > 0)
        {
            std::size_t row;       // 拼接后矩阵的行数
            std::size_t column;    // 拼接后矩阵的列数（拼接前后保持不变）

            std::size_t temp_row = 0;
            std::size_t temp_column = 0;

            std::tie(row, column) = (*(M.begin())).size(); // 第一个矩阵的行列数

            // 计算拼接后矩阵的行数, 并判断矩阵是否可拼接
            for (std::size_t i = 1; i < count; i++)
            {
                std::tie(temp_row, temp_column) = (*(M.begin() + i)).size();

                row += temp_row;

                // 列数不相同, 不能拼接
                if (temp_column != column)
                {
                    throw std::length_error("Dimensions do not match.");
                    return Matrix();
                }
            }

            // 拼接矩阵 [M1; M2; ...]
            Matrix matrix(row, column);

            // 进行矩阵拼接
            size_t p = 0;
            for (std::size_t i = 0; i < count; i++)
            {
                temp_row = std::get<0>((*(M.begin() + i)).size());

                for (std::size_t j = 0; j < temp_row; j++)
                {
                    for (std::size_t k = 0; k < column; k++)
                    {
                        matrix(j + p, k) = (*(M.begin() + i))(j, k);
                    }
                }

                p += temp_row;
            }

            return matrix;
        }

        return Matrix();
    }

    Matrix cbind(const std::initializer_list<Matrix>& M)
    {
        std::size_t count = M.size();
        if (count > 0)
        {
            std::size_t row;       // 拼接后矩阵的行数（拼接前后保持不变）
            std::size_t column;    // 拼接后矩阵的列数

            std::size_t temp_row = 0;
            std::size_t temp_column = 0;

            std::tie(row, column) = (*(M.begin())).size(); // 第一个矩阵的行列数

            // 判断矩阵是否可拼接,并计算拼接后矩阵的行数
            for (std::size_t i = 1; i < count; i++)
            {
                std::tie(temp_row, temp_column) = (*(M.begin() + i)).size();

                column += temp_column;
                // 行数不同，不能拼接
                if (temp_row != row)
                {
                    throw std::length_error("Dimensions do not match.");
                    return Matrix();
                }
            }

            Matrix matrix(row, column);

            // 拼接矩阵 [M1, M2, ...]
            size_t p = 0;
            for (std::size_t i = 0; i < count; i++)
            {
                temp_column = std::get<1>((*(M.begin() + i)).size());

                for (std::size_t j = 0; j < row; j++)
                {
                    for (std::size_t k = 0; k < temp_column; k++)
                    {
                        matrix(j, k + p) = (*(M.begin() + i))(j, k);
                    }
                }

                p += temp_column;
            }
            return matrix;
        }

        return Matrix();
    }

    //// 求解线性方程组: Ax = b
    //std::tuple<Matrix, int> solve(const Matrix& A, const Matrix& b)
    //{
    //    /*******************************************************/
    //    // 求解线性方程组: Ax = b
    //    // A为系数矩阵，b为目标向量
    //    // 采用全选主元的高斯消元法求解方程组

    //    // (1)若方程有唯一解：   返回 tuple(解向量x,  标志量0)
    //    // (2)若方程有无穷多解： 返回 tuple(一个特解x, 标志量1)
    //    // (3)若方程无精确解：   返回 tuple(近似解x*,  标志量2)
    //    //    近似解x*使得|Ax* - b|最小
    //    /********************************************************/

    //    // 判断系数矩阵的行数宇目标向量的行数是否相等
    //    if (A.row != b.row)
    //    {
    //        throw std::length_error("Size of coefficient matrix and target vector does not match.");
    //    }

    //    auto temp = gaussElimination(A, b);                        // 对系数矩阵进行全选主元的高斯消元
    //    Matrix& matrix = std::get<0>(temp);                        // 全选主元高斯消元后的矩阵
    //    std::unique_ptr<std::size_t[]>& R = std::get<1>(temp);     // 行交换信息
    //    std::unique_ptr<std::size_t[]>& C = std::get<2>(temp);     // 列交换信息
    //    std::size_t rankA = std::get<3>(temp);                     // 系数矩阵的秩rank(A)
    //    std::size_t rankAb = cbind({ A, b }).rank();               // rank(A,b)
    //    std::size_t n = A.column;                                  // 未知数个数
    //    Matrix& vec = std::get<4>(temp);                           // 高斯消元后方程组的目标向量

    //    // rank(A) == rank(A,b) && rank(A) == n: 方程有唯一解
    //    // rank(A) == rank(A,b) && rank(A) < n:  方程有无穷多解
    //    // others: 方程无解( rank(A)不可能大于n )

    //    // (1) 方程有解
    //    if (rankA == rankAb)
    //    {
    //        Matrix x(n, 1);   // 方程的解

    //        // 求解高斯消元后的上三角方程组
    //        for (std::size_t i = rankA - 1; i >= 0; i--)
    //        {
    //            for (std::size_t j = 0; j < rankA - i - 1; j++)
    //            {
    //                vec(i, 0) -= matrix(i, rankA - 1 - j) * x(rankA - 1 - j, 0);
    //            }

    //            x(i, 0) = 1.0 * vec(i, 0) / matrix(i, i);
    //        }

    //        // 方程有无穷多个解时，另解向量x部分元素为0，得到一个特解
    //        if (rankA < n)
    //        {
    //            for (std::size_t i = rankA; i < n; i++)
    //            {
    //                x(i, 0) = 0;
    //            }
    //        }

    //        // 依据列交换信息,恢复向量x中元素的顺序
    //        auto temp = x(rankA - 1, 0);
    //        for (std::size_t i = rankA - 1; i >= 0; i--)
    //        {
    //            if (C[i] != i)
    //            {
    //                temp = x(i, 0);
    //                x(i, 0) = x(C[i], 0);
    //                x(C[i], 0) = temp;
    //            }
    //        }

    //        // (1.1) 方程有唯一解
    //        if (rankA == n)
    //        {
    //            return std::make_tuple(x, 0);
    //        }

    //        // (1.2) 方程有无穷多个解
    //        if (rankA < n)
    //        {
    //            return std::make_tuple(x, 1);
    //        }
    //    }

    //    // (2)方程无解
    //    if (rankA != rankAb)
    //    {
    //        Matrix AT(A.trans());
    //        return std::make_tuple(std::get<0>(solve(AT*A, AT*b)), 2);
    //    }

    //    return std::make_tuple(Matrix(), 2);
    //}
}