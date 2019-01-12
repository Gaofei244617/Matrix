#include "matrix.h"

namespace mat
{
    // 获取矩阵行列数
    const std::pair<usize, usize> size(const Matrix& mat)
    {
        return mat.size();
    }

    // 获取行向量
    Matrix getRow(const Matrix& mat, usize n)
    {
        return mat.getRow(n);
    }

    // 获取列向量
    Matrix getColumn(const Matrix& mat, usize n)
    {
        return mat.getColumn(n);
    }

    // 获取对角元素
    std::vector<double> getDiag(const Matrix& mat)
    {
        return mat.getDiag();
    }

    // 1范数,列和范数
    double normOne(const Matrix& mat)
    {
        return mat.normOne();
    }

    // 无穷范数,行和范数
    double normInf(const Matrix& mat)
    {
        return mat.normInf();
    }

    // 矩阵转置
    Matrix trans(const Matrix& mat)
    {
        return mat.trans();
    }

    // 矩阵的秩
    usize rank(const Matrix& mat)
    {
        return mat.rank();
    }

    // 矩阵的迹
    double trace(const Matrix& mat)
    {
        return mat.trace();
    }

    // 逆矩阵
    Matrix inv(const Matrix& mat)
    {
        return mat.inv();
    }

    // 矩阵行列式
    double det(const Matrix& mat)
    {
        return mat.det();
    }

    // 矩阵QR分解
    std::pair<Matrix, Matrix> QR(const Matrix& mat)
    {
        return mat.QR();
    }

    // 矩阵LU分解,返回值{P,L,U}
    std::tuple<Matrix, Matrix, Matrix> LU(const Matrix& mat)
    {
        return mat.LU();
    }

    // 矩阵特征值
    std::vector<std::complex<double>> eig(const Matrix& mat, double e)
    {
        return mat.eig(e);
    }

    // 子阵
    Matrix subMat(const Matrix& mat, usize r1, usize c1, usize r2, usize c2)
    {
        return mat.subMat(r1, c1, r2, c2);
    }

    // 单位矩阵
    Matrix eye(const usize& m, const usize& n)
    {
        Matrix matrix(m, n);
        usize size = 0;
        m < n ? size = m : size = n;
        for (usize i = 0; i < size; i++)
        {
            matrix[i][i] = 1;
        }
        return matrix;
    }

    // 元素全为1的矩阵
    Matrix ones(const usize& m, const usize& n)
    {
        return Matrix(m, n) + 1.0;
    }

    // 元素全为0的矩阵
    Matrix zeros(const usize& m, const usize& n)
    {
        return Matrix(m, n);
    }

    // 随机矩阵, 元素取值[0, 1.0)
    Matrix rand(const usize& m, const usize& n)
    {
        Matrix M(m, n);
        for (usize i = 0; i < m; i++)
        {
            for (usize j = 0; j < n; j++)
            {
                M[i][j] = random_real();
            }
        }
        return M;
    }

    // 以向量为对角元素生成方阵
    Matrix diag(const std::initializer_list<double>& nums)
    {
        usize size = nums.size();
        Matrix matrix(size, size);
        for (usize i = 0; i < size; i++)
        {
            matrix[i][i] = *(nums.begin() + i);
        }
        return matrix;
    }

    // [M1; M2; ...], 需要列数相等
    Matrix rbind(const std::initializer_list<Matrix>& M)
    {
        usize count = M.size();
        if (count > 0)
        {
            usize row;       // 拼接后矩阵的行数
            usize column;    // 拼接后矩阵的列数（拼接前后保持不变）

            usize temp_row = 0;
            usize temp_column = 0;

            std::tie(row, column) = (*(M.begin())).size(); // 第一个矩阵的行列数

            // 计算拼接后矩阵的行数, 并判断矩阵是否可拼接
            for (usize i = 1; i < count; i++)
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
            for (usize i = 0; i < count; i++)
            {
                temp_row = std::get<0>((*(M.begin() + i)).size());

                for (usize j = 0; j < temp_row; j++)
                {
                    for (usize k = 0; k < column; k++)
                    {
                        matrix[j + p][k] = (*(M.begin() + i))[j][k];
                    }
                }
                p += temp_row;
            }
            return matrix;
        }
        return Matrix();
    }

    // [M1, M2, ...], 需要行数相等
    Matrix cbind(const std::initializer_list<Matrix>& M)
    {
        usize count = M.size();
        if (count > 0)
        {
            usize row;       // 拼接后矩阵的行数（拼接前后保持不变）
            usize column;    // 拼接后矩阵的列数

            usize temp_row = 0;
            usize temp_column = 0;

            std::tie(row, column) = (*(M.begin())).size(); // 第一个矩阵的行列数

            // 判断矩阵是否可拼接,并计算拼接后矩阵的行数
            for (usize i = 1; i < count; i++)
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
            for (usize i = 0; i < count; i++)
            {
                temp_column = std::get<1>((*(M.begin() + i)).size());

                for (usize j = 0; j < row; j++)
                {
                    for (usize k = 0; k < temp_column; k++)
                    {
                        matrix[j][k + p] = (*(M.begin() + i))[j][k];
                    }
                }

                p += temp_column;
            }
            return matrix;
        }

        return Matrix();
    }

    // 求解线性方程组: Ax = b
    std::pair<Matrix, int> solve(const Matrix& A, const Matrix& b)
    {
        /*******************************************************/
        // 求解线性方程组: Ax = b
        // A为系数矩阵，b为目标向量
        // 采用全选主元的高斯消元法求解方程组

        // (1)若方程有唯一解：   返回 pair(解向量x,  标志量0)
        // (2)若方程有无穷多解： 返回 pair(一个特解x, 标志量1)
        // (3)若方程无精确解：   返回 pair(近似解x*,  标志量2),  近似解x*使得|Ax* - b|最小
        /********************************************************/

        // 判断系数矩阵的行数与目标向量的行数是否相等
        if (A.row != b.row)
        {
            throw std::length_error("Size of coefficient matrix and target vector does not match.");
            return std::make_pair(Matrix(), 2);
        }

        auto temp = Matrix::gaussElimination(A, b);                // 对系数矩阵进行全选主元的高斯消元
        Matrix& matrix = std::get<0>(temp);                        // 全选主元高斯消元后的矩阵
        std::unique_ptr<usize[]>& R = std::get<1>(temp);           // 行交换信息
        std::unique_ptr<usize[]>& C = std::get<2>(temp);           // 列交换信息
        usize rankA = std::get<3>(temp);                           // 系数矩阵的秩rank(A)
        usize rankAb = cbind({ A, b }).rank();                     // rank(A,b)
        usize n = A.column;                                        // 未知数个数
        Matrix& vec = std::get<4>(temp);                           // 高斯消元后方程组的目标向量

        /* rank(A) == rank(A,b) && rank(A) == n: 方程有唯一解 */
        /* rank(A) == rank(A,b) && rank(A) < n:  方程有无穷多解 */
        /* others: 方程无解( rank(A)不可能大于n ) */

        // (1) 方程有解
        if (rankA == rankAb)
        {
            Matrix x(n, 1);   // 方程的解

            // 求解高斯消元后的上三角方程组
            for (int i = rankA - 1; i >= 0; i--)
            {
                for (int j = 0; j < rankA - i - 1; j++)
                {
                    vec[i][0] -= matrix[i][rankA - 1 - j] * x[rankA - 1 - j][0];
                }

                x[i][0] = 1.0 * vec[i][0] / matrix[i][i];
            }

            // 方程有无穷多个解时，另解向量x部分元素为0，得到一个特解
            if (rankA < n)
            {
                for (usize i = rankA; i < n; i++)
                {
                    x[i][0] = 0;
                }
            }

            // 依据列交换信息,恢复向量x中元素的顺序
            auto temp = x[rankA - 1][0];
            for (int i = rankA - 1; i >= 0; i--)
            {
                if (C[i] != i)
                {
                    temp = x[i][0];
                    x[i][0] = x[C[i]][0];
                    x[C[i]][0] = temp;
                }
            }

            // (1.1) 方程有唯一解
            if (rankA == n)
            {
                return std::make_pair(x, 0);
            }

            // (1.2) 方程有无穷多个解
            if (rankA < n)
            {
                return std::make_pair(x, 1);
            }
        }

        // (2)方程无解
        if (rankA != rankAb)
        {
            Matrix AT(A.trans());
            return std::make_pair(std::get<0>(solve(AT*A, AT*b)), 2);
        }

        return std::make_pair(Matrix(), 2);
    }
}