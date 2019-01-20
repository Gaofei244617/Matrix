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
    std::vector<std::complex<double>> eigVal(const Matrix& mat, double e)
    {
        return mat.eigVal(e);
    }

    // 子阵
    Matrix subMat(const Matrix& mat, usize r1, usize c1, usize r2, usize c2)
    {
        return mat.subMat(r1, c1, r2, c2);
    }

    // 矩阵元素是否全部为零
    bool isZero(const Matrix& mat, const double e)
    {
        return mat.isZero(e);
    }

    // 高阶函数-filter
    Matrix filter(const Matrix& mat, std::function<bool(double)> f)
    {
        return mat.filter(std::move(f));
    }

    // 高阶函数-map
    Matrix map(const Matrix& mat, std::function<double(double)> f)
    {
        return mat.map(std::move(f));
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

    // Hilbert矩阵(高度病态)
    Matrix hilb(const usize& m, const usize& n)
    {
        Matrix mat(m, n);
        for (usize i = 0; i < m; i++)
        {
            for (usize j = 0; j < n; j++)
            {
                mat[i][j] = 1.0 / (i + j + 1);
            }
        }
        return mat;
    }

    // 随机矩阵, 元素取值[0, 1.0)
    Matrix rand(const usize& m, const usize& n)
    {
        Matrix mat(m, n);
        for (usize i = 0; i < m; i++)
        {
            for (usize j = 0; j < n; j++)
            {
                mat[i][j] = random_real();
            }
        }
        return mat;
    }

    // 随机矩阵, 元素服从正态分布(均值u,方差t)
    Matrix randn(const usize& m, const usize& n, const double u, const double t)
    {
        Matrix mat(m, n);
        for (usize i = 0; i < m; i++)
        {
            for (usize j = 0; j < n; j++)
            {
                mat[i][j] = random_norm(u, t);
            }
        }
        return mat;
    }

    // 以向量为对角元素生成方阵
    Matrix diag(const Matrix& vec)
    {
        usize row = 0;
        usize column = 0;
        std::tie(row, column) = vec.size();
        if (row == 1 || column == 1)
        {
            const usize N = row > column ? row : column;
            Matrix MAT(N, N);
            if (row == 1)
            {
                for (usize i = 0; i < N; i++)
                {
                    MAT[i][i] = vec[0][i];
                }
            }
            else
            {
                for (usize i = 0; i < N; i++)
                {
                    MAT[i][i] = vec[i][0];
                }
            }
            return MAT;
        }
        else
        {
            throw std::length_error("argument must be a vector.");
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
        const usize rankA = std::get<3>(temp);                     // 系数矩阵的秩rank(A)
        const usize rankAb = A.cbind(b).rank();                    // rank(A,b)

        const usize N = A.column;                                  // 变量数个数
        Matrix vec = std::get<4>(temp);                            // 高斯消元后方程组的目标向量

        /* rank(A) == rank(A,b) && rank(A) == N: 方程有唯一解 */
        /* rank(A) == rank(A,b) && rank(A) < N:  方程有无穷多解 */
        /* others: 方程无解( rank(A)不可能大于N ) */

        // (1) 方程有解
        if (rankA == rankAb)
        {
            // 方程的解向量
            // (1.1) rankA == N : 方程有唯一解
            // (1.2) rankA < N : 方程有无穷多解
            const usize M = N - rankA + 1;  // 解向量的个数
            Matrix x(N, M);

            // 求解高斯消元后的上三角方程组
            for (int i = rankA - 1; i >= 0; i--)
            {
                for (int j = 0; j < rankA - i - 1; j++)
                {
                    vec[i][0] -= matrix[i][rankA - 1 - j] * x[rankA - 1 - j][0];
                }

                x[i][0] = vec[i][0] / matrix[i][i];
            }

            // 依据列交换信息,恢复向量x中元素的顺序
            double temp = 0;
            for (int i = rankA - 1; i >= 0; i--)
            {
                if (C[i] != i)
                {
                    temp = x[i][0];
                    x[i][0] = x[C[i]][0];
                    x[C[i]][0] = temp;
                }
            }

            return std::make_pair(x, 0);
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