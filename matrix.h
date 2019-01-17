#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "mat_def.h"

namespace mat
{
    const std::pair<usize, usize> size(const Matrix& mat);                                 // 获取矩阵行列数
    Matrix getRow(const Matrix& mat, usize n);                                             // 获取行向量
    Matrix getColumn(const Matrix& mat, usize n);                                          // 获取列向量
    std::vector<double> getDiag(const Matrix& mat);                                        // 获取对角元素
    double normOne(const Matrix& mat);                                                     // 1范数,列和范数(每一列元素绝对值之和的最大值)
    double normOne(const Matrix& mat);                                                     // 2范数,列和范数(AA'的最大特征值的平方根)
    double normInf(const Matrix& mat);                                                     // 无穷范数,行和范数(每一行元素绝对值之和的最大值)
    Matrix trans(const Matrix& mat);                                                       // 矩阵转置
    usize rank(const Matrix& mat);                                                         // 矩阵的秩
    double trace(const Matrix& mat);                                                       // 矩阵的迹
    Matrix inv(const Matrix& mat);                                                         // 逆矩阵
    double det(const Matrix& mat);                                                         // 矩阵行列式
    std::pair<Matrix, Matrix> QR(const Matrix& mat);                                       // 矩阵QR分解,返回值{Q,R}
    std::tuple<Matrix, Matrix, Matrix> LU(const Matrix& mat);                              // 矩阵LU分解,返回值{P,L,U}
    std::vector<std::complex<double>> eigVal(const Matrix& mat, double e = 0);                // 矩阵特征值
    std::tuple<Matrix, Matrix, Matrix> SVD(const Matrix& mat);                             // 奇异值分解，返回S、V、D三个矩阵

    Matrix eye(const usize& m, const usize& n);                                            // 单位矩阵
    Matrix ones(const usize& m, const usize& n);                                           // 元素全为1的矩阵
    Matrix zeros(const usize& m, const usize& n);                                          // 元素全为0的矩阵
    Matrix rand(const usize& m, const usize& n);                                           // 随机矩阵, 元素取值[0, 1.0)
    Matrix diag(const Matrix& vec);                                                        // 以向量为对角元素生成方阵
    template<class T, class... Args> Matrix rbind(const T& M, const Args&... args);        // [M1; M2; ...], 需要列数相等
    template<class T, class... Args> Matrix cbind(const T& M, const Args&... args);        // [M1, M2, ...], 需要行数相等
    Matrix subMat(const Matrix& mat, usize r1, usize c1, usize r2, usize c2);              // 子阵

    std::pair<Matrix, int> solve(const Matrix& A, const Matrix& b);                        // 求解线性方程组 Ax = b

    /******************************************** 函数模板 ****************************************************************/
    // [M1; M2; ...], 需要列数相等
    template<class T, class... Args>
    Matrix rbind(const T& M, const Args&... args)
    {
        return M.rbind(args...);
    }

    // [M1, M2, ...], 需要行数相等
    template<class T, class... Args>
    Matrix cbind(const T& M, const Args&... args)
    {
        return M.cbind(args...);
    }
}

#endif
