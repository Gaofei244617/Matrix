#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "mat_def.h"

namespace mat
{
    const std::pair<std::size_t, std::size_t> size(const Matrix& mat);                     // 获取矩阵行列数
    std::vector<double> getDiag(const Matrix& mat);                                        // 获取对角元素
    double normOne(const Matrix& mat);                                                     // 1范数,列和范数
    double normInf(const Matrix& mat);                                                     // 无穷范数,行和范数
    Matrix trans(const Matrix& mat);                                                       // 矩阵转置
    std::size_t rank(const Matrix& mat);                                                   // 矩阵的秩
    double trace(const Matrix& mat);                                                       // 矩阵的迹
    Matrix inv(const Matrix& mat);                                                         // 逆矩阵
    double det(const Matrix& mat);                                                         // 矩阵行列式
    std::vector<Matrix> QR(const Matrix& mat);                                             // 矩阵QR分解
    std::vector<std::complex<double>> eigs(const Matrix& mat, double e = 0);               // 矩阵特征值
    std::vector<Matrix> SVD(const Matrix& mat);                                            // 奇异值分解，返回S、V、D三个矩阵

    Matrix eye(const std::size_t& m, const std::size_t& n);                                // 单位矩阵
    Matrix ones(const std::size_t& m, const std::size_t& n);                               // 元素全为1的矩阵
    Matrix diag(const std::initializer_list<double>& nums);                                // 以向量为对角元素生成方阵
    Matrix rbind(const std::initializer_list<Matrix>& M);                                  // [M1; M2; ...], 需要列数相等
    Matrix cbind(const std::initializer_list<Matrix>& M);                                  // [M1, M2, ...], 需要行数相等

    std::tuple<Matrix, int> solve(const Matrix& A, const Matrix& b);                       // 求解线性方程组 Ax = b
}

#endif
