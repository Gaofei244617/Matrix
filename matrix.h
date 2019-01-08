﻿#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "mat_def.h"

namespace mat
{
    const std::pair<usize, usize> size(const Matrix& mat);                                 // 获取矩阵行列数
    Matrix getRow(const Matrix& mat, usize n);                                             // 获取行向量
    Matrix getColumn(const Matrix& mat, usize n);                                          // 获取列向量
    std::vector<double> getDiag(const Matrix& mat);                                        // 获取对角元素
    double normOne(const Matrix& mat);                                                     // 1范数,列和范数
    double normInf(const Matrix& mat);                                                     // 无穷范数,行和范数
    Matrix trans(const Matrix& mat);                                                       // 矩阵转置
    usize rank(const Matrix& mat);                                                         // 矩阵的秩
    double trace(const Matrix& mat);                                                       // 矩阵的迹
    Matrix inv(const Matrix& mat);                                                         // 逆矩阵
    double det(const Matrix& mat);                                                         // 矩阵行列式
    std::pair<Matrix, Matrix> QR(const Matrix& mat);                                       // 矩阵QR分解,返回值{Q,R}
    std::tuple<Matrix, Matrix, Matrix> LU(const Matrix& mat);                              // 矩阵LU分解,返回值{P,L,U}
    std::vector<std::complex<double>> eig(const Matrix& mat, double e = 0);                // 矩阵特征值
    std::tuple<Matrix, Matrix, Matrix> SVD(const Matrix& mat);                             // 奇异值分解，返回S、V、D三个矩阵
    Matrix subMat(const Matrix& mat, usize r1, usize c1, usize r2, usize c2);              // 子阵

    Matrix eye(const usize& m, const usize& n);                                            // 单位矩阵
    Matrix ones(const usize& m, const usize& n);                                           // 元素全为1的矩阵
    Matrix diag(const std::initializer_list<double>& nums);                                // 以向量为对角元素生成方阵
    Matrix rbind(const std::initializer_list<Matrix>& M);                                  // [M1; M2; ...], 需要列数相等
    Matrix cbind(const std::initializer_list<Matrix>& M);                                  // [M1, M2, ...], 需要行数相等

    std::pair<Matrix, int> solve(const Matrix& A, const Matrix& b);                        // 求解线性方程组 Ax = b

    /***********************************************************************************************************************/
    // 加法运算
    Matrix operator+(const Matrix& m);                                  // +M
    Matrix operator+(const Matrix& m1, const Matrix& m2);               // M + M
    Matrix& operator+(const Matrix& m1, Matrix&& m2);                   // M + move(M)
    Matrix& operator+(Matrix&& m1, const Matrix& m2);                   // move(M) + M
    Matrix& operator+(Matrix&& m1, Matrix&& m2);                        // move(M) + move(M)
    Matrix operator+(const Matrix& m, const double& num);               // M + num
    Matrix operator+(const double& num, const Matrix& m);               // num + M
    Matrix& operator+(const double& num, Matrix&& m);                   // num + move(M)
    Matrix& operator+(Matrix&& m, const double& num);                   // move(M) + num

    // 减法运算
    Matrix operator-(const Matrix& m);                                  // -M
    Matrix operator-(const Matrix& m1, const Matrix& m2);               // M - M
    Matrix& operator-(const Matrix& m1, Matrix&& m);                    // M - move(M)
    Matrix& operator-(Matrix&& m1, const Matrix& m2);                   // move(M) - M
    Matrix& operator-(Matrix&& m1, Matrix&& m2);                        // move(M) - move(M)
    Matrix operator-(const Matrix& m, const double& num);               // M - num
    Matrix operator-(const double& num, const Matrix& m);               // num - M
    Matrix& operator-(const double& num, Matrix&& m);                   // num - move(M)
    Matrix& operator-(Matrix&& m, const double& num);                   // move(M) - num

    // 乘法运算
    Matrix operator*(const Matrix& m1, const Matrix& m2);               // M * M
    Matrix operator*(const Matrix& m, const double& num);               // M * num
    Matrix& operator*(Matrix&& m, const double& num);                   // move(M) * num
    Matrix& operator*(Matrix&& m1, Matrix&& m2);                        // move(M) * move(M)
    Matrix operator*(const double& num, const Matrix& m);               // num * M
    Matrix& operator*(const double& num, Matrix&& m);                   // num * move(M)

    // 除法运算
    Matrix operator/(const Matrix& m, const double& num);               // M / num
    Matrix& operator/(Matrix&& m, const double& num);                   // move(M) / num

    // 逻辑运算符
    bool operator==(const Matrix& m1, const Matrix& m2)noexcept;        // 等于
    bool operator!=(const Matrix& m1, const Matrix& m2)noexcept;        // 不等于
}

#endif
