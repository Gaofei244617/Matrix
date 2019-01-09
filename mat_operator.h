#ifndef _MAT_OPERATOR_H_
#define _MAT_OPERATOR_H_

#include "mat_def.h"

namespace mat
{
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
