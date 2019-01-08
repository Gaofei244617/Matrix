﻿// 基于C++11标准
#ifndef _MAT_DEF_H_
#define _MAT_DEF_H_

#include <string>
#include <vector>
#include <tuple>
#include <memory>
#include <cmath>
#include <complex>
#include "math_funcs.h"

namespace mat
{
    using usize = std::size_t;

    class Matrix
    {
    private:
        std::size_t row;                                                           // 行数
        std::size_t column;                                                        // 列数
        double* mat_data;                                                          // 采用一维数组存储元素

    public:
        // 构造函数
        Matrix();
        Matrix(const std::size_t& m, const std::size_t& n);                        // 常规构造: m行数，n列数
        Matrix(const std::initializer_list<double>& m);                            // 列表构造
        Matrix(const std::initializer_list<std::initializer_list<double>>& m);     // 列表构造
        Matrix(const Matrix& other);                                               // 拷贝构造
        Matrix(Matrix&& other);                                                    // 移动构造
        virtual ~Matrix();                                                         // 析构函数

        // 赋值运算
        Matrix& operator=(const Matrix& other)noexcept;                            // M1 = M2
        Matrix& operator=(Matrix&& other)noexcept;                                 // M1 = move(M2)
        Matrix& operator+=(const double& num);                                     // M += num
        Matrix& operator+=(const Matrix& m);                                       // M += M
        Matrix& operator-=(const double& num);                                     // M -= num
        Matrix& operator-=(const Matrix& m);                                       // M -= M
        Matrix& operator*=(const double& num);                                     // M *= num
        Matrix& operator*=(const Matrix& m);                                       // M *= M
        Matrix& operator/=(const double& num);                                     // M /= num

        // 下标索引
        double *const operator[](const std::size_t& n)const noexcept
        {
            return mat_data + n * column;
        }
        double& at(std::size_t m, std::size_t n);
        const double& at(std::size_t m, std::size_t n)const;

        // 点乘(除)运算
        Matrix dotMult(const Matrix& m)const;                                      // M .* M
        Matrix& dotMult(Matrix&& m)const;                                          // M .* move(M)
        Matrix dotDiv(const Matrix& m)const;                                       // M ./ M
        Matrix& dotDiv(Matrix&& m)const;                                           // M ./ move(M)

        const std::pair<std::size_t, std::size_t> size()const;                     // 获取矩阵行列数
        Matrix getRow(std::size_t n)const;                                         // 获取行向量
        Matrix getColumn(std::size_t n)const;                                      // 获取列向量
        std::vector<double> getDiag()const;                                        // 获取对角元素
        double normOne()const;                                                     // 1范数,列和范数
        double normInf()const;                                                     // 无穷范数,行和范数
        Matrix trans()const;                                                       // 矩阵转置

        std::size_t rank()const;                                                   // 矩阵的秩
        double trace()const;                                                       // 矩阵的迹
        Matrix inv()const;                                                         // 逆矩阵
        double det()const;                                                         // 矩阵行列式
        std::pair<Matrix, Matrix> QR()const;                                       // 矩阵QR分解
        std::pair<Matrix, Matrix> LU()const;                                       // 矩阵LU分解
        std::vector<std::complex<double>> eigs(double e = 0)const;                 // 矩阵特征值
        std::tuple<Matrix, Matrix, Matrix> SVD()const;                             // 奇异值分解，返回S、V、D三个矩阵
        Matrix subMat(std::size_t r1, std::size_t c1, std::size_t r2, std::size_t c2)const;     // 子阵

    private:
        // 矩阵拟上三角分解(A = P*B*P'),P为正交矩阵,B为拟上三角阵
        Matrix hess2()const;

        // 计算二阶矩阵特征值(按行输入)
        std::vector<std::complex<double>> eigVal22(double a, double b, double c, double d)const;

        // QR方法计算特征值中的矩阵迭代
        Matrix& iterM(Matrix& A)const;

        // 全选主元高斯消去法
        static std::tuple<Matrix, std::unique_ptr<std::size_t[]>, std::unique_ptr<std::size_t[]>, std::size_t, Matrix>
            gaussElimination(const Matrix& A, const Matrix& b = Matrix());

        /****************************************** Friend Functions ***************************************************************/
        friend std::pair<Matrix, int> solve(const Matrix& A, const Matrix& b);

        // 加法运算
        friend Matrix operator+(const Matrix& m);                                  // +M
        friend Matrix operator+(const Matrix& m1, const Matrix& m2);               // M + M
        friend Matrix& operator+(const Matrix& m1, Matrix&& m2);                   // M + move(M)
        friend Matrix& operator+(Matrix&& m1, const Matrix& m2);                   // move(M) + M
        friend Matrix& operator+(Matrix&& m1, Matrix&& m2);                        // move(M) + move(M)
        friend Matrix operator+(const Matrix& m, const double& num);               // M + num
        friend Matrix operator+(const double& num, const Matrix& m);               // num + M
        friend Matrix& operator+(const double& num, Matrix&& m);                   // num + move(M)
        friend Matrix& operator+(Matrix&& m, const double& num);                   // move(M) + num

        // 减法运算
        friend Matrix operator-(const Matrix& m);                                  // -M
        friend Matrix operator-(const Matrix& m1, const Matrix& m2);               // M - M
        friend Matrix& operator-(const Matrix& m1, Matrix&& m);                    // M - move(M)
        friend Matrix& operator-(Matrix&& m1, const Matrix& m2);                   // move(M) - M
        friend Matrix& operator-(Matrix&& m1, Matrix&& m2);                        // move(M) - move(M)
        friend Matrix operator-(const Matrix& m, const double& num);               // M - num
        friend Matrix operator-(const double& num, const Matrix& m);               // num - M
        friend Matrix& operator-(const double& num, Matrix&& m);                   // num - move(M)
        friend Matrix& operator-(Matrix&& m, const double& num);                   // move(M) - num

        // 乘法运算
        friend Matrix operator*(const Matrix& m1, const Matrix& m2);               // M * M
        friend Matrix operator*(const Matrix& m, const double& num);               // M * num
        friend Matrix& operator*(Matrix&& m, const double& num);                   // move(M) * num
        friend Matrix& operator*(Matrix&& m1, Matrix&& m2);                        // move(M) * move(M)
        friend Matrix operator*(const double& num, const Matrix& m);               // num * M
        friend Matrix& operator*(const double& num, Matrix&& m);                   // num * move(M)

        // 除法运算
        friend Matrix operator/(const Matrix& m, const double& num);               // M / num
        friend Matrix& operator/(Matrix&& m, const double& num);                   // move(M) / num

        // 逻辑运算符
        friend bool operator==(const Matrix& m1, const Matrix& m2)noexcept;        // 等于
        friend bool operator!=(const Matrix& m1, const Matrix& m2)noexcept;        // 不等于
    };
}

#endif
