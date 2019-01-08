#include <iostream>
#include <windows.h>
#include <string>
#include <memory>
#include <cmath>
#include <complex>
#include <vector>
#include "matrix.h"

using namespace std;
using namespace mat;

class Test
{
private:
    double* num;
    int length;
public:
    Test(int n) :length(n), num(nullptr)
    {
        if (n > 0) num = new double[n]();
        for (int i = 0; i < n; i++)
            num[i] = i;
    }

    Test(const Test& other)
    {
        length = other.length;
        num = new double[other.length];
        for (int i = 0; i < length; i++)
            num[i] = other.num[i];
    }

    double& operator()(int n)
    {
        if (n <= length)
            return num[n];
    }

    void print()
    {
        for (int i = 0; i < length; i++)
            cout << num[i] << "  ";
        cout << endl;
    }

    Test operator+ (const double& a) const
    {
        Test t(length);
        for (int i = 0; i < length; i++)
            t.num[i] = num[i] + a;
        return t;
    }
    Test operator+ () const
    {
        cout << "运算符重载+" << endl;
        return Test(this->length);
    }

    friend Test operator+ (const double& a, Test& other)
    {
        return other + a;
    }
    ~Test() { delete[] num; }
};
/******************************数据结构测试**************************************/
class Data {
private:
    int row;
    int column;
    double* data;
public:
    Data(int r, int c);
    double& operator()(int r, int c);
    ~Data();
};
Data::Data(int r, int c) :row(r), column(c), data(nullptr)
{
    int size = r * c;
    data = new double[size]();
}
double& Data::operator()(int r, int c)
{
    return data[r*column + c];
}
Data::~Data()
{
    delete[] data;
}
//////////////////////////////////////////////////////////////////////////////////////////
class Data2 {
private:
    int row;
    int column;
    int* Row;
    double* data;
public:
    Data2(int r, int c);
    double& operator()(int r, int c);
    ~Data2();
};
Data2::Data2(int r, int c) :row(r), column(c), data(nullptr)
{
    int size = r * c;
    data = new double[size]();
    Row = new int[r];
    for (int i = 0; i < row; i++)
    {
        Row[i] = i * column;
    }
}
double& Data2::operator()(int r, int c)
{
    return data[Row[r] + c];
}
Data2::~Data2()
{
    delete[] data;
    delete[] Row;
}
/*****************************************************************************/
void print(const Matrix& m)
{
    int r = std::get<0>(m.size());
    int c = std::get<1>(m.size());
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            cout << m[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

int main()
{
    // 5X5矩阵, rank = 5
    Matrix A = { {1,10,3,6,7}, {-100,5,4,1,9}, {2,0,3,1,1}, {4,4,9,2,0}, {1,5,2,3,4} };
    // 5X5矩阵, rank = 4
    Matrix B = { {1,10,3,6,7}, {1,10,3,6,7}, {2,0,3,1,1}, {4,4,9,2,0}, {1,5,2,3,4} };
    // 4X5矩阵, rank = 4
    Matrix C = { {1,10,3,6,7}, {-100,5,4,1,9}, {2,0,3,1,1}, {4,4,9,2,0} };
    // 5X4矩阵, rank = 4
    Matrix D = { {1,10,3,6}, {-100,5,4,1}, {2,0,3,1}, {4,4,9,2}, {1,5,2,3} };

    // 向量
    Matrix b = Matrix({ 1,2,3,4,5 }).trans();
    /****************************************************************************************/

    //// 矩阵的秩
    //cout << "rank(A) = " << A.rank() << endl;

    //// 矩阵的迹
    //cout << "trace(A) = " << A.trace() << endl;

    //// 行列式
    //cout << "det(A) = " << det(A) << endl;

    //// 逆矩阵
    //cout << "inv(A) = " << endl;
    //print(inv(A));

    //// QR分解
    //cout << "QR Decomposition:" << endl;
    //auto QRVec = QR(A);
    //cout << "Q = " << endl;
    //print(QRVec[0]);
    //cout << "R = " << endl;
    //print(QRVec[1]);

    //// 解线性方程组
    //cout << "方程组的解：" << endl;
    //print(std::get<0>(solve(A, b)));

    //// 特征值
    //auto val = eigs(A);
    //for (auto& v : val) { cout << v << endl; }

    system("pause");
    return 0;
}