#include <iostream>
#include <windows.h>
#include <string>
#include <memory>
#include <cmath>
#include <complex>
#include <vector>
#include "matrix.h"
using namespace std;

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
void print(Matrix& m)
{
    int r = m.size()[0];
    int c = m.size()[1];
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            cout << m(i, j) << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

int main()
{
    //Matrix m = { { 1,2,3,4,1 },{ 7,6,1,5,2 },{ 3,1,1,1,3 },{ 4,6,9,2,4 },{ 7,9,0,1,3 } };
    //Matrix m2 = { { 1,2,3,4,1 },{ 7,6,1,5,2 },{ 3,1,1,1,3 },{ 4,6,9,2,4 },{ 2,4,6,8,2 } };
    //Matrix m3 = Matrix::diag({ 1,2,3,4,5 });
    ////DWORD start = GetTickCount();
    ////DWORD end = GetTickCount();
    ////cout << "The run time is:" << (end - start) << " ms" << endl;
    //m3(4, 3) = 2;
    //m3(3, 4) = 1;
    //auto a = m.eigs();

    //Matrix matrix = { {1,10,3,6,7},{-100,5,4,1,9},{2,0,3,1,1},{4,4,9,2,0},{1,5,2,3,4} };
    //Matrix b = Matrix({ 1,2,3,4,5 }).trans();
    //Matrix matrix = { {3,4},{6,7},{6,8} };
    //Matrix b = Matrix({ 1,3,2 }).trans();

    //Matrix matrix = { { 3,4,5 } };
    //Matrix b = Matrix({ 7 }).trans();

    //cout << "The matrix is:" << endl;
    //print(matrix);
    //cout << endl;
    //cout << "The b is:" << endl;
    //print(b);
    //cout << endl;
    //cout << "The solution is:" << endl;
    //print(get<0>(Matrix::solve(matrix, b)));
    //cout << get<1>(Matrix::solve(matrix, b));
    //cout << endl;

    //cout << "Rank = " << matrix.rank() << endl << endl;
    //cout << "After gauss elimination:" << endl;
    //auto temp = Matrix::gaussElimination(matrix);
    //print(get<0>(temp));
    //cout << endl;
    Matrix m = { { 1,2 },{ 3,4 } };
    cout << m.det() << endl;

    system("pause");
}