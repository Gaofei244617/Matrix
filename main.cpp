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

// 输出矩阵
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

    // 矩阵的秩
    cout << "rank(A) = " << A.rank() << endl;

    // 矩阵的迹
    cout << "trace(A) = " << A.trace() << endl;

    // 行列式
    cout << "det(A) = " << det(A) << endl;

    // 逆矩阵
    cout << "inv(A) = " << endl;
    print(inv(A));

    // QR分解
    cout << "QR Decomposition:" << endl;
    auto QRVec = QR(A);
    cout << "Q = " << endl;
    print(get<0>(QRVec));
    cout << "R = " << endl;
    print(get<1>(QRVec));

    // 解线性方程组
    cout << "方程组的解：" << endl;
    print(std::get<0>(solve(A, b)));

    // 特征值
    auto val = eigs(A);
    for (auto& v : val) { cout << v << endl; }

    system("pause");
    return 0;
}