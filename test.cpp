#include <iostream>
#include <windows.h>
#include <string>
#include <memory>
#include <cmath>
#include <complex>
#include <vector>
#include "matrix\matrix.h"

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
    //Matrix B = { {1,10,3,6,7}, {1,10,3,6,7}, {2,0,3,1,1}, {4,4,9,2,0}, {1,5,2,3,4} };
    //Matrix B = { {1,10,3,6,7}, {1,10,3,6,7}, {2,0,3,1,1}, {1,10,3,6,7}, {1,10,3,6,7} };
    Matrix B = { {1,10,3,6,7}, {1,10,3,6,7}, {2,0,3,1,1} };
    Matrix B2 = { {11, 11, 11, 11, 11} };
    Matrix B3 = { { 22, 22, 22, 22, 22 } };
    Matrix B4 = { { 33, 33, 33, 33, 33 } };

    // 4X5矩阵, rank = 4
    Matrix C = { {1,10,3,6,7}, {-100,5,4,1,9}, {2,0,3,1,1}, {4,4,9,2,0} };
    // 5X4矩阵, rank = 4
    Matrix D = { {1,10,3,6}, {-100,5,4,1}, {2,0,3,1}, {4,4,9,2}, {1,5,2,3} };

    // 向量
    Matrix b = Matrix({ {1,2,3,4,5} }).trans();
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
    //print(get<0>(QRVec));
    //cout << "R = " << endl;
    //print(get<1>(QRVec));

    // LU分解
    cout << "LU Decomposition:" << endl;
    auto QRVec = LU(B);
    cout << "P = " << endl;
    print(get<0>(QRVec));
    //cout << "L = " << endl;
    //print(get<1>(QRVec));
    //cout << "U = " << endl;
    //print(get<2>(QRVec));
    //cout << "A = " << endl;
    //print(B);
    //cout << "-------------------------------------------------" << endl;
    //print(get<0>(QRVec) * get<1>(QRVec) * get<2>(QRVec));

    //// 解线性方程组
    //cout << "方程组的解：" << endl;
    //print(std::get<0>(solve(A, b)));

    //// 特征值
    //cout << "矩阵A的特征值：" << endl;
    //auto val = eigVal(A);
    //for (auto& v : val) { cout << v << endl; }

    //// 特征值和特征向量
    //cout << "矩阵A的特征值：" << endl;
    //auto val = A.eig();
    //for (auto& v : val)
    //{
    //    cout << "特征值：" << get<0>(v) << ",  特征向量：";
    //    for (int i = 0; i < get<0>(get<1>(v).size()); i++)
    //    {
    //        cout << get<1>(v)[i][0] << "  ";
    //    }
    //    cout << endl;
    //}

    //// 矩阵拼接
    //print(A);
    //cout << "-------------------------------------------------" << endl;
    ////print(A.cbind(B2.trans(), B3.trans(), B4.trans()));
    //print(cbind(A, B2.trans(), B3.trans(), B4.trans()));
    //cout << "-------------------------------------------------" << endl;
    ////print(A.rbind(B2, B3, B4));
    //print(rbind(A, B2, B3, B4));

    //// 生成对角阵
    //print(diag(B2));

    // 矩阵的核
    //Matrix E = { {1,2,2,2},{2,4,6,8},{3,6,8,10} };
    //print(E.kernel());

    system("pause");
    return 0;
}