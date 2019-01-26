# API

## 1.简介

* 包含矩阵的常用操作，如四则运算、求逆矩阵、求特征值、解线性方程组等。
* 编译器需支持C++11及以上标准。

## 2.矩阵构造

* **构造函数（矩阵元素为double类型）**

```cpp
Matrix();                                                                // 空矩阵
Matrix(const size_t& m, const size_t& n);                                // 常规构造: m行数，n列数
Matrix(const std::initializer_list<std::initializer_list<double>>& m);   // 列表构造
Matrix(const Matrix& other);                                             // 拷贝构造
Matrix(Matrix&& other);                                                  // 移动构造
```

**Examples:**

```cpp
Matrix A;
Matrix B(2,3);
Matrix C1({{1,2,3}});
Matrix C2({{1},{2},{3}});
Matrix D({1,2,3},{4,5,6});
Matrix E(D);
```

**Results：**
$$
A = (),
\quad
B = \left(
\begin{matrix}
   0 & 0 & 0 \\
   0 & 0 & 0
  \end{matrix}
\right),
\quad
C1 = \left(
\begin{matrix}
   1 & 2 & 3
  \end{matrix}
\right),
\quad
C2 = \left(
\begin{matrix}
   1 \\
   2 \\
   3
  \end{matrix}
\right),
\quad
D = \left(
\begin{matrix}
   1 & 2 & 3 \\
   4 & 5 & 6
  \end{matrix}
\right),
\quad
E = \left(
\begin{matrix}
   1 & 2 & 3 \\
   4 & 5 & 6
  \end{matrix}
\right)
$$

## 3.元素索引

支持`[]`操作符和`at()`函数进行矩阵元素索引。

**Example:**

```cpp
Matrix M({1,2,3},{4,5,6});
M[1][1] = 0;                // []操作符
M.at(0,2) = 0.5;            // at()函数
```

**Result：** $M = \left( \begin{matrix}  1 & 2 & 0.5 \\  4 & 0 & 6 \end{matrix} \right)$

**`[]`操作符与`at()`函数的区别:** `[]`操作符不进行下标越界检查，`at()`函数会进行下标检查，若下标越界则抛出`std::out_of_range`异常。

## 4.矩阵生成

* `eye`
* `ones`
* `zeros`
* `hilb`
* `rand`
* `randn`
* `diag`
* `getDiag`
* `rbind`
* `cbind`
* `subMat`

## 5.四则运算

* **赋值运算**
  =, +=, -=, *=, /=

* **矩阵加法**

* **矩阵减法**

* **矩阵乘法**

* **矩阵除法**

## 6.矩阵属性

* `size`

**功能：** 获取矩阵行数和列数
**函数原型：** `const std::pair<usize, usize> size(const Matrix& mat)`
**输入：** 矩阵行数
**返回：** 矩阵行数和列数
**示例：** `auto a = MAT.size();` or `auto a = size(MAT);`

* `trans`
* `rank`
* `trace`
* `inv`
* `kernel`
* `det`
* `normOne, normTwo, normInf`

## 7.矩阵分解

* **QR分解**

* **LU分解(Doolittle分解)**

* **特征值分解**

* **SVD分解**

## 8.其它函数

* `getRow`

**功能：** 获取矩阵**行**向量
**函数原型：** `Matrix getRow(usize n)`
**输入：** 矩阵行数(≥0)
**返回：** 行向量
**示例：** `Matrix vec = MAT.getRow(0);` or `Matrix vec = getRow(MAT, 0);`

* `getColumn`

**功能：** 获取矩阵**列**向量
**函数原型：** `Matrix getRow(usize n)`
**输入：** 矩阵列数(≥0)
**返回：** 列向量
**示例：** `Matrix vec = MAT.getColumn(0);` or `Matrix vec = getColumn(MAT, 0);`

* `rref`
* `isZero`
* `filter`
* `solve`
* `map`