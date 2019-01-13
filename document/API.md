# API

### 矩阵构造

* **1. 构造函数（矩阵元素为double类型）**

```cpp
Matrix();                                                                // 空矩阵
Matrix(const size_t& m, const size_t& n);                                // 常规构造: m行数，n列数
Matrix(const std::initializer_list<double>& m);                          // 列表构造
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

### 常用构造

* `Matrix eye(const usize& m, const usize& n)`
* `Matrix ones(const usize& m, const usize& n)`
* `Matrix zeros(const usize& m, const usize& n)`
* `Matrix rand(const usize& m, const usize& n)`
* `Matrix diag(const Matrix& vec)`
* `Matrix subMat()`

### 元素索引

```cpp
Matrix M({1,2,3},{4,5,6});
M[1][1] = 0;                // []操作符
M.at(0,2) = 0.5;            // at()函数
```

**Result：**$
M = \left(
\begin{matrix}
   1 & 2 & 0.5 \\
   4 & 0 & 6 
  \end{matrix}
\right)
$

**`[]`操作符与`at()`函数的区别:** `[]`操作符不进行下标越界检查，`at()`函数会进行下标检查，若下标越界则抛出`std::out_of_range`异常。

### 四则运算

* **赋值运算**

* **加法**

* **减法**

* **乘法**

* **除法**

### 矩阵基本操作

* `size()`
* `Matrix getRow(usize n)`
* `Matrix getColumn(usize n)`
* `std::vector<double> getDiag()`
* `Matrix rbind(Matrix& M1, [M2, M3,...])`
* `Matrix cbind(Matrix& M1, [M2, M3,...])`
* `Matrix subMat(usize r1, usize c1, usize r2, usize c2)`
* `Matrix trans()`
* `size_t rank()`
* `double trace()`
* `Matrix inv()`
* `double det()`
* `double normOne(), double normTwo(), double normInf()`
* `std::vector<std::complex<double>> eig()`
* `solve()`
* `Matrix filter()`
* `Matrix map(std::function<double(double)> f)`

### 矩阵分解

* `std::pair<Matrix, Matrix> QR()`
* `std::tuple<Matrix, Matrix, Matrix> LU()`
* `std::tuple<Matrix, Matrix, Matrix> SVD()`
