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

* `eye`
* `ones`
* `zeros`
* `rand`
* `diag`
* `subMat`

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

* `size`
**功能：** 获取矩阵行数和列数
**函数原型：** `const std::pair<usize, usize> size(const Matrix& mat)`
**输入：** 矩阵行数
**返回：** 行向量
**示例：** `auto a = MAT.size();` or `auto a = size(MAT);`

* `getRow`
**功能：** 获取矩阵行向量
**函数原型：** `Matrix getRow(usize n)`
**输入：** 矩阵行数
**返回：** 行向量
**示例：** `Matrix vec = MAT.getRow(0);` or `Matrix vec = getRow(MAT, 0);`

* `Matrix getColumn`
* `getDiag`
* `rbind`
* `cbind`
* `subMat`
* `trans`
* `rank`
* `trace`
* `inv`
* `det`
* `normOne, normTwo, normInf`
* `eig`
* `solve`
* `filter`
* `map`

### 矩阵分解

* **QR分解**
* **LU分解(Doolittle分解)**
* **SVD分解**
