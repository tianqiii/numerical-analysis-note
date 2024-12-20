# LU分解

求解$Ax=b$时，可以先将矩阵A分解为LU，原式变为$LUx=b$再另$Ux=y$将原式变换为$Ly=b$

先求$y$再求$x$ 因为 ***doolittle*** 分解默认$L$的对角线元素均为1

求解顺序：$U$第一行，$L$第一列，$U$第二行，$L$第二列，$U$第三行……

**例：**

$$
A=\begin{bmatrix}1&2&3&4\\1&3&4&5\\2&1&4&4\\2&3&2&5\end{bmatrix}
x={\begin{bmatrix}x_1\ x_2\ x_3\ x_4\end{bmatrix}}^T
b={\begin{bmatrix}2\ 3\ 3\ 2\end{bmatrix}}^T
$$

$$
A=LU
$$

$$
L=\begin{bmatrix}1&0&0&0\\1&1&0&0\\2&-3&1&0\\2&-1&-3&1 \end{bmatrix}
U=\begin{bmatrix}1&2&3&4\\0&1&1&1\\0&0&1&-1\\0&0&0&-5 \end{bmatrix}
$$

$$
Ly=b \rightarrow y={\begin{bmatrix}2&1&2&5\end{bmatrix}}^T
$$

$$
Ux=y \rightarrow x={\begin{bmatrix}1&1&1&-1\end{bmatrix} }^T
$$

# Hermite插值

### 方法：

1. 写差商表（重复的几次导数就写几次）

2. 差商乘以节点

3. 整理结果

4. 写截断误差

**例：**

构造三次多项式$p_3(x)$，使曲线$y=p_3(x)$满足如下表格所示的插值条件，求出插值多项式$p_3(x)$的截断误差（不需证明）。

| $x$     | 0   | 1   | 3   |
|:-------:|:---:|:---:|:---:|
| $f(x)$  | 1   | 0   | -2  |
| $f'(x)$ |     | 1   |     |

写差商表如下：

| 0   | 1   | -1  | 2   | -1  |
|:---:|:---:|:---:|:---:| --- |
| 1   | 0   | 1   | -1  |     |
| 1   | 0   | -1  |     |     |
| 3   | -2  |     |     |     |

$$
P_3(x)=1+(-1)(x-0)+2(x-0)(x-1)+(-1)(x-0)(x-1)^2
$$

$$
R_3(x)=\frac{f^{(4)}(\xi)}{4!}x(x-1)^2(x-3)
$$

# 最小二乘法

### 方法

1. 整理公式

2. 写出矛盾方程组

3. 转置左乘

4. 解出未知变量

5. 代回原式

**例：**

用最小二乘法求解一个形如$y=a+bx^2$的经验公式，使它与下列的数据拟合，计算误差平方和。

| $x_i$   | -1.0 | 0    | 1.0  | 2    |
|:-------:|:----:|:----:|:----:|:----:|
| $x_i^2$ | 1    | 0    | 1    | 4    |
| $y_i$   | 1.91 | 1.05 | 2.08 | 5.21 |

$$
a+bx^2=y
$$

$$
\begin{bmatrix}1&1\\1&0\\1&1\\1&4 \end{bmatrix}
\begin{bmatrix}a\\b \end{bmatrix}=
\begin{bmatrix}1.91\\1.05\\2.08\\5.21 \end{bmatrix}
$$

转置左乘，有：

$$
\begin{bmatrix}1&1&1&1\\1&0&1&4 \end{bmatrix}
\begin{bmatrix}1&1\\1&0\\1&1\\1&4 \end{bmatrix}
\begin{bmatrix}a\\b \end{bmatrix}=
\begin{bmatrix}1&1&1&1\\1&0&1&4 \end{bmatrix}
\begin{bmatrix}1.91\\1.05\\2.08\\5.21 \end{bmatrix}
$$

$$
\begin{bmatrix}4&6\\6&18 \end{bmatrix}
\begin{bmatrix}a\\b \end{bmatrix}=
\begin{bmatrix}10.25\\24.83 \end{bmatrix}
$$

解得$a=0.9867,\ b=1.0506$，$y=0.9867+1.0506x^2$

# 不动点迭代和牛顿迭代法

### 不动点迭代法的全局收敛定理：

设方程$x=\varphi(x)$满足

1. 迭代函数$\varphi(x)$在区间 **[a,b]** 上导函数存在且连续；

2. 当$x\in[a,b]$时，$\varphi(x)\in[a,b]$；（封闭性）

3. 存在常数 $0<L<1$，使对任意$x\in(a,b)$有$|\varphi'(x)|\leq L$.（压缩性）

**不动点迭代格式：**

$$
x_{k+1}=\varphi(x_k)
$$

**收敛速度:**

$$
s=\lim_{i\to\infty}\frac{e_{i+1}}{e_i},\ e_i=|x_i-r|
$$

其中r是该不动点格式的精确解。

### Newton迭代法的全局收敛性

设有$f(x)=0$且$f(x)\in C^2[a,b]$，（二阶导数存在）且满足：

1. $f(a)f(b)<0$

2. $f'(x)\neq 0$，$f''(x)\neq0$，$\forall x\in[a,b]$（$f'(x),f''(x)$不变号）

3. 取$x_0\in[a,b]$，使$f(x_0)f''(x_0)>0$
   
   则（1）$f(x)=0$在$[a,b]$内有唯一的根$x^*$
   
   （2）由Newton迭代式产生的序列${x_k}$收敛于$x^*$，且
   
   $$
   \lim_{k\to\infty}\frac{x_{k+1}-x^*}{(x_k-x^*)^2}=\frac{f^*(x^*)}{2f'(x^*)}
   $$

**Newton迭代格式：**

$$
x_{k+1}=x_k-\frac{f(x_k)}{f'(x_k)},k=0,1,2,\cdots
$$

**例：**

设有迭代格式$x_{n+1}=\frac{3}{4}+\frac{1}{3}x_n^2$

1. 若该格式收敛，试求所收敛到的收敛值；

2. 在该极限值附近，试确定此格式的局部收敛阶是几阶。

**解：**

1. 两边取极限：$x=\frac{3}{4}+\frac{1}{3}x$，解得极限值$x=1.5$.

2. $$
   \frac{x_{n+1}-1.5}{(x_n-1.5)^p}=
\frac{\frac{1}{3}x_n^2-0.75} {(x_n-1.5)^p}=
\frac{x_n^2-2.25}{3(x_n-1.5)^p}=
\frac{(x_n-1.5)(x_n+1.5)}{3(x_n-1.5)^p}
   $$

当$p=1$时，原式$=\frac{x_n+1.5}{3}$ 取极限等于1，所以$p=1$，一阶收敛.

# 误差

1. **绝对误差**
   
   **定义：**
   
   $x^*-x=e(x^*)$
   
   在不引起混淆时，简记$e(x^*)$为$e^*$.
   
   **绝对误差限：**
   
       如果存在正数$\varepsilon^*=\varepsilon(x^*)$，使得有绝对误差$|e^*|-|x^*-x|\leq\varepsilon^*$，
   
   则称$\varepsilon^*$为$x^*$近似$x$的一个**绝对误差限**。$x\in[x^*-\varepsilon^*,x^*+\varepsilon^*],\ x=x^*\pm\varepsilon^*.$
   
   通常计算中所要求的误差，是指估计尽可能小的绝对误差限。

2. **相对误差**
   
   **定义：**
   
   设$x^*$是对准确值$x(\neq0)$的一个近似，则称
   
   $$
   e_r(x^*)=\frac{x^*-x}{x}=\frac{e(x^*)}{x}
   $$
   
   为$x^*$近似$x$的相对误差，不引起混淆时，简记$e_r(x^*)$为$e_r^*$

3. **有效数字**
   
   **定义：**
   
   设$x$的近似值$x^*$有如下标准形式
   
   $$
   x^*=\pm10^m \times 0.\underline{x_1 x_2 x_3 \cdots x_n x_{n+1} 
\cdots x_p}
   $$
   
   其中$m$为整数，$\{x_i\}\subset\{0,1,2,\cdots,9\}$且$x_i\neq0，p\geq n.$如果有
   
   $$
   |e^*|=|x-x^*|\leq\frac{1}{2}\times10^{m-n}
   $$
   
   则称$x^*$为$x$的具有$n$位有效数字的近似数，或称$x^*$准确到$10^{m-n}$位，其中数字$x_1,x_2,\cdots,x_n$分别称为$x^*$的第$1,2,\cdots,n$个有效数字。
   
   **例：**
   
   $$
   x=\pi,x_1^*=3.141,x_2^*=3.142
   $$
   
   $$
   |x-x_1^*|=0.00059\cdots\leq0.005=\frac{1}{2} \cdot 10^{1-3}
   $$
   
   3位有效数字，非有效数
   
   $$
   |x-x_2^*|=0.00040\cdots\leq0.0005=\frac{1}{2}\cdot10^{1-4}
   $$
   
   4位有效数字，有效数

4. **二元函数算术运算误差传播规律**
   
   **绝对误差限：**
   
   $$
   \varepsilon(x_1^*\pm x_2^*)\approx
\varepsilon(x_1^*)+\varepsilon(x_2^*)
   $$
   
   $$
   \varepsilon(x_1^* x_2^*)\approx
|x_2^*|\varepsilon(x_1^*)+|x_1^*|\varepsilon(x_2^*)
   $$
   
   $$
   \varepsilon(\frac{x_1^*} {x_2^*})\approx
\frac{|x_2^*|\varepsilon(x_1^*)+|x_1^*|\varepsilon(x_2^*)}
{|x_2^*|^2}\quad(x_2^*\neq0)
   $$
   
   **相对误差限：（可能要考）**
   
   $$
   \varepsilon_r(x_1^*+x_2^*)\approx max\{\varepsilon_r(x_1^*),
\varepsilon_r(x_2^*)\}\quad(x_1^*x_2^*>0)
   $$
   
   $$
   \varepsilon_r(x_1^*x_2^*)\approx \varepsilon_r(x_1^*)+
\varepsilon_r(x_2^*)\quad(x_1^*x_2^*\neq0)
   $$
   
   $$
   \varepsilon_r(\frac{x_1^*} {x_2^*})\approx
\varepsilon_r(x_1^*)+\varepsilon_r(x_2^*)
\quad(x_1^*x_2^*\neq0)
   $$

**例1：**

设近似数$x_1^*=2.1567$是有效数，$x_2^*=-3.5017$是由准确数四舍五入所得的数，则相对误差$|e_r(x_1^*x_2^*)|\leq$ ？

$$
|e_r(x_1^*x_2^*)|=|e_r(x_1^*)|+|e_r(x_2^*)\leq
|\frac{0.5\times10^{1-5}} {2.1567}|+
|\frac{0.5\times10^{1-5}}{-3.5017}|\approx3.75\times10^{-5}
$$

**例2:**

近似数$x^*=2.453800$具有7位有效数字？

近似数代表最后一位数字是近似出来的，所以有 $0.000005=\frac{1}{2}\times10^{1-6}$

所以是6位有效数字

# 迭代法解线性方程组

# 数值积分

# 拉格朗日插值
