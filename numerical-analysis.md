# LU分解

求解$Ax=b$时，可以先将矩阵$A$分解为$LU$，原式变为$LUx=b$再另$Ux=y$将原式变换为$Ly=b$

先求$y$再求$x$ 因为***doolittle***分解默认$L$的对角线元素均为1

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

### 方法:

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

**不动点迭代法的误差估计：**

$$\begin{aligned}|x^*-x_k|&\leq\frac{L}{1-L}|x_k-x_{k-1}|\\
&\leq\frac{L^k}{1-L}|x_1-x_0|\end{aligned}$$
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

**Newton迭代法的误差估计：**

$$
x_{k+1}-x^*=\frac{f''(\xi)}{f'(x_k)}(x_k-x^*)^2
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
   $$\begin{aligned}
   &\varepsilon(x_1^*\pm x_2^*)\approx
\varepsilon(x_1^*)+\varepsilon(x_2^*)
   \quad\\
   &\varepsilon(x_1^* x_2^*)\approx
|x_2^*|\varepsilon(x_1^*)+|x_1^*|\varepsilon(x_2^*)
   \quad\\
   &\varepsilon(\frac{x_1^*} {x_2^*})\approx
\frac{|x_2^*|\varepsilon(x_1^*)+|x_1^*|\varepsilon(x_2^*)}
{|x_2^*|^2}\quad(x_2^*\neq0)
   \end{aligned}$$
   **相对误差限：**
   
   $$\begin{aligned}
   &\varepsilon_r(x_1^*+x_2^*)\approx max\{\varepsilon_r(x_1^*),
\varepsilon_r(x_2^*)\}\quad&(x_1^*x_2^*>0)
   \\
   &\varepsilon_r(x_1^*x_2^*)\approx \varepsilon_r(x_1^*)+
\varepsilon_r(x_2^*)\quad&(x_1^*x_2^*\neq0)
   \\
   &\varepsilon_r(\frac{x_1^*} {x_2^*})\approx
\varepsilon_r(x_1^*)+\varepsilon_r(x_2^*)
\quad&(x_1^*x_2^*\neq0)
   \end{aligned}$$

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

### 简单迭代法的收敛性（适用于所有迭代法）

**定理：** 简单迭代法$\overline{x}^{(k+1)}=B\overline{x}^{(k)}+\overline{g},\quad k=0,1,2,\cdots$对任意初始向量$\overline{x}^{(0)}$都收敛的充要条件是：**迭代矩阵的谱半径$\rho(B)<1\ or\ \lim_{k\to\infty}B^k=0$** 

**$\rho(B)<1$是判定收敛的根本法则；**

**当$\rho(B)\geq1$时，有可能存在某个初始向量$\overline{y}^{(0)}$使简单迭代法收敛。**

$$
A=\begin{bmatrix}
a_{11}&a_{12}&a_{13}\\
a_{21}&a_{22}&a_{23}\\
a_{31}&a_{32}&a_{33}\\
\end{bmatrix},x=\begin{bmatrix}x_1\\ x_2\\ x_3\end{bmatrix},
b=\begin{bmatrix}b_1\\ b_2\\ b_3\end{bmatrix}.
$$

$$
Ax=b\\
A=D+L+U\\
(D+U+L)x=b\\

$$

### 1. Jacobi迭代法

$$
Dx_{k+1}=-(L+U)x_k+b\\
x_{k+1}=-D^{-1}(L+U)x_k+D^{-1}b
$$

能够写成如下格式：

$$\begin{aligned}
x_1^{(k+1)}=\frac{1}{a_{11}}(b_1-a_{12}x_2^{(k)}-a_{13}x_3^{(k)})\\
\quad \\
x_2^{(k+1)}=\frac{1}{a_{22}}(b_2-a_{21}x_1^{(k)}-a_{23}x_3^{(k)})\\
\quad \\
x_3^{(k+1)}=\frac{1}{a_{33}}(b_3-a_{31}x_1^{(k)}-a_{32}x_2^{(k)})\\
\end{aligned}$$

### 2.Gauss-Seidel迭代法

$$
(D+L)x_{k+1}=-Ux_k+b\\
x_{k+1}=-(D+L)^{-1}Ux_k+(D+L)^{-1}b
$$

能够写成如下格式：

$$\begin{aligned}
x_1^{(k+1)}&=\frac{1}{a_{11}}(b_1-a_{12}x_2^{(k)}-a_{13}x_3^{(k)})\\
\quad \\
x_2^{(k+1)}&=\frac{1}{a_{22}}(b_2-a_{21}x_1^{(k+1)}-a_{23}x_3^{(k)})\\
\quad \\
x_3^{(k+1)}&=\frac{1}{a_{33}}(b_3-a_{31}x_1^{(k+1)}-a_{32}x_2^{(k+1)})\\
\end{aligned}$$

**Gauss-Seidel迭代法与Jacobi迭代法的区别：**

GS迭代法逐次更新变量，Jacobi迭代法完成一轮才更新变量。

### 3.Sor迭代法

超松弛迭代法是由Gauss-Seidel迭代法改写得到的：

$$\begin{aligned}
&(D+L)x_{k+1}=-Ux_k+b\\
\quad \\
&x_{k+1}=(1-\omega)x_k+\omega[-(D+L)^{-1}Ux_k+(D+L)^{-1}b]\\
\quad\\
& \hspace{5cm} \text{or} \hspace{2cm} \\
\quad\\
&Dx_{k+1}+Lx_{k+1}=-Ux_k+b\\
\quad \\
&Dx_{k+1}=(1-\omega)Dx_k+\omega(b-Lx_{k+1}-Ux_{k})
\end{aligned}$$

即：

$$
x_{k+1}=(D+\omega L)^{-1}((1-\omega)D-\omega U)x_k+
\omega(D+\omega L)^{-1}b
$$

|         | Jacobi | Gauss-Seidel | Sor                |
|:-------:|:------:|:------------:|:------------------:|
| A严格对角占优 | 收敛     | 收敛           | $0<\omega\leq1$时收敛 |
| A对称正定   | 未必收敛   | 收敛           | $0<\omega<2$时收敛    |

严格对角占优：$|a_{ii}|>\sum_{j\neq i}|a_{ij}|$

对称正定：$A_{ij}=A_{ji}$，$\lambda_i>0,i=1,2,\dots,n$

# 数值积分

### 求积公式

$$
I=\int_a^bf(x)dx\approx\sum_{k=0}^nA_kf(x_k)+R[f]=I_n+R[f]
$$

其中$R[f]$称为**求积公式的余项，**$x_k(k=1,2,\cdots,n)$称为**求积节点，**$A_k(k=0,1,2,\cdots,n)$称为**求积系数。**$A_k$仅与求积节点$x_k$的选取有关，而不依赖与被积函数$f(x)$的具体形式。

### 代数精度

衡量一个求积公式好坏的标准。

**定义：** 如果求积公式$\int_a^bf(x)dx\approx\sum_{k=0}^nA_kf(x_k)$对于一切不高于$m$次的代数多项式准确成立，而对于某个$m+1$次多项式并不准确成立，则称上述求积公式具有 **$m$次代数精确度，** 或称为具有 **$m$次代数精度。** 

1. 代数精度越高，求积公式的适应性越强。

2. 凡至少具有零次代数精度的求积公式
   
   $$
   \int_a^bf(x)dx\approx\sum_{k=0}^nA_kf(x_k)
   $$
   
   一定满足
   
   $$
   \int_a^b1dx=\sum_{k=0}^nA_k \cdot 1
   $$
   
   从而有
   
   $$
   \sum_{k=0}^nA_k=b-a
   $$
   
   即求积系数之和等于积分区间长度，这是求积系数的一个**基本特性**。

**例：**

设用函数值$f(0),f(1)$和导数值$f'(0),f'(1)$的线性组合近似表达积分$\int_0^1f(x)dx$，试求：(1) 求四个组合系数$A_0,A_1,B_0,B_1$的值

(2) 求该数值求积公式的截断误差$E(f)$.

1. $$
   \int_0^1f(x)dx\approx A_0f(0)+A_1f(1)+B_0f'(0)+B_1f'(1)
   $$
   
   另$f(x)$分别等于$1;x;x^2;x^3$，可以解得
   
   $$
   A_0=\frac{1}{2},A_1=\frac{1}{2};
B_0=\frac{1}{12},B_1=-\frac{1}{12}.
   $$
   
   可以知道其有3次代数精度，因为当$f(x)=x^4$时，该求积公式不再准确。

2. $$\begin{aligned}
   E(f)&=\int_0^1f(x)dx-[A_0f(0)+A_1f(1)+B_0f'(0)+B_1f'(1)]\\
   &=\int_0^1f(x)dx-[A_0H_3(0)+A_1H_3(1)+B_0H_3'(0)+B_1H_3'(1)]\;(使用插值条件)\\
   &=\int_0^1f(x)dx-\int_0^1H_3(x)dx=\int_0^1(f-H_3)dx\;(三次代数精度)\\
   &=\int_0^1\frac{f^{(4)}(\xi)}{4!}x^2(x-1)^2dx=\frac{1}{720}f^{(4)}(\eta)\;(使用积分中值定理)
   \end{aligned}$$

### 收敛性与稳定性

如果$\lim_{n\to\infty}\sum_{k=0}^nA_kf(x_k)=\int_a^bf(x)dx,\,(\lim_{h\to\infty\\ n\to\infty}R[f]=0)$其中$h=\max_{1\leq i\leq n}(x_i-x_{i-1})$，则称该求积公式是**收敛的。** 如果求积公式对舍入误差不敏感（误差能够控制），则称该求积公式是**稳定的。**

一个求积公式首先应该是收敛的，其次应该是稳定的。

当系数$A_k$全为正时，求积公式是稳定的。

### Newton-Cotes公式

节点等间距分布的插值型求积公式即为**Newton-Cotes求积公式。** 当$n\leq7$时，N-C公式是稳定的；$n\geq7$时则是不稳定的。

对于$n+1$个节点的N-C公式，当$n$为奇数时，其代数精度至少为$n$次；n为偶数时，其代数精度至少为$n+1$次。

该求积公式可以写为

$$
\int_a^bf(x)dx\approx(b-a)\sum_{k=0}^nc_k^{(n)}f(x_0+kh)
$$

$$
c_k^{(n)}=\frac{(-1)^{n-k}}{nk!(n-k)!}\int_0^nt(t-1)\cdots(t-k+1)(t-k)
\cdots(t-n)dt
$$

其中$c_k^{(n)}$称为**柯特斯系数，** 上式称为**Newton-Cotes公式。** 可以证明$c_k^{(n)}=c_{n-k}^{(n)}$。

### 梯形公式（n=1）

$$
\int_a^bf(x)dx\approx\frac{(b-a)}{2}[f(a)+f(b)]
$$

上式称为**梯形公式。**

### Simpson公式（n=2）

$$
\int_a^bf(x)dx\approx(b-a)[\frac{1}{6}f(a)+\frac{4}{6}f(\frac{a+b}{2})+
\frac{1}{6}f(b)]
$$

上式称为**Simpson公式或抛物线公式。**

### Cotes公式（n=4）

$$
\int_a^bf(x)dx\approx(b-a)[\frac{7}{90}f(x_0)+\frac{32}{90}f(x_1)+
\frac{12}{90}f(x_2)+\frac{32}{90}f(x_3)+\frac{7}{90}f(x_4)]
$$

其中$x_k=a_0+kh\;(k=0,1,2,3,4)$，这个公式特别称为**柯特斯公式。**

### 复化梯形公式

$$
T(n)=\frac{h}{2}[f(a)+2\sum_{k=1}^{n-1}f(x_k)+f(b)]
$$

截断误差：

$$
R_{T_{(n)}}[f]=I-T(n)=-\frac{h^3}{12}\sum_{k=0}^{n-1}f^n(\eta_k)
=-\frac{b-a}{12}h^2f''(\eta)
$$

若取9个节点，其步长以及与9个节点所对应的求积系数分别是多少？

复化梯形公式：$n=8,h=\frac{b-a}{8}$，对应的求积系数为：

1，2，2，2，2，2，2，2，1。

### 复化Simpson公式

$$
S(n)=\frac{h}{6}[f(a)+4\sum_{k=0}^{n-1}f(x_{k+\frac{1}{2}})+
2\sum_{k=1}^{n-1}f(x_k)+f(b)]
$$

截断误差：

$$
R_{S(n)}[f]=I-S(n)=-\frac{b-a}{2880}h^4f^{(4)}(\eta)
$$

若取9个节点，其步长以及与9个节点所对应的求积系数分别是多少？

复化Simpson公式：$n=4,h=\frac{b-a}{4}$，对应的求积系数为：

1，4，2，4，2，4，2，4，1。

**复化公式的收敛性**：当$h\to0$时上述复化公式均收敛到所求积分值$I$。

### Gauss型求积公式

把具有$2n+1$次代数精度的求积公式

$$
\int_a^b\rho(x)f(x)dx\approx\sum_{k=0}^nA_kf(x_k)
$$

称为**Gauss型求积公式，** 节点$x_k\ (k=0,1,\dots,n)$称为**高斯点，** 系数$A_k$称为**Gauss求积系数。**

构造Gauss型求积公式的关键在于确定**高斯点，** 再由$n+1$个高斯点构造**基函数**进行插值求积，从而得到高斯系数$A_k$.

在求积节点个数一定的情形下，高斯型求积公式一定是代数精度最高的求积公式$(2n+1)$次，且一定稳定。

**例：**

定义内积$(f,g)=\int_{-1}^1|x|f(x)g(x)dx$.
(1) 试建立首项系数为1的正交多项式系$\{g_i\}_{i=0}^{+\infty}$中的前三项$\{g_i\}_{i=0}^2$
(2) 建立两点Gauss型求积公式

$$
\int_{-1}^{1}|x|f(x)dx\approx A_0f(x_0)+A_1f(x_1)
$$

(1) 建立首项为1的正交多项式的前三项：(令$g_0=1$)

$$\begin{aligned}
&(g_0,g_1)=\int_{-1}^1|x|g_0g_1dx=\int_{-1}^1|x|(x+a)dx=0
\\
&(g_0,g_2)=\int_{-1}^1|x|g_0g_2dx=\int_{-1}^1|x|(x^2+bx+c)dx=0
\\
&(g_1,g_2)=\int_{-1}^1|x|g_1g_2dx=\int_{-1}^1|x|(x+a)(x^2+bx+c)dx=0
\end{aligned}$$

可以解得$a=0,b=0,c=-\frac{1}{2},$所以$g_1=x,g_2=x^2-\frac{1}{2}$ 

(2) 令最高次数多项式等于零，零点即为高斯点。令$g_2=0$,得

$$
x_0=-\frac{\sqrt{2}}{2},\, x=\frac{\sqrt{2}}{2}
$$

设所求Gauss型求积公式为$\int_{-1}^1|x|f(x)dx\approx A_0f(-\frac{\sqrt{2}}{2})+A_1f(\frac{\sqrt{2}}{2})$，其中$A_0,A_1$为待求Gauss系数，取$f(x)=1,x$使上式精确相等，有

$$
A_0+A_1=\int_{-1}^1|x|dx\\
\quad \\
-\frac{\sqrt{2}}{2}A_0+\frac{\sqrt{2}}{2}A_1=\int_{-1}^1|x|xdx
$$

解得$A_0=A_1=\frac{1}{2}$，故所求Gauss型求积公式为

$$
\int_{-1}^1|x|f(x)dx\approx \frac{1}{2}f(-\frac{\sqrt{2}}{2})+
\frac{1}{2}f(\frac{\sqrt{2}}{2})
$$

# 拉格朗日插值

$$
L_n(x)=\sum_{i=0}^nf(x_i)l_i(x_i)=\sum_{i=0}^nf(x_i)
\frac{\omega_{n+1}(x)}{(x-x_i)\omega_{n+1}(x_i)}
$$

上式是不超过$n$次的多项式，且满足所有的插值条件，因而就是我们需要构造的插值多项式，称之为**Lagrange插值多项式。**

当$n=1$时，有

$$
L_1(x)=\frac{x-x_0}{x_0-x_1}y_0+\frac{x-x_0}{x_1-x_0}y_1
$$

当$n=2$时，有

$$
L_2(x)=\frac{(x-x_1)(x-x_2)}{(x_0-x_1)(x_0-x_2)}y_0+
\frac{(x-x_0)(x-x_2)}{(x_1-x_0)(x_1-x_2)}y_1+
\frac{(x-x_0)(x-x_1)}{(x_2-x_0)(x_2-x_1)}y_2
$$

其插值余项为：

$$
R_n(x)=f(x)-L_n(x)
$$

**定理：** 如果$f^{(n)}(x)$在区间$[a,b]$上连续，$f^{(n+1)}(x)$在$(a,b)$内存在，$L_n(x)$为在节点$0\leq x_0<x_1<\dots<x_n\leq b$上满足插值条件的$n$次**Lagrange插值多项式，** 则对任一$x\in(a,b)$，其插值余项为：

$$
R_n(x)=f(x)-L_n(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_{n+1}(x)
$$

其中$\xi\in(a,b)$且依赖于$x$。上式给出的余项通常称为 **Lagrange型余项。**

一般情况下，余项表达式中的$\xi\in(a,b)$的具体数值无法知道，但是如果能求出$\max_{a\leq x\leq b}|f^{(n+1)}(x)|=M_{n+1}$，则可以得出插值多项式的截断误差限为：

$$
|R_n(x)|=|f(x)-L_n(x)|\leq\frac{M_{n+1}}{(n+1)!}\max_{a\leq x\leq b}
|\omega_{n+1}(x)|
$$

**例：**

设函数$f(x)$在区间$[a,b]$上具有二阶连续导数，且$f(a)=f(b)=0$. 证明如下不等式成立：

$$
\max_{a\leq x\leq b}|f(x)|\leq\frac{(b-a)^2}{8}\max_{a\leq x\leq b}|f''(x)|
$$

1. 写出插值函数：
   
   $$
   L(x)=\frac{x-a}{b-a}f(a)+\frac{x-b}{a-b}f(b)=0
   $$

2. 考虑余项：
   
   $$
   R(x)=f(x)-L(x)=\frac{f''(x)}{2!}(x-a)(x-b)
   $$

3. 求余项的最大值：
   
   $$\begin{aligned}
   \max_{x\in[a,b]}|f(x)|&=\max_{x\in[a,b]}|\frac{f''(x)}{2!}(x-a)(x-b)|\\
\quad \\
&\leq \frac{1}{2}\max_{x\in[a,b]}|f''(x)|\frac{(b-a)^2}{4}\\
\quad \\
&=\frac{(b-a)^2}{8}\max_{x\in[a,b]}|f''(x)|
   \end{aligned}$$
   
   证毕。

