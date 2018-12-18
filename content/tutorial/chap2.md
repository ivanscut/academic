+++
title = "Chapter 2"

date = 2018-09-09T00:00:00
# lastmod = 2018-09-09T00:00:00

draft = false  # Is this a draft? true/false
toc = true  # Show table of contents? true/false
type = "docs"  # Do not modify.

# Add menu entry to sidebar.
linktitle = "Chapter 2"
[menu.tutorial]
  parent = "Mathematical Statistics"
  weight = 2

# Featured image.
# Uncomment below lines to use.
# [header]
# image = "headers/getting-started.png"
# caption = "Image credit: [**Academic**](https://github.com/gcushen/hugo-academic/)"
+++

Let us consider fitting a straight line, $y = \beta_0+\beta_1x$, to points $(x_i,y_i)$, where $i=1,\dots,n$. 

1. Write down the normal equations for the simple linear model via the matrix formalism.

2. Solve the normal equations by tha matrix approach and see whether the solutions agree with the earlier calculation derived in the simple linear models.

`Solution`:

1. Let $Y=(y_1,\dots,y_n)^\top$, $\beta=(\beta_0,\dots,\beta_{p-1})^\top$,
$\epsilon=(\epsilon_1,\dots,\epsilon_n)^\top$, and let $X$ be the $n\times 2$ matrix
$$
X=
\left[
\begin{matrix}
1 & x_1\\
1 & x_2\\
\vdots & \\
1 & x_n\\
\end{matrix}
\right].
$$

The model can be rewritten as $$Y=X\beta+\epsilon.$$

The normal equations are $(X^\top X) \hat\beta = X^\top Y$. By simple algebra, we have

$$X^\top X = \left[
\begin{matrix}
1 & 1 & \dots & 1\\
x_1 & x_2 &\dots & x_n\\
\end{matrix}
\right]\left[
\begin{matrix}
1 & x_1\\
1 & x_2\\
\vdots & \\
1 & x_n\\
\end{matrix}
\right]=\left[
\begin{matrix}
n & \sum_{i=1}^n x_i\\
\sum_{i=1}^n x_i & \sum_{i=1}^n x_i^2
\end{matrix}
\right]$$

$$X^\top Y=\left[\begin{matrix}
1 & 1 & \dots & 1\\
x_1 & x_2 &\dots & x_n\\
\end{matrix}
\right]\left[\begin{matrix}
y_1\\
y_2\\
\vdots\\
y_n
\end{matrix}
\right]=\left[\begin{matrix}
\sum_{i=1}^n y_i\\
\sum_{i=1}^n x_iy_i
\end{matrix}
\right].
$$

The normal equations turn out to be

$$\left[
\begin{matrix}
n & \sum_{i=1}^n x_i\\
\sum_{i=1}^n x_i & \sum_{i=1}^n x_i^2
\end{matrix}
\right]\left[\begin{matrix}
\hat\beta_0\\
\hat\beta_1
\end{matrix}
\right]=\left[\begin{matrix}
\sum_{i=1}^n y_i\\
\sum_{i=1}^n x_iy_i
\end{matrix}
\right]$$

As a result,

$$
\begin{align}
\left[\begin{matrix}
\hat\beta_0\\
\hat\beta_1
\end{matrix}
\right]&=\left[\begin{matrix}
n & \sum_{i=1}^n x_i\\
\sum_{i=1}^n x_i & \sum_{i=1}^n x_i^2
\end{matrix}
\right]^{-1}\left[\begin{matrix}
\sum_{i=1}^n y_i\\
\sum_{i=1}^n x_iy_i
\end{matrix}
\right]\\
&=\frac{1}{n\sum_{i=1}^n x_i^2-(\sum_{i=1}^n x_i)^2}\left[\begin{matrix}
\sum_{i=1}^n x_i^2 & -\sum_{i=1}^n x_i\\
-\sum_{i=1}^n x_i & n
\end{matrix}
\right]\left[\begin{matrix}
\sum_{i=1}^n y_i\\
\sum_{i=1}^n x_iy_i
\end{matrix}
\right]\\
&=\frac{1}{n\sum_{i=1}^n x_i^2-(\sum_{i=1}^n x_i)^2}\left[\begin{matrix}
\sum_{i=1}^n x_i^2\sum_{i=1}^n y_i-\sum_{i=1}^n x_i\sum_{i=1}^n x_iy_i\\
n\sum_{i=1}^n x_iy_i-\sum_{i=1}^n x_i\sum_{i=1}^n y_i
\end{matrix}
\right]\\
&=\frac{1}{\ell_{xx}}\left[\begin{matrix}
\bar y\sum_{i=1}^n x_i^2-\bar x\sum_{i=1}^n x_iy_i\\
\ell_{xy}
\end{matrix}
\right]\\
&=\frac{1}{\ell_{xx}}\left[\begin{matrix}
\bar y\ell_{xx}-\bar x\ell_{xy}\\
\ell_{xy}
\end{matrix}
\right]=\left[\begin{matrix}
\bar y-\bar x\ell_{xy}/\ell_{xx}\\
\ell_{xy}/\ell_{xx}
\end{matrix}
\right],
\end{align}
$$

where $\ell_{xx} = \sum_{i=1}^n(x_i-\bar x)^2,$ $\ell_{yy} = \sum_{i=1}^n(y_i-\bar y)^2,$ $\ell_{xy} = \sum_{i=1}^n(x_i-\bar x)(y_i-\bar y).$

The solutions agree with the earlier calculation derived in the simple linear models.


---

Prove that the projection matrix $P=X(X^\top X)^{-1} X^\top$ has an eigenvalue 1, and 
$(1,\dots,1)^\top$ is one of the associated eigenvectors.

`Proof`: Note that $PX = X(X^\top X)^{-1} X^\top X  = X$. The first column of $X$ is $\mathbf{1}:=(1,\dots,1)^\top$. This implies that $P \mathbf{1} = \mathbf{1}$, which completes the proof.


---

(The QR Method) This problem outlines the basic ideas of an alternative method,
the QR method, of finding the least squares estimate $\hat \beta$. An advantage of the
method is that it does not include forming the matrix $X^\top X$, a process that tends
to increase rounding error. The essential ingredient of the method is that if $X_{n\times p}$
has $p$ linearly independent columns, it may be factored in the form

$$
\begin{align}
X\quad &=\quad Q\quad \quad R\\
n\times p &\quad  n\times p\quad p\times p
\end{align}
$$

where the columns of $Q$ are orthogonal ($Q^\top Q = I_p$) and $R$ is upper-triangular
($r_{ij} = 0$, for $i > j$) and nonsingular. For a discussion of this decomposition and
its relationship to the Gram-Schmidt process, see [https://en.wikipedia.org/wiki/QR_decomposition](https://en.wikipedia.org/wiki/QR_decomposition).

Show that $\hat \beta = (X^\top X)^{-1}X^\top Y$ may also be expressed as $\hat \beta = R^{-1}Q^\top Y$,
or $R\hat \beta  = Q^\top Y$. Indicate how this last equation may be solved for $\hat \beta$ by back-substitution, using that $R$ is upper-triangular, and show that it is thus unnecessary
to invert $R$.


`Solution`: Since $X=QR$, $$(X^\top X)^{-1}X^\top = (R^\top Q^\top  QR)^{-1} R^\top Q^\top =(R^\top R)^{-1} R^\top Q^\top =R^{-1}  Q^\top.$$

Therefore, $\hat \beta = R^{-1}  Q^\top Y$, or equivalently, $R\hat\beta = Q^\top Y=:b=(b_1,\dots,b_p)^\top$. Since $R$ is upper-triangular, then the normal equations are

$$
\begin{align}
r_{pp} \hat\beta_{p-1} &= b_p\\
r_{p-1,p-1} \hat\beta_{p-2}+ r_{p-1,p}\hat\beta_{p-1} &= b_{p-1}\\
\vdots&\\
r_{11} \hat\beta_{0}+ r_{12}\hat\beta_{1} +\dots +r_{1p}\hat\beta_{p-1}&= b_{1}
\end{align}.
$$

This can be sloved by back-substitution:

$$
\begin{align}
\hat\beta_{p-1} &= \frac{b_p}{r_{pp}}\\
 \hat\beta_{i} &=\frac{b_{i+1}}{r_{i+1,i+1}} - \frac{1}{r_{i+1,i+1}}\sum_{j=i+2}^p r_{i+1,j}\hat\beta_{j-1},\ i=p-2,\dots,0
\end{align}.
$$

---


Consider fitting the curve $y = \beta_0x+\beta_1x^2$ to points ($x_i,y_i$), where $i = 1,\dots,n$.

1. Use the matrix formalism to find expressions for the least squares estimates
of $\beta_0$ and $\beta_1$.

2. Find an expression for the covariance matrix of the estimates.

`Solution`: Let $Y=(y_1,\dots,y_n)^\top$, $\beta=(\beta_0,\dots,\beta_{p-1})^\top$, $\epsilon=(\epsilon_1,\dots,\epsilon_n)^\top$, and let $X$ be the $n\times 2$ matrix
$$
X=
\left[
\begin{matrix}
x_1 & x_1^2\\
x_2 & x_2^2\\
\vdots & \\
x_n & x_n^2\\
\end{matrix}
\right].
$$

The model can be rewritten as $$Y=X\beta+\epsilon.$$

The normal equations are $(X^\top X) \hat\beta = X^\top Y$. By simple algebra, we have

$$X^\top X = \left[
\begin{matrix}
x_1 & x_2 & \dots & x_n\\
x_1^2 & x_2^2 &\dots & x_n^2\\
\end{matrix}
\right]\left[
\begin{matrix}
x_1 & x_1^2\\
x_2 & x_2^2\\
\vdots & \\
x_n & x_n^2\\
\end{matrix}
\right]=\left[
\begin{matrix}
\sum_{i=1}^n x_i^2 & \sum_{i=1}^n x_i^3\\
\sum_{i=1}^n x_i^3 & \sum_{i=1}^n x_i^4
\end{matrix}
\right]$$

$$X^\top Y=\left[\begin{matrix}
x_1 & x_2 & \dots & x_n\\
x_1^2 & x_2^2 &\dots & x_n^2\\
\end{matrix}
\right]\left[\begin{matrix}
y_1\\
y_2\\
\vdots\\
y_n
\end{matrix}
\right]=\left[\begin{matrix}
\sum_{i=1}^n x_iy_i\\
\sum_{i=1}^n x_i^2y_i
\end{matrix}
\right].
$$

The normal equations turn out to be

$$\left[
\begin{matrix}
\sum_{i=1}^n x_i^2 & \sum_{i=1}^n x_i^3\\
\sum_{i=1}^n x_i^3 & \sum_{i=1}^n x_i^4
\end{matrix}
\right]\left[\begin{matrix}
\hat\beta_0\\
\hat\beta_1
\end{matrix}
\right]=\left[\begin{matrix}
\sum_{i=1}^n x_iy_i\\
\sum_{i=1}^n x_i^2y_i
\end{matrix}
\right]$$

Let $s_x^k = \sum_{i=1}^n x_i^k$, $s_y^k = \sum_{i=1}^n y_i^k$, $s_{xy}^{jk}=\sum_{i=1}^n x_i^jy_i^k$.
As a result,

$$
\begin{align}
\left[\begin{matrix}
\hat\beta_0\\
\hat\beta_1
\end{matrix}
\right]&=\left[
\begin{matrix}
s_x^2 & s_x^3\\
s_x^3 & s_x^4
\end{matrix}
\right]^{-1}\left[\begin{matrix}
s_{xy}^{11}\\
s_{xy}^{21}
\end{matrix}
\right]\\
&=\frac{1}{s_x^2s_x^4-(s_x^3)^2}\left[\begin{matrix}
s_x^4 & -s_x^3\\
-s_x^3 & s_x^2
\end{matrix}
\right]\left[\begin{matrix}
s_{xy}^{11}\\
s_{xy}^{21}
\end{matrix}
\right]\\
&=\left[\begin{matrix}
\frac{s_x^4s_{xy}^{11}-s_x^3s_{xy}^{21}}{s_x^2s_x^4-(s_x^3)^2}\\
\frac{s_x^2s_{xy}^{21}-s_x^3s_{xy}^{11}}{s_x^2s_x^4-(s_x^3)^2}
\end{matrix}
\right].
\end{align}
$$

Note that
$$Var[\hat\beta] = Var[(X^\top X)^{-1} X^\top Y]=(X^\top X)^{-1} X^\top Var[\epsilon] X(X^\top X)^{-1}.$$

For the usual assumption $Var[\epsilon] =\sigma^2 I_n$, then $$Var[\hat\beta]=\sigma^2 (X^\top X)^{-1} =\frac{\sigma^2}{s_x^2s_x^4-(s_x^3)^2}\left[\begin{matrix}
s_x^4 & -s_x^3\\
-s_x^3 & s_x^2
\end{matrix}
\right].$$

