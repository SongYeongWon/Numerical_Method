`#include "myNM.h"`

# Non-linear Solver

### BisectonNL()

The bisection method is a bracketing method that uses an interval which includes the true solution **Xt** in [a,b], where a,b are the upper & lower bound for the solution. Here, we have to assume or check that the function f(x) is continuous and it has a solution within the interval [a, b]

![image](https://user-images.githubusercontent.com/84532205/122443355-76ed1080-cfda-11eb-8dca-8c9347ad9b4f.png)


```c
double bisectionNL(double _a0, double _b0, double _tol)
```

#### Parameters

* continuous and it has a solution within the interval [_a0, _b0]
* tolerance _tol 

### Newton-Raphson()

This method also solves for the root of an equation of the form f(x)=0. In contrast to bisection method, this uses the slope of the equation to converge to the solution at a much faster rate.

![image](https://user-images.githubusercontent.com/84532205/122443367-7bb1c480-cfda-11eb-9d4b-c4babaefb41c.png)


This method requires the f(x) is 

* continuous 
* differentiable 
* start at the initial point near the true solution

```c
double newtonRaphson(double _x0, double _tol)
```

#### Parameters

* initial point near the true solution _x0
* tolerance _tol 

### hybridmethod()

This method combines bisection and Newton-Raphson for a fail-safe routine.

​	If Newton-Raphson gives the solution out of bounds, use bisection method for the next estimation of x(n+1)

​	If Newton-Raphson is not decreasing fast enough, or seems to be diverging, use bisection method for x(n+1)

​	Otherwise use Newton-Raphson

```c
double hybridmethod(float _x0, float _tol)
```

#### Parameters

* initial point near the true solution _x0
* tolerance _tol 

### example code

```c
    float tol = 0.00001;
    float a0 = 0;
    float b0 = 3;
    float x0 = 1;
    double BM_result, NM_result, HBM_result;

    BM_result = bisectionNL(a0, b0, tol);
    printf("Final Solution: %f \t", BM_result);

    NM_result = newtonRaphson(x0, tol);
    printf("Final Solution: %f \t", NM_result);

    HBM_result = hybridmethod(x0, tol);
    printf("Final Solution: %f \t", HBM_result);
```

See full example code:  [Assignment1.cpp]([NumericalMethods/Assignment1.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment1/Assignment1.cpp))



# linear Solver

### gaussEilm()

Solve the system of linear equation  **Ax=b**

The purpose of Gaussian elimination is to transform the linear system into an equivalent form with an upper or lower triangular Matrix. 

If matrix A is mxn(m=n) square matrix and dimension of Matrix A and B are same, we can find solution **x**.

```c
void gaussEilm(Matrix _A, Matrix _b, Matrix _U, Matrix _vecUd)
```

#### Parameters

* **Matrix A** and **B** can be expressed equation form of Ax=b
* Upper triangle **matrix U**
* **vecUd** when Ux=d

### gaussEilmJordan()

Solve the following linear systems of Ax=b

The purpose of Gaussian Jordan elimination is to make the matrix as the reduced row-echelon form of diagonal matrix with ones in the pivot elements

```c
void gaussEilmJordan(Matrix _A, Matrix _b, Matrix _U, Matrix _vecUd)
```

#### Parameters

* **Matrix A** and **B** can be expressed equation form of Ax=b
* Upper triangle **matrix U**
* **vecUd** when Ux=d

### example code

```c
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = createMat(matA.cols,matA.rows);
Matrix vecUb = createMat(vecb.rows,vecb.cols); 

gaussEilm(matA, vecb, matU, vecUb);
gaussEilmJordan(matA, vecb, matU, vecUb);
```

* #### backward substitution

![image](https://user-images.githubusercontent.com/84532205/122443089-34c3cf00-cfda-11eb-8ba2-fc2d01fcad33.png)

```c
Matrix	backSub(Matrix _A, Matrix _b)
{
	Matrix Out = createMat(_b.rows, 1);
	for (int i = _A.rows - 1; i >= 0; i--) {
		double temp = 0;
		for (int j = i + 1; j < _A.cols; j++)
			temp += _A.at[i][j] * Out.at[j][0];
		Out.at[i][0] = (_b.at[i][0] - temp) / _A.at[i][i];
	}
	return Out;
}
```

See full example code:  [Assignment2.cpp]([NumericalMethods/Assignment2.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment2/Assignment2.cpp))



### LUdecomp()

Solve the following linear systems of Ax=b, but LUdecomposition is more efficient method for finding the solution for multiple **b** vectors. Using the LU decomposition, we can easily solve for Ax=b, without finding the inverse of A.

![image](https://user-images.githubusercontent.com/84532205/122443472-97b56600-cfda-11eb-925f-733744681bd9.png)


```c
void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P)
```

#### Parameters

* Low triangle matrix **L**
* Upper triangle matrix **U**
* Permutation matrix **P**

### SolveLU()

Solve the following linear systems of Ax=b

First step :  **Ly=d** we can find **y** by **forward substitution** 

second step :  **Ux=y ** by using **backward substitution** find solution **x**

```c
Matrix solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b)
```

#### Parameters

* multiple **b** vector
* Low triangle matrix **L**
* Upper triangle matrix **U**
* Permutation matrix **P**

### Inv()

If matrix A is invertible , using inverse of matrix A, we can solve the following linear systems of Ax=b

First step : finding matrix P,L,U through LUdecomposition

Second step : using solveLU function, find element of each column of inverse matrix A

```c
Matrix inv(Matrix _A, Matrix _invA)
```

#### Parameters

* Square matrix **A**
* inverse matrix of **A**

### example code

```c
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = createMat(matA.rows,matA.cols);
Matrix vecUb = createMat(vecb.rows,vecb.cols);  
Matrix matL = createMat(matA.rows, matA.cols);
Matrix matP = createMat(matA.rows, matA.cols);
Matrix X = createMat(vecb.rows, vecb.cols);		
Matrix invA = createMat(matA.rows, matA.cols);	

gaussEilm(matA, vecb, matU, vecUb);		
LUdecomp(matA, matL, matU, matP);	
X=solveLU(matL, matU, matP, vecb);	
printMat(X, "Final solution X = \n");			
inv(matA, invA);				
```

* #### forward substitution

```c
void fwdsub(Matrix _L, Matrix _d, Matrix _y) { 
	double sub;

	for (int i = 0; i < _L.rows; i++) {
		sub = 0;
		for (int j = i-1 ; j >= 0; j--) {
			sub += _y.at[j][0] * _L.at[i][j];
		}
		_y.at[i][0] = (_d.at[i][0] - sub)/ _L.at[i][i];
	}
}
```

See full example code:  [Assignment3.cpp]([NumericalMethods/Assigment3.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment3/Assigment3.cpp))



### QRdecomp()

Solve the eigenvalues and eigenvectors of a system from the basic form of **Av=λv**

basic concept of QR decomposition method of finding eigenvalues is transforming **A** into a similar triangular Matrix **U**, through QR factorization and Iteration.

![image](https://user-images.githubusercontent.com/84532205/122443752-d77c4d80-cfda-11eb-84e6-9fe809785dcc.png)


```c
Matrix QRdecomp(Matrix _A, Matrix _Q, Matrix _R)
```

#### Parameters

* Square matrix **A** must have independent columns
* Orthogonal matrix **Q** 
* Upper triangular **R**

### eigenvalue()

Solve the following linear systems of Ax=b

First step : QRdecomposition using household matrix or Gram-Schmit

Second step : make into similar matrix 

It is repeated N times or until Ai becomes a closely approximation of a triangular matrix 

Ai is an upper triangular matrix, similar to A, having the same eigenvalues which are the diagonal elements.

```c
void eigenvalue(Matrix _A, Matrix _Q, Matrix _R)
```

#### Parameters

* Square matrix **A** must have independent columns
* Orthogonal matrix **Q**   
* Upper triangular **R**

### Cond()

Find the condition number of matrix **A**

![image](https://user-images.githubusercontent.com/84532205/122443249-5c1a9c00-cfda-11eb-8d71-6c94daf43326.png)

```c
double Cond(Matrix _A)
```

#### Parameters

* Matrix A (mxn) can be a rectangular matrix 

### example code

```c
QRdecomp(matA, matQ, matR);
eigenvalue(matA, matQ, matR);

double CondNumberA, CondNumberC = 0;
CondNumberA = Cond(matA);
CondNumberC = Cond(matC);
printf("CondNumberA = %f\n", CondNumberA);
printf("CondNumberC = %f\n", CondNumberC);
```

### Jacobian()

Newton's method for a system of non-linear equations.

Newton's method is extended from scalar variable to a vector variable as 


![image](https://user-images.githubusercontent.com/84532205/122438906-270c4a80-cfd6-11eb-9860-a73397ea5373.png)


where the matrix **J** is t he Jacobian matrix for the vector-valued functions with multiple variables, <img src="https://user-images.githubusercontent.com/84532205/122416944-1a7ef680-cfc4-11eb-8ca2-bd94c88b4ded.png" alt="image" style="zoom:67%;" />

Jacobian matrix is the matrix of first order partial derivatives

```c
Matrix Jacobian(double _x0,double _y0)
```

#### Parameters

* Initial value **x** and **y**

### example code

```c
Matrix matF = createMat(2, 1);
matF = Fmatrix(2, 3);

Matrix matdF = createMat(2, 2);
matdF = Fdmatrix(2, 3);

double a = 4 ,b = 1;
Matrix matJaco = createMat(2, 1);
matJaco = Jacobian(a, b);
printMat(matJaco, "Jacobian solution  \n");
```

See full example code:  [Assignment4.cpp]([NumericalMethods/Assignment4.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment4/Assignment4.cpp))

# CurveFitting and Interpolation

### Linearfit()

Curve Fitting with a line (first degree polynomial) is finding the best parameters (a0,a1) for the linear equation

<img src="https://user-images.githubusercontent.com/84532205/122417066-31254d80-cfc4-11eb-8e9e-aeb9d8ce37d4.png" alt="image" style="zoom:80%;" />

```c
void linearFit(Matrix _x, Matrix _y) 
```

#### Parameters

* Input dataset xi, yi vector or array

### LinearInterp()

Interpolation is used for estimating a value between the known data points.

Interpolation with a single polynomial can be used for low-order polynomials to estimate a value between points.
                                                             <img src="https://user-images.githubusercontent.com/84532205/122420128-659a0900-cfc6-11eb-944f-7b0efc898356.png" alt="image" style="zoom:80%;" />

#### Lagrange Polynomial

<img src="https://user-images.githubusercontent.com/84532205/122420168-6d59ad80-cfc6-11eb-83b3-fb69c4321fc0.png" alt="image" style="zoom:80%;" />

```c
void linearInterp(Matrix _x, Matrix _y, Matrix _xq)
```

#### Parameters

* Input dataset xi, yi vector or array
* query input **xq** , vector or array

### example code

```c
int M = 21;
double Xq_array[] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,100};
Matrix vecXq = arr2Mat(Xq_array, M, 1);
Matrix vecT = txt2Mat(path, "vecT");
Matrix vecP = txt2Mat(path, "vecP");
linearFit(vecT, vecP);
linearInterp(vecT, vecP, vecXq);
```

See full example code:  [Assignment5.cpp]([NumericalMethods/Assignment5.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment5/Assignment5.cpp))

# Numerical Differentiation

### gradient()

Find the slope of given point x and y at least two point and calculate a numerical approximation of the derivative from the given discrete datasets. 

Return the dy/dx results for the input data. We can use Two-point or Three-point and other methods depends on the number of datasets.

<img src="https://user-images.githubusercontent.com/84532205/122420816-ea852280-cfc6-11eb-845c-c8f579a5668a.png" alt="image" style="zoom:80%;" />

```c
Matrix	gradient(Matrix _x, Matrix _y)
```

#### Parameters

* Input dataset **x** and **y**

### gradient1D()

Find the slope of given point x and y at least two point using 1D-array vector 

Return the dy/dx results for the input data. 

```c
void	gradient1D(double x[], double y[], double dydx[], int m)
```

#### Parameters

* Input dataset **x** and **y**
* the number of dataset **m**

### gradientFunc()

Return the **dy/dx** results for the target equation.

```c
Matrix	gradientFunc(double func(const double x), Matrix xin) 
```

#### Parameters	

* target equation **func**
* input data **xin**

### newtonRaphsonFunc()

Modified newtonRaphson function to pass functions as input

```c
double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), double _x0, double _tol)
```

#### Parameters

* input functions **func** and **dfunc (first derivative function)** 
* initial point near the true solution **x**
* tolerance **tol **

### example code

```c
double myFunc(const double x) {
	return  x * x * x;
}
double mydFunc(const double x) {
	return 3 * x * x;
}
	double t_array[21] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
	double pos_array[21] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 	5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59};
	double vel_array[21] = {0};
	double acc_array[21] = {0};
	int m = 21;
//gradient 
	Matrix t = arr2Mat(t_array, m, 1);
	Matrix x = arr2Mat(pos_array, m, 1);
	Matrix vel = gradient(t, x);
	Matrix acc = gradient(t, vel);
	
	printMat(t, "t");
	printMat(x, "x");
	printMat(vel, "vel");
	printMat(acc, "acc");

//gradient 1D
	gradient1D(t_array, pos_array, vel_array, m);
	gradient1D(t_array, vel_array, acc_array, m);
	
	printf("t=\n");
	printarr(t_array,m);
	printf("pos=\n");
	printarr(pos_array,m);
	printf("vel=\n");
	printarr(vel_array,m);
	printf("acc=\n");
	printarr(acc_array,m);

//gradientFunc
	Matrix xin = arr2Mat(t_array, m, 1);
	Matrix dydx = gradientFunc(myFunc, xin);
	printMat(xin, "xin");
	printMat(dydx, "dydx");

//newtonRaphsonFunc
	double tol = pow(10,-20);
	double x0 = 2;
	double NM_result;
	NM_result = newtonRaphsonFunc(myFunc,mydFunc, x0, tol);
	printf("final solution: %f \t", NM_result);
```

See full example code:  [Assignment6.cpp]([NumericalMethods/Assignment6.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment6/Assignment6.cpp))

# Numerical Integration

Calculate a numerical approximation of the integration

<img src="https://user-images.githubusercontent.com/84532205/122420915-fa9d0200-cfc6-11eb-8e3b-c200ca14d422.png" alt="image" style="zoom: 80%;" />

### trapz()

Integration using trapezoidal method for discrete data inputs

for each subinterval is interpolation with a line (1st order polynomial)

<img src="https://user-images.githubusercontent.com/84532205/122421030-0e486880-cfc7-11eb-88fc-f2a6ba0d5140.png" alt="image" style="zoom:80%;" />

```c
double trapz(double x[], double y[], int m)
```

#### Parameters

* input dataset **x[]** and **y[]**
* the number of dataset **m**

### integralMid()

Integration using Mid-point method for discrete data inputs


<img src="https://user-images.githubusercontent.com/84532205/122421055-143e4980-cfc7-11eb-9c87-2d537781b330.png" alt="image" style="zoom:80%;" />



```c
double integralMid(double x[], double y[], int m) 
```

#### Parameters

* input dataset **x[]** and **y[]**
* the number of dataset **m**

### integral() //simpson13

Integration using Simpson1/3 method for discrete data inputs

we can apply the interpolation with quadratic(2nd order polynomial)

<img src="https://user-images.githubusercontent.com/84532205/122421152-27511980-cfc7-11eb-9b67-41bc5e628cc3.png" alt="image" style="zoom:80%;" />

```c
double integral(double func(const double x), double a, double b, int n)
```

#### Parameters

* Given function **func**
* Range from **a** to **b**, **n** intervals

### integral38() //simpson38

Integration using Simpson3/8 method for discrete data inputs

we can apply the interpolation with quadratic(3rd order polynomial)

<img src="https://user-images.githubusercontent.com/84532205/122421356-46e84200-cfc7-11eb-9158-be51a2488380.png" alt="image" style="zoom:80%;" />

```c
double integral38(double func(const double x), double a, double b, int n)
```

#### Parameters

* Given function **func**
* Range from **a** to **b**, **n** intervals

### example code

```c
double myFunc(const double x) {
	return sqrt(1 - x * x);
}
	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
	int M = sizeof(x) / sizeof(x[0]);
	double I_rect = IntegrateRect(x, y, M);
	double I_trapz = trapz(x, y, M);
	double l_midpoint = integralMid(x, y, M);

	printf("I_rect  = %f\n", I_rect);
	printf("I_trapz = %f\n", I_trapz);
	printf("interalMid = %f\n", l_midpoint);

	double a = -1;
	double b = 1;
	int n = 12;
	double I_simpson13 = integral(myFunc, a, b, n);
	double l_simson38 = integral38(myFunc, a, b, n);
	printf("Simpson13 = %f\n", I_simpson13);
	printf("Simpson38 = %f\n", l_simson38);
```

See full example code:  [Assignment7.cpp]([NumericalMethods/Assignment7.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment7/Assignment7.cpp))

# ODE Numerical Method

### Basic mechanism of ODE

<img src="https://user-images.githubusercontent.com/84532205/122421729-8c0c7400-cfc7-11eb-8f1b-41ba9e0a07f7.png" alt="image" style="zoom: 80%;" />

There are some popular ODE method Euler's method, Midpoint method, Runge-Kutta method depends on how we use numerical integration.

### odeEU()

Euler's Explicit Method 

We can find next point yi+1,by using initial condition and the function of slope  at f(x,y)=dy/dx , y(x0)=y0

<img src="https://user-images.githubusercontent.com/84532205/122421978-bceca900-cfc7-11eb-85b3-8395869d9117.png" alt="image" style="zoom:80%;" />

```c
void odeEU(double func(const double x, const double y), double x[], double y[], double t0, double tf, double h)
```

#### Parameters

* **dy/dx** function **func**
* the independent variable **x** is in the range of **[t0, tf]**
* stepsize **h**

### odeEM()

Modified Euler's Method 

We can find next point yi+1 by using initial condition and the average of the two slopes. Slope at the beginning of the slope and estimate of the slope at the end of the interval. 

<img src="https://user-images.githubusercontent.com/84532205/122422104-d988e100-cfc7-11eb-97e9-8925f2203a4b.png" alt="image" style="zoom:80%;" />

```c
void odeEM(double func(const double x, const double y), double x[], double y[], double t0, double tf, double h)
```

#### Parameters

* **dy/dx** function **func**
* the independent variable **x** is in the range of **[t0, tf]**
* stepsize **h**

### ode()

Calls different ODE method

```c
void ode(double func(const double x, const double y), double x[], double y[], double t0, double tf, double h, int method)
```

#### Parameters

* **dy/dx** function **func**
* the independent variable **x** is in the range of **[t0, tf]**
* stepsize **h**
* defined method number **method**

### example code

```c
double myFunc(const double x, const double y) {
	double tau = 1;
	double Vm = 1;
	double w = 2 * Pi * 10;
	double f = 10;
	return (-1 / tau) * y + (1 / tau) * Vm * cos(2 * Pi * f * x);
}
//odeEU
	double h = 0.001;
	double t0 = 0;
	double tf = 0.1;
	double x_array[101] = { 0 };
	double y_array[101] = { 0 };
	odeEU(myFunc, x_array, y_array, t0, tf, h);
	printf("x \n");
	printarr(x_array, 101);
	printf("y \n");
	printarr(y_array, 101);

//odeEM
	double x1_array[101] = { 0 };
	double y1_array[101] = { 0 };
	odeEM(myFunc, x1_array, y1_array, t0, tf, h);
	printf("y \n");
	printarr(y1_array, 101);

//ode
	ode(myFunc, x_array, y_array, t0, tf, h, Em);
```

See full example code:  [Assignment8.cpp]([NumericalMethods/Assignment8.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment8/Assignment8.cpp))



### odeRK2()

Second-order Runge-Kutta Method  

Find next value with present value and slope 

<img src="https://user-images.githubusercontent.com/84532205/122422164-e4dc0c80-cfc7-11eb-8a61-69093ea03e8c.png" alt="image" style="zoom:80%;" />

General form of second order Runge-Kutta

<img src="https://user-images.githubusercontent.com/84532205/122422178-e9a0c080-cfc7-11eb-84a1-86d08043bfb4.png" alt="image" style="zoom:80%;" />

```c
void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
```

#### Parameters

* **dy/dx** function **odefunc**
* the independent variable **x** is in the range of **[t0, tf]**
* stepsize **h**
* initial condition **y0**

### odeRK4()

Fourth-order Runge-Kutta Method  

Find next value with present value and slope 

<img src="https://user-images.githubusercontent.com/84532205/122423516-ece87c00-cfc8-11eb-90d7-54d8f39c0f85.png" alt="image" style="zoom:80%;" />

Classical Fourth-order Runge-Kutta

<img src="https://user-images.githubusercontent.com/84532205/122423461-e5c16e00-cfc8-11eb-9394-d86e9918ad2d.png" alt="image" style="zoom:80%;" />

```c
void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
```

#### Parameters

* **dy/dx** function **odefunc**
* the independent variable **x** is in the range of **[t0, tf]**
* stepsize **h**
* initial condition **y0**



# Higher-order IVP

### sys2RK2()

we can solve one of 2nd order ODE in the form of two of 1st order ODE using RK2 method 

```c
void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
```

#### Parameters

* two dependent variables **y1[]** and **y2[]**
* the independent variable **x** is in the range of **[t0, tf]**
* stepsize **h**
* initial condition **y1_init,y2_init**

### sys2RK4()

we can solve one of 2nd order ODE in the form of two of 1st order ODE using RK4 method 

```c
void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
```

#### Parameters

* two dependent variables **y1[]** and **y2[]**
* the independent variable **x** is in the range of **[t0, tf]**
* stepsize **h**
* initial condition **y1_init,y2_init**

### example code

```c
double odeFunc_rc(const double t, const double v) {
	double tau = 1;
	double Vm = 1;
	double w = 2 * PI * 10;
	double f = 10;
	return (-1 / tau) * v + (1 / tau) * Vm * cos(2 * PI * f * t);
}
void odeFunc_mck(const double t, const double Y[], double dYdt[])
{
	double m = 1;
	double c = 7;
	double k = 6.9;
	double f = 5;
	double Fin = 2 * cos(2 * PI * f * t);  //y2dot

	dYdt[0] = Y[1];
	dYdt[1] = (-k * Y[0] - c * Y[1] + Fin);
}

	//Parameter Definitions
	double a = 0;
	double b = 0.1;
	double h = 0.001;
	unsigned int N = (b - a) / h + 1;
	double y_EU[200] = { 0 };				
	double y_EM[200] = { 0 };
	double y_RK2[200] = { 0 };
	double y_RK4[200] = { 0 };

	// Initial value
	double v0 = 0;

	// ODE solver
	odeEU(odeFunc_rc, y_EU, a, b, h, v0);
	odeEM(odeFunc_rc, y_EM, a, b, h, v0);
	odeRK2(odeFunc_rc, y_RK2, a, b, h, v0);
	odeRK4(odeFunc_rc, y_RK4, a, b, h, v0);

	// Print outputs
	printf("/*-----------------------*/\n");
	printf("/ Single of 1st Order ODE /\n");
	printf("/*-----------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++)
		printf("t= %f\tyEU= %f\tyEM= %f\tyRK2= %f\tyRK4= %f\n", a + i * h, y_EU[i], y_EM[i], y_RK2[i], y_RK4[i]);


	//Parameter Definitions
	double t0 = 0;
	double tf = 1;
	h = 0.01;
	N = (tf - t0) / h + 1;
	double y[200] = { 0 };
	double v[200] = { 0 };
	double y1[200] = { 0 };
	double v1[200] = { 0 };
	// Initial values
	double y0 = 0;
	v0 = 0.2;
	// ODE solver: RK2 RK4
	sys2RK2(odeFunc_mck, y, v, t0, tf, h, y0, v0);
	sys2RK4(odeFunc_mck, y1, v1, t0, tf, h, y0, v0);

	// Print outputs
	printf("/*---------------------------*/\n");
	printf("/ 2nd Order ODE : MCK example /\n");
	printf("/*---------------------------*/\n");
	printf(" - Total number of data N=%d \n", N);

	for (int i = 0; i < N; i++)
		printf("|RK2| t= %f\ty= %f\tv= %f\t |RK4| t= %f\ty= %f\tv= %f\n", t0 + i * h, y[i], v[i], t0 + i * h, y1[i], v1[i]);
	printf("\n\n");
```

See full example code [Assignment9]([NumericalMethods/Assignment9.cpp at main · wonhg1446/NumericalMethods (github.com)](https://github.com/wonhg1446/NumericalMethods/blob/main/sources/Assignment9/Assignment9.cpp))
