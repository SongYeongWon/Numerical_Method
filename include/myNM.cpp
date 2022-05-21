/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Song-Yeong-Won
Created          : 26-03-2018
Modified         : 29-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"
//matrix multiply
Matrix multMat(Matrix _A, Matrix _B) {
	Matrix _Out = createMat(_A.rows, _B.cols);
	int i, j, k;
	initMat(_Out, 0);
	if (_A.cols != _B.rows)
		printf("Can't multiply these two Matricis\n");
	else {
		for (i = 0; i < _A.rows; i++) {
			for (j = 0; j < _B.cols; j++) {
				for (k = 0; k < _A.cols; k++)
					_Out.at[i][j] += _A.at[i][k] * _B.at[k][j];
			}
		}
	}
	return _Out;
}


void   multiMat2(Matrix _A, Matrix _B, Matrix _C)
{
	initMat(_C, 0);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _B.cols; j++) {
			for (int k = 0; k < _A.cols; k++)
				_C.at[i][j] += _A.at[i][k] * _B.at[k][j];
		}
	}
}
// Matrix addition
Matrix	addMat(Matrix _A, Matrix _B)
{
	if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] + _B.at[i][j];

	return Out;
}

Matrix	subMat(Matrix _A, Matrix _B)
{
	/*if (_A.rows != _B.rows || _A.cols != _B.cols) {
		printf("\n*************************************************");
		printf("\n  ERROR!!: dimension error at 'addMat' function");
		printf("\n*************************************************\n");
		return createMat(0, 0);
	}*/

	Matrix Out = createMat(_A.rows, _B.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			Out.at[i][j] = _A.at[i][j] - _B.at[i][j];

	return Out;
}


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

void gaussEilm(Matrix _A, Matrix _b, Matrix _U, Matrix _vecUd) // definition of Guass elimination function
															   // find the solution by using guass elimination 
{
	for (int i = 0; i < _A.rows; i++) {          //mat A 값을 mat U 에 Copy
		for (int j = 0; j < _A.cols; j++) {
			_U.at[i][j] = _A.at[i][j];
		}
	}

	for (int i = 0; i < _b.rows; i++) {          //vecb 값을 vecUd 에 Copy
		for (int j = 0; j < _b.cols; j++) {
			_vecUd.at[i][j] = _b.at[i][j];
		}
	}
	double term = 0;

	if (_A.rows == _A.cols) { //consider matrix A is squre matrix
		printMat(_A, "matA");
		printMat(_b, "vector b");
		if (_A.rows == _b.rows) { //consider dimension of A,b are appropriate			
			if (_A.at[0][0] != 0) { //consider division by zero 
				for (int k = 0; k < _U.rows - 1; k++) {
					for (int i = k + 1; i < _U.rows; i++)
					{
						term = _U.at[i][k] / _U.at[k][k];
						for (int j = 0; j < _U.cols; j++)
						{
							_U.at[i][j] = _U.at[i][j] - (term * _U.at[k][j]);
							if (j < _vecUd.cols) {
								_vecUd.at[i][j] = _vecUd.at[i][j] - term * _vecUd.at[k][j];
							}
						}
					}
				}
				printMat(_U, "Upper matrix U"); //show Upper triangle matrix U
				printMat(_vecUd, "vector d ");  //show vecd when Ux=d
				Matrix matX = createMat(_vecUd.rows, _vecUd.cols); //creat solution matrix X
				matX = backSub(_U, _vecUd);
				printMat(matX, "solution matrix X"); //show Final solution matrix X
			}
			else
			{
				printf("\n");
				printf(" we can't find solution, beacase matrix element division by zero \n ");
			}
		}
		else
		{
			printf("\n");
			printf("dimension of Matrix A and b are not appropriate\n"); 
		}
	}
	else
	{
		printMat(_A, "matrix A");
		printf("\n");
		printf("Matrix A is not square matrix\n");
	}
}

void gaussEilmJordan(Matrix _A, Matrix _b, Matrix _U, Matrix _vecUd) // definition of Guass elimination function
															   // find the solution by using guass elimination 
{
	for (int i = 0; i < _A.rows; i++) {          //mat A 값을 mat U 에 Copy
		for (int j = 0; j < _A.cols; j++) {
			_U.at[i][j] = _A.at[i][j];
		}
	}

	for (int i = 0; i < _b.rows; i++) {          //vecb 값을 vecUd 에 Copy
		for (int j = 0; j < _b.cols; j++) {
			_vecUd.at[i][j] = _b.at[i][j];
		}
	}
	double term = 0;
	double max, maxrow, firstrow, pivotterm, temp, initialmax = 0;
	Matrix newcol = createMat(_U.rows, 1);							//newcol 은 각 partial pivot 의 값을 저장하는 행렬
	initMat(newcol, 0);
	if (_A.rows == _A.cols) { //consider matrix A is squre matrix
		//printMat(_A, "matA");
		//printMat(_b, "vector b");
		if (_A.rows == _b.rows) { //consider dimension of A,b are appropriate			
				for (int k = 0; k < _U.rows - 1; k++) {
					int e = 0;								// pivot 행의 번호를 체크하기 위해서 e 라는 변수 선언
					for (int c = k; c < _U.cols; c++)
					{
						firstrow = fabs(_U.at[c][k]);		//firstrow 변수는 비교기준이 되는 첫번째 행렬값 
						for (int p = 1; p < _U.cols; p++)
						{
							maxrow = fabs(_U.at[c][p]);		//firstrow 가 있는 행에서 열을 이동하면서 값을 읽음
							if (maxrow > firstrow) {		//maxrow 가  firstrow 보다 크면 firstrow 값을 maxrow 로 update
								firstrow = maxrow;
							}
						}
						max = fabs(_U.at[c][k]) / firstrow;		//partial pivot 크기
						newcol.at[c][0] = max;					//그 값을 newcol 행에 저장
					}
					initialmax = fabs(newcol.at[k][0]);		//newcol 행렬의 첫번째 행을 pivot 비교 기준으로 함

					for (int p = k + 1; p < _U.cols; p++) {			//newcol 행렬에 저장된 각 행의 partial pivot 값을 비교해서 pivoting 될 행의 번호를 결정
						pivotterm = fabs(newcol.at[p][0]);		 //다음 행의 값을 읽음.

						if (pivotterm > initialmax) {
							e = p;									//pivot 되는 행의 번호를 e 에 저장
							initialmax = fabs(newcol.at[e][0]);		//pivotterm 이 initialmax 보다 크면 initialmax 값을 pivotterm으로 바꿔주고 계속 비교 진행 
						}
					}
					if (e >= 1) {

						for (int c = 0; c < _U.cols; c++)			//기준이 되는 행과 pivot 행의 열을 하나씩 변경
						{
							// Matrix L,U,P 모두다 partial pivot 시행
							temp = _U.at[e][c];						//pivot 되어야 하는 행의 열을 term 에 저장
							_U.at[e][c] = _U.at[k][c];				//pivot 행에다가 기준이 되는 행의 값을 저장
							_U.at[k][c] = temp;						//기준이 되는 행에 pivot 행을 저장

							temp = _vecUd.at[e][c];
							_vecUd.at[e][c] = _vecUd.at[k][c];
							_vecUd.at[k][c] = temp;
						}
					}		
					// L 의 형태로 만드는 과정
					for (int i = k + 1; i < _U.rows; i++)
					{
						term = _U.at[i][k] / _U.at[k][k];
						for (int j = 0; j < _U.cols; j++)
						{
							_U.at[i][j] = _U.at[i][j] - (term * _U.at[k][j]);
							if (j < _vecUd.cols) {
								_vecUd.at[i][j] = _vecUd.at[i][j] - term * _vecUd.at[k][j];
							}
						}
					}
					
				}
				//여기서 부터는 반대로 L행태에서 올라가면서 윗행을 0으로 만들어주는 code
					double temp = 0;
					for (int k = _U.rows - 1; k >= 0; k--) {
						for (int i = k - 1; i >= 0; i--) {
							temp = _U.at[i][k] / _U.at[k][k];
							for (int j = 0; j < _U.cols; j++)
							{
								_U.at[i][j] = _U.at[i][j] - (temp * _U.at[k][j]);
								if (j < _vecUd.cols) {
									_vecUd.at[i][j] = _vecUd.at[i][j] - temp * _vecUd.at[k][j];
								}
							}
						}
					}
					//각 diagonal의 크기로 나눠주어 I 를 만들어줌 
					double diagonal = 0;
					for (int i = 0; i < _U.rows; i++) {
						
						for (int j = 0; j < _U.cols; j++)
						{
							diagonal = _U.at[i][i];
							_vecUd.at[i][j] = _vecUd.at[i][j] / diagonal;
							_U.at[i][i] = _U.at[i][i] / diagonal;
						}
						
					}
					//여기 백슬래쉬 풀어야 출력값 볼수 있음
					printMat(_U, "Upper matrix U"); //show Upper triangle matrix U
					printMat(_vecUd, "vector d ");  //show vecd when Ux=d

			}
			else
			{
				printf("\n");
				printf(" we can't find solution, beacase matrix element division by zero \n ");
			}
		}
		else
		{
			printMat(_A, "matrix A");
			printf("\n");
			printf("Matrix A is not square matrix\n");
		}
}

void LUdecomp(Matrix _A, Matrix _L, Matrix _U, Matrix _P) // definition of Guass elimination function
														  // find the solution by using guass elimination 
{
	Matrix DI = createMat(_A.rows, _A.cols); 
	for (int i = 0; i < _A.rows; i++)				//creat diagonal matrix 
	{						
		for (int j = 0; j < _A.cols; j++) {
			if (i == j)
			{
				DI.at[i][j] = 1;
			}
			else
			{
				DI.at[i][j] = 0;
			}
		}
	}

	for (int i = 0; i < _A.rows; i++)		// matrix L 0 으로 초기화
	{ 
		for (int j = 0; j < _A.cols; j++) {
				_L.at[i][j] = 0;
		}
	}

	for (int i = 0; i < _A.rows; i++) {				// permuation matrix P 초기화
		for (int j = 0; j < _A.cols; j++) {
			if (i == j)
			{
				_P.at[i][j] = 1;
			}
			else
			{
				_P.at[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < _A.rows; i++) {          //mat A 값을 mat U 에 Copy
		for (int j = 0; j < _A.cols; j++) {
			_U.at[i][j] = _A.at[i][j];
		}
	}
	double term,max,maxrow,firstrow,pivotterm,temp,initialmax= 0;  
	Matrix newcol = createMat(_U.rows, 1);							//newcol 은 각 partial pivot 의 값을 저장하는 행렬
	initMat(newcol, 0);
																	// pivot 행의 번호를 체크하기 위해서 e 라는 변수 선언
	if (_A.rows == _A.cols)											//consider matrix A is squre matrix
	{
		//printMat(_A, "matA");		
				for (int k = 0; k < _U.rows - 1; k++) {
					int e = 0;								// pivot 행의 번호를 체크하기 위해서 e 라는 변수 선언
					for (int c = k; c < _U.cols; c++)		
					{
						firstrow = fabs(_U.at[c][k]);		//firstrow 변수는 비교기준이 되는 첫번째 행렬값 
						for (int p = 1; p < _U.cols; p++)
						{			
							maxrow = fabs(_U.at[c][p]);		//firstrow 가 있는 행에서 열을 이동하면서 값을 읽음
							if (maxrow > firstrow) {		//maxrow 가  firstrow 보다 크면 firstrow 값을 maxrow 로 update
								firstrow = maxrow;			
							}
						}
						max = fabs(_U.at[c][k]) / firstrow;		//partial pivot 크기
						newcol.at[c][0] = max;					//그 값을 newcol 행에 저장
					}
						initialmax = fabs(newcol.at[k][0]);		//newcol 행렬의 첫번째 행을 pivot 비교 기준으로 함
										
					for (int p = k+1; p < _U.cols; p++) {			//newcol 행렬에 저장된 각 행의 partial pivot 값을 비교해서 pivoting 될 행의 번호를 결정
							pivotterm = fabs(newcol.at[p][0]);		 //다음 행의 값을 읽음.
							
							if (pivotterm > initialmax) {								
								e = p;									//pivot 되는 행의 번호를 e 에 저장
								initialmax = fabs(newcol.at[e][0]);		//pivotterm 이 initialmax 보다 크면 initialmax 값을 pivotterm으로 바꿔주고 계속 비교 진행 
							}
						}	
					if (e >=1) {
						
						for (int c = 0; c < _L.cols; c++)			//기준이 되는 행과 pivot 행의 열을 하나씩 변경
						{
																	// Matrix L,U,P 모두다 partial pivot 시행
							temp = _U.at[e][c];						//pivot 되어야 하는 행의 열을 term 에 저장
							_U.at[e][c] = _U.at[k][c];				//pivot 행에다가 기준이 되는 행의 값을 저장
							_U.at[k][c] = temp;						//기준이 되는 행에 pivot 행을 저장

							temp = _L.at[e][c];
							_L.at[e][c] = _L.at[k][c];
							_L.at[k][c] = temp;

							temp = _P.at[e][c];
							_P.at[e][c] = _P.at[k][c];
							_P.at[k][c] = temp;
						}
					}
					for (int i = k + 1; i < _U.rows; i++)
					{
						term = _U.at[i][k] / _U.at[k][k];
						
						_L.at[i][k] = _U.at[i][k] / _U.at[k][k];
						
						for (int j = k; j < _U.cols; j++)
						{
							_U.at[i][j] = _U.at[i][j] - (term * _U.at[k][j]);
						}
					}
					//printMat(_U, "UUUU\n");
				}
				Matrix LI = createMat(_A.rows, _A.cols);
				LI=addMat(_L, DI);  // L = L + I
				//printMat(_U, "Upper matrix U"); //show Upper triangle matrix U
				//printMat(LI, "Lower matirx L");
				//printMat(_P, "permutation matrix P");
	}else
	{
		printMat(_A, "matrix A");
		printf("\n");
		printf("Matrix A is not square matrix\n");
	}
}

Matrix solveLU(Matrix _L, Matrix _U, Matrix _P, Matrix _b) {
	Matrix DI = createMat(_L.rows, _L.cols);
	for (int i = 0; i < _L.rows; i++) {
		for (int j = 0; j < _L.cols; j++) {
			if (i == j)
			{
				DI.at[i][j] = 1;
			}
			else
			{
				DI.at[i][j] = 0;
			}
		}
	}
	_L = addMat(_L, DI);						

	Matrix vecy = createMat(_b.rows, _b.cols);  //Ly=d 에 사용될 vecy 
	initMat(vecy, 0);

	Matrix d=createMat(_b.rows, _b.cols);		//Pb = d 행렬 d 생성
	d = multMat(_P, _b);						//d = Pb , b는 처음 주어진 Ax=b
	fwdsub(_L, d, vecy);						//Forward sub 을 통해 y 값 구함 
	
	Matrix X = createMat(_b.rows, _b.cols);
	X= backSub(_U, vecy);						//Ux=y , Backsub 을 통해서 최종 x 값 도출
	return X;
}


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

void inv(Matrix _A, Matrix _invA) {						
	Matrix matX = createMat(_A.rows, _A.cols); 
	Matrix b= createMat(_A.rows, 1);
	double term,temp = 0;
	Matrix matX2 = createMat(_A.rows, _A.cols); //creat solution matrix X2
	initMat(matX2, 0);
	Matrix matL = createMat(_A.rows, _A.cols);
	Matrix matU = createMat(_A.rows, _A.cols);
	Matrix matP = createMat(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++) {
			if (i == j)
			{
				matP.at[i][j] = 1;
			}
			else
			{
				matP.at[i][j] = 0;
			}
		}
	}
	LUdecomp(_A, matL, matU, matP);				//주어진 Matrix A 를 먼저 LU decomp 함.
	double det=1;
	for (int i = 0; i < matU.rows; i++) {	   //det U 계산
			det *=  matU.at[i][i];
	}
	if (_A.rows == _A.cols)											//consider matrix A is squre matrix
	{
		if (det != 0) {									//det = 0 이 아니면 실행
			for (int k = 0; k < _A.rows; k++) {
				for (int i = 0; i < _A.rows; i++) {
					if (i == k) {
						b.at[i][0] = 1;
					}
					else
					{
						b.at[i][0] = 0;
					}
				}
				matX = solveLU(matL, matU, matP, b); //solveLU 를 통해 구한 백터 값을 matX 에 저장
				for (int i = 0; i < matX2.cols; i++) {
					temp = matX.at[i][0];			//matX 에 저장된 행의 값을 temp 에 저장
					matX2.at[i][k] = temp;			//matX2 의 각 행에 temp 값 저장
				}
			}
			for (int i = 0; i < _A.rows; i++) {          //mat A 값을 mat U 에 Copy
				for (int j = 0; j < _A.cols; j++) {
					_invA.at[i][j] = matX2.at[i][j];
				}
			}
			
			//printMat(matX2, "Inv A is =\n");	 //show Final solution matrix X

		}
		else
		{
			printf("det = %f\n", det);
			printf("determinant is zero\n");
		}
	}
	else
	{
		printMat(_A, "matrix A");
		printf("\n");
		printf("Matrix A is not square matrix\n");
	}
}

Matrix QRdecomp(Matrix _A, Matrix _Q, Matrix _R) {
	Matrix vecC = createMat(_A.rows, 1);
	Matrix vecV = createMat(_A.rows, 1);
	Matrix vecU = createMat(_A.rows, 1);
	Matrix vecUT = createMat(1, _A.cols); 
	Matrix vece = createMat(_A.rows, 1);
	Matrix matUUT = createMat(_A.rows, _A.cols);
	Matrix matH = createMat(_A.rows, _A.cols);
	Matrix DI = createMat(_A.rows, _A.rows);
	
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.rows; j++) {
			if (i == j)
			{
				DI.at[i][j] = 1;
			}
			else
			{
				DI.at[i][j] = 0;
			}
		}
	}

	for (int i = 0; i < _A.rows; i++) {          //mat DI 값을 mat Q 에 Copy
		for (int j = 0; j < _A.rows; j++) {
			_Q.at[i][j] = DI.at[i][j];
		}
	}
	for (int i = 0; i < _A.rows; i++) {          //mat A 값을 mat R 에 Copy
		for (int j = 0; j < _A.cols; j++) {
			_R.at[i][j] = _A.at[i][j];
		}
	}
	for (int k = 0; k < _A.rows-1; k++) {
		initMat(vecC, 0);
		initMat(vecV, 0);
		initMat(vece, 0);
		double normvalueC, normvalueV= 0;
		for (int i = k; i < _A.rows; i++)
		{
			vecC.at[i][0] = _R.at[i][k];
		}
		for (int e = 0; e < _R.rows; e++) {
			if (vecC.at[k][0] >= 0) {

				vece.at[k][0] = 1;
			}
			else
			{
				vece.at[k][0] = -1;
			}
		}
		normvalueC = norm(vecC);

		for (int i = 0; i < _R.rows; i++)
		{
			vecV.at[i][0] = vecC.at[i][0] + (normvalueC * vece.at[i][0]);
		}
		normvalueV = norm(vecV);
		for (int i = 0; i < _R.rows; i++)
		{
			vecU.at[i][0] = vecV.at[i][0] / normvalueV;
		}
		vecUT = transpose(vecU);
		matUUT= multMat(vecU, vecUT);
		double num = 0;
		for (int i = 0; i < matUUT.rows ; i++)
		{
			for (int j = 0; j < matUUT.cols ; j++)
			{
				num = 2 * matUUT.at[i][j];
				matUUT.at[i][j] = num;
			}
			
		}
		matH = subMat(DI, matUUT);

		_Q = multMat(_Q, matH);
		_R = multMat(matH, _R);
		//printf("k = %d\n", k);
	
	}
	//printMat(matH, "matH\n");
	//printMat(_Q, "matrix Q\n");
	//printMat(_R, "matrix R\n");
	Matrix QR = createMat(_Q.rows, _R.cols);
	QR = multMat(_Q, _R);
	//printMat(QR, "matrix A = QR\n");
	Matrix RQ = createMat(_Q.rows, _R.cols);
	RQ = multMat(_R, _Q);
	//printMat(RQ, "matrix A = RQ\n");
	return RQ;
}
//eigenvalue diagonal matrix 찾는 함수
void eigenvalue(Matrix _A, Matrix _Q, Matrix _R) {
	Matrix A1 = createMat(_A.rows, _A.cols);
	Matrix RQ = createMat(_A.rows, _A.cols);
	Matrix Aeigenvalue = createMat(_A.rows, _A.cols);
	initMat(A1, 0);
	int k = 0;
	double num = 0;
	double Nmax = 100;
	double lamda = 0;

	for (int i = 0; i < _A.rows; i++) {          //_A 값을 mat A1 에 Copy
		for (int j = 0; j < _A.rows; j++) {
			A1.at[i][j] = _A.at[i][j];
		}
	}
	do
	{
		for (int i = 0; i < _A.rows; i++) {          //A1 값을 _A 에 Copy
			for (int j = 0; j < _A.rows; j++) {
				_A.at[i][j] = A1.at[i][j];
			}
		}
		RQ = QRdecomp(_A, _Q, _R);
		for (int i = 0; i < _A.rows; i++) {          //trans 값을 mat A1 에 Copy  A1 = RQ;
			for (int j = 0; j < _A.rows; j++) {
				A1.at[i][j] = RQ.at[i][j];
			}
		}
		k++;
	} while (k < Nmax);
	Aeigenvalue = A1;
	printMat(Aeigenvalue, "eigenvalue matrix A = \n");
	for (int i = 0; i < Aeigenvalue.rows; i++)
	{
		lamda = Aeigenvalue.at[i][i];
		printf(" eigenvalue of matrix A is %f \n", lamda);
	}
}


double norm(Matrix _vecC) {  // norm 값 계산 함수
	double sum = 0;
	for (int i = 0; i < _vecC.rows; i++)
	{
		sum += _vecC.at[i][0] * _vecC.at[i][0];
	}
	sum = sqrt(sum);
	return sum;
}

Matrix transpose(Matrix _vecU) {  // transpose	
	Matrix out = createMat(_vecU.cols, _vecU.rows);

	for (int i = 0; i < _vecU.cols; i++)
	{
		for (int j = 0; j < _vecU.rows; j++)
		{
			out.at[i][j] = _vecU.at[j][i];
		}
	}
	return out;
}
//--------------------------Cond 함수-------------------------

double Cond(Matrix _A) {
	// Matrix at = createMat(A.cols, A.rows);
	Matrix matAT = createMat(_A.rows, _A.cols); // mat A 의 transpose 행렬
	Matrix mulA = createMat(_A.rows, _A.cols);  // mat AT 곱 mat A
	double eigenvalmax, eigenvalnext, eigenvalmin, condnum = 0;

	if (_A.rows  != _A.cols)
	{
		printf(" Input matrix is not square matrix.  \n");
		printf(" Calculate condition number by finding sigular value  \n");
	}
	matAT = transpose(_A);
	printMat(_A, "matA \m");
	printMat(matAT, "matAT \m");
	mulA = multMat(matAT, _A);
	printMat(mulA, "matAT * mat A\n");

	Matrix matQ = createMat(mulA.rows, mulA.rows);
	Matrix matR = createMat(mulA.rows, mulA.cols);
	
	eigenvalue(mulA, matQ, matR);

		eigenvalmax =fabs(mulA.at[0][0]);
		//printf("eigenvalmax = %f\n", eigenvalmax);
		for (int i = 1; i < mulA.cols; i++)
		{
			eigenvalnext = fabs(mulA.at[i][i]);
			//printf("eigenvalnext = %f\n", eigenvalnext);
			if (eigenvalnext > eigenvalmax) {
				eigenvalmax = eigenvalnext;
			}
		}
		eigenvalmin= fabs(mulA.at[0][0]);

		for (int i = 1; i < mulA.cols; i++)
		{
			eigenvalnext = fabs(mulA.at[i][i]);
			//printf("eigenvalnext = %f\n", eigenvalnext);
			if ((eigenvalnext !=0) && (eigenvalnext < eigenvalmin)) { // 0 을 최소값으로 가지면 connum 은 발산함으로 0 을 제외한 최소값
				eigenvalmin = eigenvalnext;
			}
		}
		eigenvalmax = sqrt(eigenvalmax);
		eigenvalmin = sqrt(eigenvalmin);
		//printf("eigenvalmax = %f\n", eigenvalmax);
		//printf("eigenvalmin = %f\n", eigenvalmin);
		condnum = eigenvalmax / eigenvalmin;
		//printf("condnum = %f\n", condnum);
		return condnum;
}


//-----------------------------Jacobian matrix-------------------------------------
//---------------------------------------------------------------------------------

//----------------------------첫번째 함수----------------------------
double func1(double _x, double _y) // hybrid function f(x)=(1/x)-2
{
	double e = 2.718282;
	double F = _y + (pow(e,_x/2)- pow(e, _x / -2))/-2;
	return F;
}

double funcx1d(double _x, double _y) // hybrid function f(x)=(1/x)-2
{
	/*double small = pow(10, -2);
	double F = (func1(_x,0) - func1(_x-small,0)) / small;*/
	double e = 2.718282;
	double F = (pow(e, _x / 2) - pow(e, _x / -2)) / -4;
	return F;
}

double funcy1d(double _x, double _y) // hybrid function f(x)=(1/x)-2
{
	double small = pow(10, -10);
	double F = (func1(0, _y) -func1(0, _y - small)) / small;
	return F;
}
//------------------------두번째 함수 ------------------------
double func2(double _x, double _y) // hybrid function f(x)=(1/x)-2
{
	double F = 9 * pow(_x, 2) + 25 * pow(_y, 2) - 225;
	return F;
}

double funcx2d(double _x, double _y) // hybrid function f(x)=(1/x)-2
{
	double small = pow(10, -10);
	double F = (func2(_x, 0) - func2(_x - small, 0)) / small;
	return F;
}

double funcy2d(double _x, double _y) // hybrid function f(x)=(1/x)-2
{
	double small = pow(10, -10);
	double F = (func2(0, _y) - func2(0, _y - small)) / small;
	return F;
}

Matrix Fmatrix(double _x, double _y) {
	Matrix funcF = createMat(2, 1);
	double a, b= 0;
	a = func1(_x, _y);
	b = func2(_x, _y);
	funcF.at[0][0] = a;
	funcF.at[1][0] = b;
	return funcF;

}

Matrix Fdmatrix(double _x, double _y) {
	Matrix funcdF = createMat(2, 2);
	double a, b,c,d = 0;
	a = funcx1d(_x, _y);
	b = funcy1d(_x, _y);
	c = funcx2d(_x, _y);
	d = funcy2d(_x, _y);
	funcdF.at[0][0] = a;
	funcdF.at[0][1] = b;
	funcdF.at[1][0] = c;
	funcdF.at[1][1] = d;

	return funcdF;

}

Matrix Jacobian(double _x0,double _y0) // // definitnion of newtonRaphson
{
	Matrix vecX = createMat(2, 1);
	Matrix matXn = createMat(2, 1);
	Matrix matH = createMat(2, 1);
	Matrix matJ = createMat(2, 2);
	Matrix matJU = createMat(2, 2);
	Matrix matF = createMat(2, 1);
	

	initMat(matH, 0);
	initMat(matJ, 0);

	vecX.at[0][0] = _x0;
	vecX.at[1][0] = _y0;
	//printMat(vecX, "vecX\n");

	double a, b,f1,f2 = 0;
	int Nmax = 100;
	int k = 0;
	double ep1 ,ep2= 0;
	double tol = 0.0000001;
	matF = Fmatrix(_x0, _y0);
	matJ = Fdmatrix(_x0, _y0);
	printMat(matJ, "matJ \n");
	do
	{
		//printf(" k = %d\n", k);
	double num = 0;
	for (int i = 0; i < matF.rows; i++)
	{
		for (int j = 0; j < matF.cols; j++)
		{
			num = -1 * matF.at[i][j];
			matF.at[i][j] = num;
		}

	}
	//printMat(matF, "matF\n");
	gaussEilmJordan(matJ, matF, matJU, matH);  //가우스 조던 함수를 통해서 H를 구한다.
	//printMat(matH, "matH \n");
	matXn = addMat(vecX ,matH);

	//printMat(matXn, "matX1 \n");
	a = matXn.at[0][0];
	b = matXn.at[1][0];
	//XN+1백터값 저장
	vecX.at[0][0] = a;
	vecX.at[1][0] = b;

	//새로 찾은 XN+1의 값 (a,b) 를 통해서 새로운 F J matrix 생성
	matF = Fmatrix(a, b);
	matJ = Fdmatrix(a, b);
	f1 = func1(a, b);
	f2 = func2(a, b);
	//printMat(matF, "matF == \n");
	k++;
	ep1=fabs(func1(a,b));
	ep2 = fabs(func2(a, b));
	//printf(" abs value of f1 at a,b = %f\n", ep1);
	//printf(" abs value of f2 at a,b = %f\n", ep2);
	//f1 과 f2 에 Xn+1을 통해 찾은 (a,b)의 값을 넣었을때 함숫값이 tol 보다 작아질때까지 반복
	} while (k<Nmax && (ep1>tol || ep2>tol)); 
	printf(" abs value of f1 at a,b = %f\n", ep1);
	printf(" abs value of f2 at a,b = %f\n", ep2);
	printf("Total Iteration: % d \n", k);
	return matXn;
}


//-------------BISECTION AND NEWTON-------------------------//
double funch(double _x) // hybrid function f(x)=(1/x)-2
{
	double F = pow(_x, -1) - 2;
	return F;
}
double funchd(double _x) // hybrid function f'(x)
{
	double F = -1 * pow(_x, -2);
	return F;
}
double funcd(double _x) // non-linear eqution first derivative
{
	//double  w0 = 20 * pow(10, 3), E = 70 * pow(10, 9), L = 4, I = 52.9 * pow(10, -6);
	double  F = _x * _x * _x;
	//double F = w0 * (7 * pow(L, 4) - 30 * pow(L, 2) * pow(_x, 2) + 15 * pow(_x, 4)) / (360 * L * E * I);
	return F;
}
double func2d(double _x) // non-linear eqution second derivative
{
	//double  w0 = 20 * pow(10, 3), E = 70 * pow(10, 9), L = 4, I = 52.9 * pow(10, -6);
	double  F = 3 * _x * _x;
	//double F = w0 * (-60 * pow(L, 2) * _x + 60 * pow(_x, 3)) / (360 * L * E * I);
	return F;
}
double bisectionNL(double _a0, double _b0, double _tol) // definitnion of bisectionNL 
{
	int k = 0;
	int Nmax = 30;
	double a = _a0;
	double b = _b0;
	double xn = 0;
	double ep = 1000;


	do {
		xn = (a + b) / 2;
		ep = fabs(funcd(xn));
		printf("Iteration:%d \t", k);
		printf("X(n): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);

		if (funcd(a) * funcd(xn) < 0) // check solution between  a and xn
			b = xn;
		else
			a = xn;
		k++;
	} while (k<Nmax && ep>_tol);


	return xn;
}

double newtonRaphson(double _x0, double _tol) // // definitnion of newtonRaphson
{
	int k = 0;
	int Nmax = 30;
	double x0 = _x0;
	double xn = 0;
	double ep = 0;


	do {
		xn = x0 - (funcd(x0) / func2d(x0));
		 ep = fabs(xn - x0);
		//ep = fabs(funcd(xn));
		printf("Iteration:%d \t", k);
		printf("X(n): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);

		k++;
		x0 = xn;
	} while (k<Nmax && ep>_tol);

	return xn;
}
double newtonRaphson_hybrid(float _x0, float _tol) // f(x)=(1/x)-2 일 때 newtonRaphson 
{
	int k = 0;
	int Nmax = 7;
	float x0 = _x0;
	float xn = 0;
	float ep = 1000;


	do {
		xn = x0 - (funch(x0) / funchd(x0));
		ep = fabs(xn - x0);
		printf("Iteration:%d \t", k);
		printf("X(n): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);

		k++;
		x0 = xn;
	} while (k<Nmax && ep>_tol);

	return xn;
}

double hybridmethod(float _x0, float _tol) // definition of hybridmethod
{
	int k = 0;
	int Nmax = 20;
	float xn = 0;
	double x0 = _x0;
	float a = 0.01;
	float b = 0.99;
	float c = 0;
	double ep = 1000;
	float tol = 0.00001;
	double hybridresult;

	do {
		if (x0 - (funch(x0) / funchd(x0)) > a && x0 - (funch(x0) / funchd(x0)) < b) // xn 값이 boundadry [ a b ]에 존재하는 지 확인  
		{                                                                   // 존재 한다면 newton-Raphson method 로 해를 구함    
			xn = x0 - (funch(x0) / funchd(x0));
			ep = fabs(xn - x0);
			printf("newtonRaphson\t");
			printf("Iteration:%d \t", k);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f\n", ep);

			x0 = xn;
		}
		else {                                                                  // xn 값이 boundadry [ a b ]에 하지 않는다면
			c = x0;                                                             // bisection-method 을 통해서 xn값이 [a b]에 존재 할 때까지 xn을 구함
			xn = (a + c) / 2;
			ep = fabs(funch(xn));
			printf("bisection\t");
			printf("Iterationd:%d \t", k);
			printf("X(n): %f \t", xn);
			printf("Tolerance: %.10f\n", ep);
			if (funch(a) * funch(xn) < 0)
				x0 = xn;
			else
				a = xn;
		}
		k++;
	} while (k<Nmax && ep>_tol);
	return xn;
}


//-----------------------Curve Fitting and Interpolation-----------------------//
void linearFit(Matrix _x, Matrix _y) {

	int mx = _x.rows;
	int my = _y.rows;
	double a1 = 0;
	double a0 = 0;
	double Sx = 0;
	double Sy = 0;
	double Sxx = 0;
	double Sxy = 0;
	double P = 0;
	double x = 100;
	// x 와 y 의 정보 data 수가 다르면 오류 검출문 작성
	// x 와 y 의 정보 data 수는 최소 2개 이상
	if ((mx == my)) {
		if ((mx >= 2) && (my >= 2)) {
			for (int i = 0; i < mx; i++)
			{
				Sx += _x.at[i][0];
				Sy += _y.at[i][0];
				Sxx += pow(_x.at[i][0], 2);
				Sxy += _x.at[i][0] * _y.at[i][0];
			}
			a0 = (Sxx * Sy - Sxy * Sx) / (mx * Sxx - Sx * Sx);
			a1 = (mx * Sxy - Sx * Sy) / (mx * Sxx - Sx * Sx);
			P = a1 * x + a0;
			//Matlab problem2 y=b*exp(mx) function
			//-------------------------------//
			double m = 0;
			double b = 0;
			m = a1;
			b = exp(a0);
			printf("m= %f \t b= %f \n",m,b);
			//------------------------------//
			double z_array[] = { a1, a0 };
			printf("f(x) = %f + %fx \n", a0, a1);
			printf("Pressure at T= %f, P = %f\n", x, P);
		}
		else
		{
			printf("the number of dataset is less than 2.\n");
		}
	}
	else
	{
		printf("The number of element in x must be the same as in y.\n");
	}
}

void planeFit(Matrix _x, Matrix _y, Matrix _z) {

	Matrix Plane = createMat(3, 3);
	Matrix invPlane = createMat(3, 3);
	Matrix vecy = zeros(3, 1);
	Matrix matq = createMat(3, 1);
	initMat(Plane, 0);
	initMat(vecy, 0);
	initMat(matq, 0);
	int mx = _x.rows;
	int my = _y.rows;
	double p0 = 0;
	double p1 = 0;
	double p2 = 0;
	double Sx = 0;
	double Sz = 0;
	double Sxx = 0;
	double Sxz = 0;
	double Szz = 0;
	double Sy = 0;
	double Sxy = 0;
	double Szy = 0;
	// x 와 y 의 정보 data 수가 다르면 오류 검출문 작성
	// x 와 y 의 정보 data 수는 최소 2개 이상
	if ((mx == my)) {
		if ((mx >= 2) && (my >= 2)) {
			for (int i = 0; i < mx; i++)
			{
				Sx += _x.at[i][0];
				Sz += _z.at[i][0];
				Sxx += pow(_x.at[i][0], 2);
				Sxz += _x.at[i][0] * _z.at[i][0];
				Szz += pow(_z.at[i][0], 2);
				Sy += _y.at[i][0];
				Sxy += _x.at[i][0] * _y.at[i][0];
				Szy += _z.at[i][0] * _y.at[i][0];
			}
		}
		else
		{
			printf("the number of dataset is less than 2.\n");
		}
	}
	else
	{
		printf("The number of element in x must be the same as in y.\n");
	}
	Plane.at[0][0] = mx;
	Plane.at[0][1] = Sx;
	Plane.at[0][2] = Sz;
	Plane.at[1][0] = Sx;
	Plane.at[1][1] = Sxx;
	Plane.at[1][2] = Sxz;
	Plane.at[2][0] = Sz;
	Plane.at[2][1] = Sxz;
	Plane.at[2][2] = Szz;

	vecy.at[0][0] = Sy;
	vecy.at[1][0] = Sxy;
	vecy.at[2][0] = Szy;

	inv(Plane, invPlane);
	//printMat(vecy, "vecy\n");
	matq = multMat(invPlane, vecy);
	//printMat(matq, "matq\n");
	
	p0 = matq.at[0][0];
	p1 = matq.at[1][0];
	p2 = matq.at[2][0];
	printf("Plane equation %fx + %fy + %fz %f = 0 \n", -1*p1, 1.0, -1*p2, -1*p0);
	double alpha = 0;
	double Beta = 0;
	double gamma = 0;
	double aiming_angle = 0;
	alpha = atan(-1*p2)*180/PI;
	Beta = atan(p1/ p2) * 180 / PI;
	gamma = atan(-1 / p1) * 180 / PI;
	aiming_angle = 90 - gamma;
	printf("alpha = %f\nBeta = %f\ngamma = %f\n", alpha,Beta,gamma);
	printf("aiming_angle = %f\n", aiming_angle);
	double vecnorm_array[200] = { -1*p1, 1,-1*p2 };
	double normvalue = 0;
	Matrix Vecnorm = arr2Mat(vecnorm_array, 3, 1);
	Matrix Normalvec = createMat(3, 1);
	
	normvalue = norm(Vecnorm);
	
	for (int i = 0; i < Vecnorm.rows; i++)
	{
		Normalvec.at[i][0] = Vecnorm.at[i][0] / normvalue;
	}
	printMat(Normalvec, "Normal vector \n");
}


void CurvelinearHO(Matrix _x, Matrix _y, Matrix _xq,int N) {

	int mx = _x.rows; // 행 
	int my = N;		// 열 차수
	//printf("mx = %d\n", mx);
	Matrix A = createMat(mx, my);
	Matrix TA = createMat(my, mx);
	Matrix ATA = createMat(mx, mx);
	Matrix invATA = createMat(mx, mx);

	Matrix TAT = createMat(mx, mx);
	initMat(A, 0);
	initMat(ATA, 0);
	Matrix b = createMat(my, 1);
	Matrix z = createMat(my, 1);

	
	// x 와 y 의 정보 data 수가 다르면 오류 검출문 작성
	//Matrix A
	for (int k = 0; k < mx; k++)
	{
		double Sxk = 0;
		for (int j = 0; j < my; j++)
		{
			Sxk = pow(_x.at[k][0], j);
			A.at[k][j] = Sxk;
		}
	}
	printMat(A, "MATRIX A\n");
	for (int k = 0; k < my; k++)
	{
		double Sxi = 0;
		for (int j = 0; j < _y.rows; j++)
		{
			printf("_y.at[%d][0]= %f\n", j, _y.at[j][0]);
			Sxi += _y.at[j][0]*pow(_x.at[j][0], k);
			printf("Sxi= %f\n", Sxi);
		}
		b.at[k][0] = Sxi;
		printf("b.at[%d][0]= %f\n", k, b.at[k][0]);
	}
	TA = transpose(A);
	printMat(TA, "MATRIX TA\n");
	ATA = multMat(TA, A);
	//invATA = inv(ATA, invATA);
	printMat(invATA, "AVA\n");
	//b = multMat(TA, _y);
	printMat(b, "b\n");
	z = multMat(invATA, b);
	printMat(z, "z\n"); // 첫번째 결과부터 a0 a1 a2 ....an
		/*for (int k = 0; k < mx; k++) //행 증가
		{
			double Sxy = 0;
			for (int j = 0; j < my; j++)  //각각의 열 증가
			{
				Sxy += _y.at[j][0] * pow(_x.at[j][0], k);
				double Sx = 0;
				for (int a = 0; a < mx; a++)  //각 열에서의 Sum
				{
					Sx += pow(_x.at[a][0], j + k);
					//printf("j+k= %d\n", j + k);
					//printf("Sx= %f\n", Sx);
				}
				ATA.at[k][j] = Sx;
			}
			b.at[k][0] = Sxy;
		}*/
	/*double Sxy = 0;
	for (int i = 0; i < my; i++)
	{
		for (int j = 0; j < my; j++)
		{
			Sxy = _y.at[j][0] * pow(_x.at[j][0], i);
		}
		b.at[i][0] = Sxy;
	}
	printMat(_y, "y\n");
	printMat(b, "b\n");*/
	//printMat(_y, "y\n");
	//printMat(b, "b\n");
	//printMat(ATA, "ATA\n");
	//Matrix y = createMat(21, 1);
	//invATA = inv(ATA, invATA);
	//printMat(invATA, "invATA\n");
	//printMat(b, "b\n");
	//b = multMat(TA, _y);
	//printMat(_y, "y\n");
	//z = multMat(invATA, b);
	//printMat(z, "z\n");

/*	int j = 0;
for (int k = 0; k <21; k ++)
	//for (int k = 0; k <= 100; k+=5)
	{
		double sum = 0;
		for (int i = 0; i < mx; i++)
		{
			//printf("k = %d\n", k);
			//printf("z.at[i][0] = %f\n", z.at[i][0]);
			//printf("_xq.at[%d][0]=%f\n",k,_xq.at[k][0]);
			//printf("pow(k,i)= %f\n", pow(_xq.at[k][0], i));
			//printf("pow(k,i)= %f\n", pow(k, i));
			//printf("a = %f\n", a);
			sum += (z.at[i][0]) * pow(_xq.at[k][0], i);
			//sum += (z.at[i][0]) * pow(k, i);
			//printf("sum 과정  = %f\n", sum);
		}
		//printf("sumout = %f\n", sum);
		y.at[j][0] = sum;
		j++;
		//printf("sumout 222= %f\n", sum);
	}
	//printMat(y, "y\n");*/
}

void linearInterp(Matrix _x, Matrix _y, Matrix _xq) {

	//10개의 구간에 대한 각각의 a0, a1 의 값을 저장할 백터
	Matrix fx = createMat(2, _x.rows - 1);
	//return 되는 yq 
	Matrix yq = createMat(_xq.rows, 1);
	initMat(fx, 0);
	initMat(yq, 0);
	double a1, a0, y, output = 0;
	int mx = _x.rows;
	int my = _y.rows;
	//11개의 dataset 에 대한 10 개의 a0, a1 의 값을 fx 에 저장
	for (int i = 0; i < _x.rows - 1; i++)
	{
		a0 = ((_y.at[i + 1][0] * _x.at[i][0]) - (_y.at[i][0] * _x.at[i + 1][0])) / (_x.at[i][0] - _x.at[i + 1][0]);
		a1 = (_y.at[i][0] - _y.at[i + 1][0]) / (_x.at[i][0] - _x.at[i + 1][0]);
		fx.at[0][i] = a0;
		fx.at[1][i] = a1;
	}
	if ((mx == my)) {
		if ((mx >= 2) && (my >= 2)) {
			//xq 의 각 data 값이 해당되는 구간의 line에 대입하여 yq 에 함숫값 저장
			for (int i = 0; i < _xq.rows; i++)
			{
				for (int j = 0; j < _x.rows - 1; j++)
				{
					if ((_xq.at[i][0] >= _x.at[j][0]) && (_xq.at[i][0] <= _x.at[j + 1][0]))
					{
						y = fx.at[1][j] * _xq.at[i][0] + fx.at[0][j];
						yq.at[i][0] = y;
						if (_xq.at[i][0] == 75) {
							printf("Estimate the interppolated value P when T = %f is %f\n", _xq.at[i][0], yq.at[i][0]);
						}
						break;
					}
				}
			}
			printMat(yq, "yq \n");
		}
		else
		{
			printf("the number of dataset is less than 2.\n");
		}
	}
	else
	{
		printf("The number of element in x must be the same as in y.\n");
	}
}

void InterpolateLinear(Matrix _x, Matrix _y, Matrix _z) {

	//구간에 대한 각각의 a0, a1 의 값을 저장할 백터
	Matrix fyz = createMat(2, _x.rows - 1);
	//구간에 들어가는 query
	Matrix yzq = createMat(60, 1);
	for (int i = 0; i < _x.rows-30; i++)
	{
		yzq.at[i][0] = (_y.at[i][0] + _y.at[i + 30][0]) / 2; //결국 yzq 넣었을때 나오는 값은 z좌표
	}
	// return 되는 값
	Matrix x1 = createMat(60, 1);
	Matrix y1 = createMat(60, 1);
	Matrix z1 = createMat(60, 1);
	initMat(fyz, 0);
	initMat(z1, 0);
	double a1yz,a0yz=0;
	double yz = 0;
	int mx = _x.rows;
	int my = _y.rows;
	//dataset 에 대한 각각의 a0, a1 의 값을 fx 에 저장
	for (int i = 0; i < _x.rows -30; i++)
	{
		a0yz = ((_z.at[i + 30][0] * _y.at[i][0]) - (_z.at[i][0] * _y.at[i + 30][0])) / (_y.at[i][0] - _y.at[i + 30][0]);
		a1yz = (_z.at[i][0] - _z.at[i + 30][0]) / (_y.at[i][0] - _y.at[i + 30][0]);
		fyz.at[0][i] = a0yz;
		fyz.at[1][i] = a1yz;
	}
	if ((mx == my)) {
		if ((mx >= 2) && (my >= 2)) {
			for (int i = 0; i < yzq.rows; i++) {
				yz = fyz.at[1][i] * yzq.at[i][0] + fyz.at[0][i];
				z1.at[i][0] = yz;
			}
			for (int i = 0; i < 60; i++)
				printf("x1= %f\t y1= %f\t z1= %f\n", _x.at[i][0],yzq.at[i][0], z1.at[i][0]);
			printf("\n");
		}
		else
		{
			printf("the number of dataset is less than 2.\n");
		}
	}
	else
	{
		printf("The number of element in x must be the same as in y.\n");
	}
}

void linearInterp2nd(Matrix _x, Matrix _y, Matrix _xq) {

	//10개의 구간에 대한 각각의 a0, a1 의 값을 저장할 백터
	Matrix fx = createMat(3, _x.rows - 2);
	//return 되는 yq 
	Matrix yq = createMat(_xq.rows, 1);
	initMat(fx, 0);
	initMat(yq, 0);
	double a2, a1, a0, y, output = 0;
	int mx = _x.rows;
	int my = _y.rows;
	//11개의 dataset 에 대한 10 개의 a0, a1 의 값을 fx 에 저장
	for (int i = 0; i < _x.rows - 2; i++)
	{
		a0 = _y.at[i][0] / ((_x.at[i][0] - _x.at[i + 1][0])* (_x.at[i][0] - _x.at[i + 2][0]));
		a1 = _y.at[i+1][0] / ((_x.at[i+1][0] - _x.at[i][0]) * (_x.at[i+1][0] - _x.at[i + 2][0]));
		a2 = _y.at[i+2][0] / ((_x.at[i+2][0] - _x.at[i][0]) * (_x.at[i+2][0] - _x.at[i + 1][0]));
		fx.at[0][i] = a0;
		fx.at[1][i] = a1;
		fx.at[2][i] = a2;
	}
	//printMat(fx, "fx\n");
	if ((mx == my)) {
		if ((mx >= 2) && (my >= 2)) {
			//xq 의 각 data 값이 해당되는 구간의 line에 대입하여 yq 에 함숫값 저장
			for (int i = 0; i < _xq.rows; i++)
			{
				for (int j = 0; j < _x.rows - 2; j++)
				{
					if ((_xq.at[i][0] >= _x.at[j][0]) && (_xq.at[i][0] <= _x.at[j + 1][0]))
					{
						y = fx.at[0][j] * ((_xq.at[i][0] - _x.at[j + 1][0]) * (_xq.at[i][0] - _x.at[j + 2][0])) + fx.at[1][j] * ((_xq.at[i][0] - _x.at[j][0]) * (_xq.at[i][0] - _x.at[j + 2][0])) + fx.at[2][j] * ((_xq.at[i][0] - _x.at[j][0]) * (_xq.at[i][0] - _x.at[j + 1][0]));
						yq.at[i][0] = y;
						if (_xq.at[i][0] == 75) {
							printf("Estimate the interppolated value P when T = %f is %f\n", _xq.at[i][0], yq.at[i][0]);
						}
						break;
					}
				}
			}
			printMat(yq, "yq \n");
		}
		else
		{
			printf("the number of dataset is less than 2.\n");
		}
	}
	else
	{
		printf("The number of element in x must be the same as in y.\n");
	}
}

// Create a matrix from 1D-array
Matrix arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}	


//------------------Gradient Fuction---------------------//

// Return the dy/dx results for the input data. (truncation error: O(h^2))
Matrix	gradient(Matrix _x, Matrix _y) {
	int mx = _x.rows;
	Matrix funcd = createMat(_x.rows, 1);
	initMat(funcd, 0);
	double h = _x.at[1][0] - _x.at[0][0];

	// For 2 point of data, use 2-point forward or backward
	if ((_x.rows == 2) && (_y.rows == 2)) {
		funcd.at[0][0] = (_y.at[1][0] - _y.at[0][0]) / h;
		funcd.at[1][0] = (_y.at[1][0] - _y.at[0][0]) / h;
	}
	else {
		funcd.at[0][0] = (-3 * _y.at[0][0] + 4 * _y.at[1][0] - _y.at[2][0]) / (2 * h);  // Three-point forward 
		funcd.at[mx - 1][0] = (_y.at[mx - 3][0] - 4 * _y.at[mx - 2][0] + 3 * _y.at[mx - 1][0]) / (2 * h); // Three-point backward
		for (int i = 1; i < _x.rows - 1; i++)
		{
			funcd.at[i][0] = (_y.at[i + 1][0] - _y.at[i - 1][0]) / (2 * h);
		}
	}
	return funcd;
}


void	gradient1D(double x[], double y[], double dydx[], int m) {
	int mx = m;
	double h = x[1] - x[0];
	// For 2 point of data, use 2-point forward or backward
	if (m == 2) {
		dydx[0] = (y[1] - y[0]) / h;
		dydx[1] = (y[1] - y[0]) / h;
	}
	else
	{
		dydx[0] = (-3 * y[0] + 4 * y[1] - y[2]) / (2 * h);  //Three-point forward
		dydx[mx - 1] = (y[mx - 3] - 4 * y[mx - 2] + 3 * y[mx - 1]) / (2 * h);  // Three-point backward
		for (int i = 1; i < mx - 1; i++)
		{
			dydx[i] = (y[i + 1] - y[i - 1]) / (2 * h);
		}
	}
}



// Return the dy/dx results for the target equation. (truncation error: O(h^2))
Matrix	gradientFunc(double func(const double x), Matrix xin) {

	int n = xin.rows;
	Matrix y = createMat(n, 1);
	Matrix df = createMat(n, 1);
	for (int i = 0; i < n; i++)
	{
		y.at[i][0] = func(xin.at[i][0]);
	}
	df = gradient(xin, y);
	return df;
}

double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), double _x0, double _tol)
{
	int k = 0;
	int Nmax = 38;
	double x0 = _x0;
	double xn = 0;
	double ep = 0;


	do {
		xn = x0 - (func(x0) / dfunc(x0));
		ep = fabs(func(xn));
		printf("Iteration:%d \t", k);
		printf("X(n): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);

		k++;
		x0 = xn;
	} while (k<Nmax && ep>_tol);

	return xn;
}

//-----------------Integral Function---------------//
double trapz(double x[], double y[], int m) {
	int mx = m;
	double sum = 0;
	double h = x[1] - x[0];

	for (int i = 1; i < mx - 1; i++)
	{
		sum += y[i];
	}
	sum = (h / 2) * (y[0] + y[m - 1]) + h * sum;

	return sum;
}

double IntegrateRect(double _x[], double _y[], int _m) {
	int N = _m - 1;
	double I = 0;
	for (int i = 0; i < N; i++)
		I += _y[i] * (_x[i + 1] - _x[i]);

	return I;
}
//Simpson 1/3
double integral(double func(const double x), double a, double b, int n) { 
	double sum = 0;
	double	term1 = 0;
	double	term2 = 0;
	double h = (b - a) / n;
	double x0 = a;
	double x1 = a + h;
	for (int i = 1; i < n; i += 2)
	{
		double xi = x0 + h * i;
		term1 += func(xi);
	}
	for (int i = 1; i < n-1; i += 2)
	{
		double xj = x1 + h * i;
		term2 += func(xj);
	}
	sum = (h / 3) * (func(a) + 4 * term1 + 2 * term2 + func(b));
	return sum;
}

double integralMid(double x[], double y[], int m) {
	int mx = m;
	double sum = 0;
	double h = x[1] - x[0];

	for (int i = 0; i < mx - 1; i++)
	{
		sum += y[i] + (y[i + 1] - y[i]) / 2;
	}
	sum = h * sum;
	return sum;
}
//Simpson 3/8
double integral38(double func(const double x), double a, double b, int n) {
	double sum = 0;
	double	term1 = 0;
	double	term2 = 0;
	double h = (b - a) / n;
	double x0 = a;
	double x1 = a + 2 * h;
	for (int i = 1; i < n-1; i += 3)
	{
		double xi = x0 + h * i;
		//double xj = x1 + h * i;
		term1 += func(xi) + func(xi + h);
		//term2 += func(xj);
	}
	for (int i = 1; i < n-2; i += 3)
	{
		double xj = x1 + h * i;
		term2 += func(xj);
	}
	sum = (3 * h / 8) * (func(a) + 3 * term1 + 2 * term2 + func(b));

	return sum;
}

//------------ODE Method--------------//
void odeEU(double func(const double x, const double y), double y[], double t0, double tf, double h, double y0) {

	double a = t0;
	double b = tf;
	double y_E = 0;
	int N = (tf - t0) / h + 1;
	y[0] = y0;
	for (int i = 0; i < N - 1; i++)
	{
		y[i + 1] = y[i] + func(t0+h*i, y[i]) * h;
	}
}

void odeEM(double func(const double x, const double y), double y[], double t0, double tf, double h, double y0) {

	double a = t0;
	double b = tf;
	double y_E = 0;
	int N = (tf - t0) / h + 1;
	y[0] = y0;
	for (int i = 0; i < N-1; i ++)
	{
		y_E = y[i] + func(t0+h*i, y[i]) * h;
		y[i + 1] = y[i] + (func(t0+h*i, y[i]) + func((t0+h)+h*i, y_E)) * h / 2;
	}
}
void odePC(double func(const double x, const double y), double y[], double t0, double tf, double h, double y0) {

	double a = t0;
	double b = tf;
	double y_E = 0;
	double y_E1 = 0;
	int N = (tf - t0) / h + 1;
	y[0] = y0;
	double ep = 0;
	double tol = pow(10, -5);
	int i = 0;

	do {
		y_E = y[i] + func(t0 + h * i, y[i]) * h;
	    y_E1= y[i] + (func(t0 + h * i, y[i]) + func((t0 + h) + h * i, y_E)) * h / 2;
		y[i + 1] = y[i] + (func(t0 + h * i, y[i]) + func((t0 + h) + h * i, y_E1)) * h / 2;
		ep = fabs((y[i + 1] - y[i]) / y[i]);
		i++;
	} while (i<N && ep>tol);
}

void ode(double func(const double x, const double y), double y[], double t0, double tf, double h, double y0,int method) {

	if (method == 0)
	{
		odeEU(func, y, t0, tf, h,y0);
	}
	else if (method == 1)
	{
		odeEM(func, y, t0, tf, h,y0);
	}
}


//--------------------------------RK2 RK4 SYSRK2 SYSRK4-----------------------------//
void odeRK2(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{
	double C1 = 0.5;
	double C2 = 0.5;
	double alpha = 1;
	double beta = alpha;  // alpha=beta

	int N = (tf - t0) / h + 1;
	double ti = t0;
	double y_EU;
	double K1 = 0, K2 = 0;

	// Initialization 
	y[0] = y0;

	for (int i = 0; i < N - 1; i++)
	{

		K1 = odeFunc(ti, y[i]);
		K2 = odeFunc(ti + alpha * h, y[i] + beta * K1 * h);
		y[i + 1] = y[i] + (C1 * K1 + C2 * K2) * h;
		ti += h;
	}
}

void odeRK3(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{

	double alpha2 = 0.5, alpha3 = 1, Betha21 = 0.5, Betha31 = -1, Betha32 = 2;
	int N = (tf - t0) / h + 1;
	double ti = t0;
	double y_EU;
	double K1 = 0, K2 = 0, K3 = 0;

	// Initialization 
	y[0] = y0;

	for (int i = 0; i < N - 1; i++)
	{

		K1 = odeFunc(ti, y[i]);
		K2 = odeFunc(ti + alpha2 * h, y[i] + Betha21 * K1 * h);
		K3 = odeFunc(ti + alpha3 * h, y[i] + Betha31 * K1 * h + Betha32 * K2 * h);

		y[i + 1] = y[i] + (K1 + 4 * K2 +  K3 ) * h / 6;
		ti += h;
	}
}

void odeRK4(double odeFunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{

	double a = 0.5;
	int N = (tf - t0) / h + 1;
	double ti = t0;
	double y_EU;
	double K1 = 0, K2 = 0, K3 = 0, K4 = 0;

	// Initialization 
	y[0] = y0;

	for (int i = 0; i < N - 1; i++)
	{

		K1 = odeFunc(ti, y[i]);
		K2 = odeFunc(ti + a * h, y[i] + a * K1 * h);
		K3 = odeFunc(ti + a * h, y[i] + a * K2 * h);
		K4 = odeFunc(ti + h, y[i] + K3 * h);

		y[i + 1] = y[i] + (K1 + 2 * K2 + 2 * K3 + K4) * h / 6;
		ti += h;
	}
}


void sys2RK2(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double z1[], double t0, double tf, double h, double y1_init, double z1_init)
{
	int N = (tf - t0) / h + 1;
	double ti = 0;

	double K1[2] = { 0 };
	double K2[2] = { 0 };
	double Yin[2] = { 0 };
	double K_Y1 = 0, K_Z1 = 0, K_Y2 = 0, K_Z2 = 0;

	// Initial condition
	y1[0] = y1_init;
	z1[0] = z1_init;

	for (int i = 0; i < N - 1; i++) {

		// Slope 1 : K1
		Yin[0] = y1[i];		// z 즉, y dot 
		Yin[1] = z1[i];		// dzdt = z dot 즉, y 2dot 	
		odeFunc_sys2(ti, Yin, K1);
		K_Y1 = K1[0];
		K_Z1 = K1[1];

		// Slope 2 : K2
		Yin[0] = y1[i] + K_Y1 * h;
		Yin[1] = z1[i] + K_Z1 * h;
		odeFunc_sys2(ti + h, Yin, K2);
		K_Y2 = K2[0];
		K_Z2 = K2[1];

		// Update
		y1[i + 1] = y1[i] + (K_Y1 + K_Y2) * 0.5 * h;
		z1[i + 1] = z1[i] + (K_Z1 + K_Z2) * 0.5 * h;
		ti += h;
	}
}


// Classical RK4
void sys2RK4(void odeFunc_sys2(const double t, const double Y[], double dYdt[]), double y1[], double z1[], double t0, double tf, double h, double y1_init, double z1_init)
{
	int N = (tf - t0) / h + 1;
	double ti = 0;

	double K1[2] = { 0 };
	double K2[2] = { 0 };
	double K3[2] = { 0 };
	double K4[2] = { 0 };
	double Yin[2] = { 0 };
	double K_Y1 = 0, K_Z1 = 0, K_Y2 = 0, K_Z2 = 0, K_Y3 = 0, K_Z3 = 0, K_Y4 = 0, K_Z4= 0;

	// Initial condition
	y1[0] = y1_init;
	z1[0] = z1_init;

	for (int i = 0; i < N - 1; i++) {

		// Slope 1 : K1
		Yin[0] = y1[i];		// z 즉, y dot 
		Yin[1] = z1[i];		// dzdt	 즉, y 2dot 	//initial y and z 
		odeFunc_sys2(ti, Yin, K1);
		K_Y1 = K1[0];		//  First slope of y and z
		K_Z1 = K1[1];

		// Slope 2 : K2
		Yin[0] = y1[i] + 0.5 * K_Y1 * h; //partial update yi+1 and zi+1 
		Yin[1] = z1[i] + 0.5 * K_Z1 * h;
		odeFunc_sys2(ti + 0.5 * h, Yin, K2);  //updated Yin 
		K_Y2 = K2[0];					// Second slope of y and z 
		K_Z2 = K2[1];

		// Slope 3 : K3
		Yin[0] = y1[i] + 0.5 * K_Y2 * h;  //partial update yi+1 and zi+1
		Yin[1] = z1[i] + 0.5 * K_Z2 * h;
		odeFunc_sys2(ti + 0.5 * h, Yin, K3);
		K_Y3 = K3[0];					// Third slope of y and z
		K_Z3 = K3[1];

		// Slope 4 : K4
		Yin[0] = y1[i] + K_Y3 * h;		// 
		Yin[1] = z1[i] + K_Z3 * h;
		odeFunc_sys2(ti + h, Yin, K4);
		K_Y4 = K4[0];
		K_Z4 = K4[1];

		// Update
		y1[i + 1] = y1[i] + (K_Y1 + 2 * K_Y2 + 2 * K_Y3 + K_Y4) * h / 6;
		z1[i + 1] = z1[i] + (K_Z1 + 2 * K_Z2 + 2 * K_Z3 + K_Z4) * h / 6;
		ti += h;
	}
}