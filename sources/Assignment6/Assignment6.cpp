/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Song Yeong Won
Created          : 10-05-2021
Modified         : 15-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Differentiation_student.cpp
-------------------------------------------------------------------------------*/

#include "../../include/myNM.h"

#define Assignment	6		// enter your assignment number
#define eval		0		// set 0

double myFunc(const double x) {
	return  x * x * x;
}
double mydFunc(const double x) {
	return 3 * x * x;
}

int main(int argc, char* argv[])
{
	// PART 1
	printf("\n**************************************************");
	printf("\n|                     PART 1.                    |");
	printf("\n**************************************************\n");

	double t_array[21] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
	double pos_array[21] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59};
	double vel_array[21] = {0};
	double acc_array[21] = {0};
	int m = 21;

	//1)Matirx 
	Matrix t = arr2Mat(t_array, m, 1);
	Matrix x = arr2Mat(pos_array, m, 1);

	Matrix vel = gradient(t, x);
	Matrix acc = gradient(t, vel);

	printMat(t, "t");
	printMat(x, "x");
	printMat(vel, "vel");
	printMat(acc, "acc");

	//2) Using 1D array
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

	//Matlab compare heat flux problem
	double x_array[11] = {0 ,0.01 ,0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1 };
	double T_array[11] = {473, 446.3 ,422.6 ,401.2 ,382 ,364.3 ,348.0 ,332.7, 318.1 ,304.0, 290.1};
	double dTdx_array[11] = { 0 };
	gradient1D(x_array, T_array, dTdx_array, 11);
	double qx0 = 0;
	double qxL = 0;
	qx0 = -240*dTdx_array[0];
	qxL = -240*dTdx_array[10];
	printf("qx0 = %f \tqxL=%f\n", qx0, qxL);
	printf("0 = %f \t 10 = %f\n", dTdx_array[0], dTdx_array[10]);
	
	
// PART 2
	printf("\n**************************************************");
	printf("\n|                     PART 2.                    |");
	printf("\n**************************************************\n");

	Matrix xin = arr2Mat(t_array, m, 1);
	Matrix dydx = gradientFunc(myFunc, xin);
	printMat(xin, "xin");
	printMat(dydx, "dydx");

	double tol = pow(10,-20);
	double x0 = 2;
	double NM_result;

	printf("------------------------------------------------------------------------------------\n");
	printf("         newtonRaphsonFunc method results             \n");
	printf("------------------------------------------------------------------------------------\n");

	printf("newtonRaphson method:\n");
	NM_result = newtonRaphsonFunc(myFunc,mydFunc, x0, tol);
	printf("final solution: %f \t", NM_result);
	printf("\n");

	system("pause");
	return 0;
}