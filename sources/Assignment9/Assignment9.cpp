/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Jan Park, YKKIM
Created          : 2021-06-03
Modified         : 2021-06-03  by YKKIM
Language/ver     : C++ in MSVS2017

Description      : [Tutorial]ODE_IVP_student.c
-------------------------------------------------------------------------------*/

#include "../../include/myNM.h"
#include <stdio.h>
#include <math.h>

#define ODE_EU 0
#define ODE_EM 1
#define ODE_RK2 2
#define ODE_RK4 3

//  PI is defined in  myNM.h
#define PI 3.14159265368979323846264338327950288412

// Problem1: Single equation of 1st order ODE
double odeFunc_rc(const double t, const double v) {
	double tau = 1;
	double Vm = 1;
	double w = 2 * PI * 10;
	double f = 10;
	return (-1 / tau) * v + (1 / tau) * Vm * cos(2 * PI * f * t);
}
//void odeFunc_mck(const double t, const double Y[], double dYdt[])
//{
//	double m = 1;
//	double c = 7;
//	double k = 6.9;
//	double f = 5;
//
//	double Fin = 2 * cos(2 * PI * f * t);  //y2dot
//
//	dYdt[0] = Y[1];
//
//	// zdot= (-k*Y - c*Z + Fin)/m;
//	dYdt[1] = (-k * Y[0] - c * Y[1] + Fin);
//}

//Matlab ODE-IVP problem compare
void odeFunc_mck(const double t, const double Y[], double dYdt[])
{
	double m = 0.5;
	double c = 0.16;
	double L = 1.2;
	double g = 9.8;

	double Fin = (-c * Y[1] / m) - (g * sin(Y[0])) / L; //y2dot

	dYdt[0] = Y[1];

	// zdot= (-k*Y - c*Z + Fin)/m;
	dYdt[1] = Fin;
}

int main(int argc, char* argv[])
{

	/*-------------------------------------------------------------------*/
	// Single of 1st Order ODE
	/*-------------------------------------------------------------------*/

	//Parameter Definitions
	double a = 0;
	double b = 0.1;
	double h = 0.001;
	unsigned int N = (b - a) / h + 1;
	double y_EU[200] = { 0 };				//Cannot use y_EU[N]
	double y_EM[200] = { 0 };
	double y_RK2[200] = { 0 };
	double y_RK3[200] = { 0 };
	double y_RK4[200] = { 0 };

	// Initial value
	double v0 = 0;

	// ODE solver
	odeEU(odeFunc_rc, y_EU, a, b, h, v0);
	odeEM(odeFunc_rc, y_EM, a, b, h, v0);

	// Exercise 1: Create a general form for RK2
	odeRK2(odeFunc_rc, y_RK2, a, b, h, v0);

	//RK3
	odeRK3(odeFunc_rc, y_RK3, a, b, h, v0);
	// Exercise 2: Create the standard form  for RK4
	odeRK4(odeFunc_rc, y_RK4, a, b, h, v0);

	// Print outputs
	printf("/*-----------------------*/\n");
	printf("/ Single of 1st Order ODE /\n");
	printf("/*-----------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++)
		printf("t= %f\tyEU= %f\tyEM= %f\tyRK2= %f\tRK3= %f\tyRK4= %f\n", a + i * h, y_EU[i], y_EM[i], y_RK2[i],y_RK3[i], y_RK4[i]);
	printf("\n");

	/*-------------------------------------------------------------------*/
	// 2nd Order ODE : MCK example
	/*-------------------------------------------------------------------*/

	//Parameter Definitions
	double t0 = 0;
	double tf = 9;
	h = 0.1;
	N = (tf - t0) / h + 1;
	double y[200] = { 0 };
	double v[200] = { 0 };
	double y1[200] = { 0 };
	double v1[200] = { 0 };
	// Initial values
	double y0 = PI/2;
	v0 = 0.2;
	// ODE solver: RK2
	sys2RK2(odeFunc_mck, y, v, t0, tf, h, y0, v0);

	// Exercise 3: Create the standard form  for RK4 for 2nd order	
	sys2RK4(odeFunc_mck, y1, v1, t0, tf, h, y0, v0);

	// Print outputs
	printf("/*---------------------------*/\n");
	printf("/ 2nd Order ODE : MCK example /\n");
	printf("/*---------------------------*/\n");
	printf(" - Total number of data N=%d \n", N);
	for (int i = 0; i < N; i++)
		printf("|RK2| t= %f\ty= %f\tv= %f\t |RK4| t= %f\ty= %f\tv= %f\n", t0 + i * h, y[i], v[i], t0 + i * h, y1[i], v1[i]);
	printf("\n\n");

	system("pause");
	return 0;
}