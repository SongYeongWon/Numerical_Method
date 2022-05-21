/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [Yeong-Won Song]
Created          : 17-05-2021
Modified         : 31-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Integration_student.cpp
-------------------------------------------------------------------------------*/

#include "../../include/myNM.h"
#define Assignment	8		// enter your assignment number
#define Pi	3.14159265358979323846264338327950288419716939937510
#define Eu  0
#define Em  1
#define eval		0		// set 0

double myFunc(const double x, const double y) {
	double tau = 1;
	double Vm = 1;
	double w = 2 * Pi * 10;
	double f = 10;
	return (-1 / tau) * y + (1 / tau) * Vm * cos(2 * Pi * f * x);
}

int main(int argc, char* argv[])
{
	printf("\n**************************************************");
	printf("\n                Euler's Method        ");
	printf("\n**************************************************\n");

	double h = 0.001;
	double x_array[101] = { 0 };
	double y_array[101] = { 0 };
	double t0 = 0;
	double tf = 0.1;
	odeEU(myFunc, x_array, y_array, t0, tf, h);
	printf("x \n");
	printarr(x_array, 101);
	printf("y \n");
	printarr(y_array, 101);

	printf("\n**************************************************");
	printf("\n               Modified Euler's Method        ");
	printf("\n**************************************************\n");
	double x1_array[101] = { 0 };
	double y1_array[101] = { 0 };
	odeEM(myFunc, x1_array, y1_array, t0, tf, h);
	printf("y \n");
	printarr(y1_array, 101);

	//ode(myFunc, x_array, y_array, t0, tf, h, Em);


	system("pause");
	return 0;
}

