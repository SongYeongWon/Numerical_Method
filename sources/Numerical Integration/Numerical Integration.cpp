/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : [Yeong-Won Song]
Created          : 17-05-2021
Modified         : 19-05-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]Integration_student.cpp
-------------------------------------------------------------------------------*/

#include "../../include/myNM.h"
#define Assignment	7		// enter your assignment number
#define eval		0		// set 0
// Integration using rectangular method for discrete data inputs


double myFunc(const double x) {
	//return sqrt(1 - x * x);
	double w = 10*pow(10, 3);
	double E = 200*pow(10, 6);
	double I = 2.1*pow(10, -4);
	double x1 = pow(x, 2);
	double k = (-w * 0.5 * x1);
	return  pow(k,2) / (2 * E * I);
}

int main(int argc, char* argv[])
{
	// PART 1. Integration from Datasets
	printf("\n**************************************************");
	printf("\n        PART 1. Integration from Datasets         ");
	printf("\n**************************************************\n");

	double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
	double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
	int M = sizeof(x) / sizeof(x[0]);

	double I_rect = IntegrateRect(x, y, M);
	double I_trapz = trapz(x, y, M);
	double l_midpoint = integralMid(x, y, M);
	printf("I_rect  = %f\n", I_rect);
	printf("I_trapz = %f\n", I_trapz);
	printf("interalMid = %f\n", l_midpoint);

	// PART 2. Integration from a Function
	printf("\n**************************************************");
	printf("\n        PART 2. Integration from a Function       ");
	printf("\n**************************************************\n");

	double a = 0;
	double b = 1;
	int n = 300;
	double I_simpson13 = integral(myFunc, a, b, n);
	double l_simson38 = integral38(myFunc, a, b, n);
	printf("Simpson13 = %f\n", I_simpson13);
	printf("Simpson38 = %f\n", l_simson38);
	system("pause");
	return 0;
}

