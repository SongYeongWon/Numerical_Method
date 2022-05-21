/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Song Yeong Won
Created          : 03-05-2021
Modified         : 25-06-2021
Language/ver     : C++ in MSVS2019

Description      : [Tutorial]curve_fitting.cpp
-------------------------------------------------------------------------------*/

#include "../../include/myNM.h"
#define Assignment	10		// enter your assignment number
#define eval		0		// set 0

int main(int argc, char* argv[])
{
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";
#if eval
path += "eval/";
#endif 

	Matrix X_DATA = txt2Mat(path, "X data");
	Matrix Y_DATA = txt2Mat(path, "Y data");
	Matrix Z_DATA = txt2Mat(path, "Z data");
	
	printf("---------------------------------------------------------------------------------\n");
	printf("                    Linenar Least Square Regression            \n");
	printf("---------------------------------------------------------------------------------\n");
	planeFit(X_DATA,Y_DATA,Z_DATA);

	printf("---------------------------------------------------------------------------------\n");
	printf("                    Linear Spline Interpolation                   \n");
	printf("---------------------------------------------------------------------------------\n");
	InterpolateLinear(X_DATA, Y_DATA, Z_DATA);
	system("pause");
	return 0;
}