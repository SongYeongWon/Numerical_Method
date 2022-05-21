/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Song-Yeong-Won
Created          : 26-03-2018
Modified         : 26-04-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	4		// enter your assignment number
#define eval		0		// set 0

#include "../../include/myNM.h"

int main(int argc, char* argv[])
{
	/*	 [¡Ø DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
	std::string path = "C:/NM_data_2021/Assignment" + std::to_string(Assignment) + "/";

#if eval
	path += "eval/";
#endif 

	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/
	Matrix matA = txt2Mat(path, "matA");
	Matrix matC = txt2Mat(path, "matC");
	Matrix matQ = createMat(matA.rows, matA.rows);	
	Matrix matR = createMat(matA.rows, matA.cols);

	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	
	printf("--------------------------------------------------------\n");
	printf("              QR decomposition method  \n");
	printf("--------------------------------------------------------\n");

	QRdecomp(matA, matQ, matR);
;
	printf("--------------------------------------------------------\n");
	printf("                      eigenvalue   \n");
	printf("--------------------------------------------------------\n");
	
	eigenvalue(matA, matQ, matR);



	printf("--------------------------------------------------------\n");
	printf("                      Cond Number \n");
	printf("--------------------------------------------------------\n");
	
	double CondNumberA, CondNumberC = 0;
	CondNumberA = Cond(matA);
	CondNumberC = Cond(matC);
	printf("CondNumberA = %f\n", CondNumberA);
	printf("CondNumberC = %f\n", CondNumberC);
	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA);		
	
	/*double tol = 0.00001;
	// set the bound [a ,b]
	double a0 = 0;
	double b0 = 3;
	double x0 =1;
	double BM_result, NM_result;

	printf("------------------------------------------------------------------------------------\n");
	printf("         bisection method results             \n");
	printf("------------------------------------------------------------------------------------\n");

	printf("bisection method:\n");
	BM_result = bisectionNL(a0, b0, tol);

	printf("final solution: %f \t", BM_result);
	printf("\n");

	printf("------------------------------------------------------------------------------------\n");
	printf("         newtonRaphson method results             \n");
	printf("------------------------------------------------------------------------------------\n");

	printf("newtonRaphson method:\n");
	NM_result = newtonRaphson(x0, tol);

	printf("final solution: %f \t", NM_result);

	printf("\n");*/
	system("pause");
	return 0;
}