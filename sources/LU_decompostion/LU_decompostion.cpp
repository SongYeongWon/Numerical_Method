/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Song-Yeong-Won
Created          : 26-03-2018
Modified         : 02-04-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	3		// enter your assignment number
#define eval		0		// set 0

#include "../../include/myNM.h"

int main(int argc, char* argv[])
{
	/*	 [�� DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
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
	Matrix matA = txt2Mat(path, "prob1_matA");
	Matrix vecb = txt2Mat(path, "prob1_vecb");
	Matrix matU = createMat(matA.rows,matA.cols);	 //matA �� ���� ũ���� ������� matrix U ����
	Matrix vecUb = createMat(vecb.rows,vecb.cols);   //vecb �� ���� ũ���� ������� matrix vecUb ���� 
													 // vecUb �� matrix A �� matrix U ���� ������� �� vecb ���� �ǹ�
	Matrix matL = createMat(matA.rows, matA.cols);
	Matrix matP = createMat(matA.rows, matA.cols);
	Matrix X = createMat(vecb.rows, vecb.cols);		//solveLU �Լ� ���� ����� ����ϱ� ���� matrix X 
	Matrix invA = createMat(matA.rows, matA.cols);	//A�� ����� ���� �ϱ� ���� Matrix invA
	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	printf("--------------------------------------------------------\n");
	printf("              Guass elimination method \n");
	printf("--------------------------------------------------------\n");
	//gaussEilm(matA, vecb, matU, vecUb);			//guass elimination �Լ� ����
	
	printf("--------------------------------------------------------\n");
	printf("              LUdecompostion  method \n");
	printf("--------------------------------------------------------\n");
	LUdecomp(matA, matL, matU, matP);				//LUdecomposition �Լ� ����
	
	
	printf("--------------------------------------------------------\n");
	printf("             solveLU by fwd,back  method \n");
	printf("--------------------------------------------------------\n");

	
	X=solveLU(matL, matU, matP, vecb);			//solveLU �Լ� ����
	printMat(X, "Final solution X = \n");			//solveLU ���� ��� ���

	printf("--------------------------------------------------------\n");
	printf("                    Inverse matrix  \n");
	printf("--------------------------------------------------------\n");
	//inv(matA, invA);								//invA �Լ� ����


	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA);		freeMat(vecb);
	
	system("pause");
	return 0;
}