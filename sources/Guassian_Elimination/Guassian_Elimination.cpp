/*-------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Song-Yeong-Won
Created          : 26-03-2018
Modified         : 29-03-2021
Language/ver     : C++ in MSVS2019

Description      : NM_MainTemplate.cpp
-------------------------------------------------------------------------------*/

#define Assignment	2		// enter your assignment number
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


	Matrix matU = createMat(matA.cols,matA.rows); //matA �� ���� ũ���� ������� matrix U ����

	Matrix vecUb = createMat(vecb.rows,vecb.cols); //vecb �� ���� ũ���� ������� matrix vecUb ���� 
												   // vecUb �� matrix A �� matrix U ���� ������� �� vecb ���� �ǹ�


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	gaussEilm(matA, vecb, matU, vecUb); //guass elimination �Լ� ����

	

	/*==========================================================================*/
	/*							  Print your results							*/
	/*==========================================================================*/
	
	

	/*==========================================================================*/
	/*							  Deallocate memory 							*/
	/*==========================================================================*/
	freeMat(matA);		freeMat(vecb);
	
	system("pause");
	return 0;
}