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
	/*	 [※ DO NOT EDIT IT !!!]   Resources file path setting for evaluation	*/
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


	Matrix matU = createMat(matA.cols,matA.rows); //matA 와 같은 크기의 행과열의 matrix U 생성

	Matrix vecUb = createMat(vecb.rows,vecb.cols); //vecb 와 같은 크기의 행과열의 matrix vecUb 생성 
												   // vecUb 는 matrix A 를 matrix U 으로 만들었을 때 vecb 값을 의미


	/*==========================================================================*/
	/*					Apply your numerical method algorithm					*/
	/*==========================================================================*/

	gaussEilm(matA, vecb, matU, vecUb); //guass elimination 함수 실행

	

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