
#include "../../include/myNM.h"
//¼öÁ¤
void main() {

    /************      Variables declaration & initialization      ************/

    float tol = 0.00001;
    // set the bound [a ,b]
    float a0 = 0;
    float b0 = 3;
    float x0 = 1;
    double BM_result, NM_result, HBM_result; //BM_result is result of bisectionNL, NM_result is result of newtonRaphson

    /************      Test  Functions & Show Output            ************/
    printf("------------------------------------------------------------------------------------\n");
    printf("         Bisection Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");

    printf("Bisection Method:\n");
    BM_result = bisectionNL(a0, b0, tol);

    printf("Final Solution: %f \t", BM_result);
    printf("\n");

    printf("------------------------------------------------------------------------------------\n");
    printf("         newtonRaphson Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");

    printf("newtonRaphson Method:\n");
    NM_result = newtonRaphson(x0, tol);

    printf("Final Solution: %f \t", NM_result);
    printf("\n");

    printf("------------------------------------------------------------------------------------\n");
    printf("         hybrid Method Results             \n");
    printf("------------------------------------------------------------------------------------\n");

    printf("hybrid Method:\n");

    HBM_result = hybridmethod(x0, tol);

    printf("Final Solution: %f \t", HBM_result);
    printf("\n");
    system("pause");
}






