#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fenv.h>


double factorial(int n, unsigned int* flops){
    if(n == 0)
        return 1;
    else{
        *flops += 1;
        return n * factorial(n - 1, flops);
    }
}

double exponentiation(int base, int exponent, unsigned int* flops){
    if(exponent == 0)
        return 1;
    else{
        *flops += 1;
        return base * exponentiation(base, exponent - 1, flops);
    }
}
double somatorio(double* pi, double tolerance, unsigned int* n, double* aproximate_absolut_error, unsigned int* ulps, unsigned int* flops){
    double last_value = 0;
// inicial erro absoluto
    printf("ABS INICIAL: %.15e\n", *aproximate_absolut_error);

    while( tolerance >= *aproximate_absolut_error){
        *pi += (exponentiation(2, *n, flops) * exponentiation(factorial(*n, flops), 2, flops)) / factorial(2* *n + 1, flops);
        *n += 1;
        if(*n > 2){
        *aproximate_absolut_error = fabs(*pi - last_value);
        }
        last_value = *pi;
        printf("ABS: %.15e\n", *aproximate_absolut_error);
    }
}

void pi_calc(double* pi, double tolerance, unsigned int* n, double* aproximate_absolut_error, unsigned int* ulps, unsigned int* flops){

    fesetround(FE_DOWNWARD);
    somatorio(&pi[0], tolerance, n, aproximate_absolut_error, ulps, flops);
    fesetround(FE_UPWARD);
    *aproximate_absolut_error = 0;
    *n = 1;
    somatorio(&pi[1], tolerance, n, aproximate_absolut_error, ulps, flops);
    pi[0] *= 2;
    pi[1] *= 2;
}


int main(int argc, char *argv[]){
/*
IMPRIMIR: 
Número de interações
Erro absoluto aproximado obtido
Erro absoluto exato
Apróximado de pi arredondodamento para baixo
Apróximado de pi arredondodamento para cima
Diferença de ULPs entre arredoamento p/baixo e p/cima
Número de operações de ponto flutuante
*/
    double tolerance = 0;
    unsigned int n = 1;
    double aproximate_absolut_error = 0;
    double exact_absolut_error = 0;
    // double pi_down = 0;
    // double pi_up = 0;
    double pi[2] = {0, 0};
    unsigned int ulps = 0;
    unsigned int flops = 0; 

    if(argc > 1){
        printf("INPUT: %s!\n", argv[1]);
    }else{
        printf("No input detected!\n");
        return 1;
    }

    tolerance = strtod(argv[1], NULL);
    printf("Tolerance: %.15e\n", tolerance);
    pi_calc(pi, tolerance, &n, &aproximate_absolut_error, &ulps, &flops);
    exact_absolut_error = abs(M_PI - pi[1]);


    // printf("%.15e\n", tolerance);
    printf("%u\n", n);
    printf("%.15e\n", aproximate_absolut_error);
    printf("%.15e\n", exact_absolut_error);
    printf("%.15e\n", pi[0]); // pi_down
    printf("%.15e\n", pi[1]); // pi_up
    // printf("%u\n", ulps);
    printf("%u\n", flops);




    printf("%.15e\n", M_PI);

    return 0;
}
