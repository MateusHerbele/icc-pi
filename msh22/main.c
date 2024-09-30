#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <fenv.h>


double factorial(int n, long long int* flops){
    if(n == 0)
        return 1;
    else{
        *flops += 1;
        return n * factorial(n - 1, flops);
    }
}

double exponentiation(int base, int exponent, long long int* flops){
    if(exponent == 1)
        return base;
    else{
        *flops += 1;
        return base * exponentiation(base, exponent - 1, flops);
    }
}
double somatorio(double* pi, double tolerance, unsigned int* n, double* aproximate_absolut_error, unsigned int* ulps, long long int* flops){
    double last_value = 0;
    double exp = 1;
    double fact_up = 1; // Fatorial de cima da fração
    double fact_down = 1; // Fatorial de baixo da fração
    double up_part = 0; // Parte de cima da fração
    
    // tentar fazer as variáveis já estarem todas alinhadas a partir da execução 1 (2)

    // inicial erro absoluto
    //printf("ABS INICIAL: %.15e\n", *aproximate_absolut_error);
    do{
        if(*n == 0){
            *pi = 2.0;
        }
        else{
            
            up_part = (exp * (fact_up*fact_up));
            *pi += 2.0 * (up_part / fact_down);
            *flops += 3;
        }
        *aproximate_absolut_error = fabs(*pi - last_value); // talvez conte como 1 operação de ponto flutuante
        *n += 1;
        last_value = *pi;
        // printf("Iteração: %d, PI: %.15e, ABS: %.15e\n", *n, *pi, *aproximate_absolut_error);
        exp *= 2;
        fact_up = fact_up * (*n);
        fact_down = fact_down * ((2 * *n + 1) * (2 * *n));
        *flops += 3; // considerando abe
    }while(*aproximate_absolut_error >= tolerance);
    
}

void pi_calc(double* pi, double tolerance, unsigned int* n, double* aproximate_absolut_error, unsigned int* ulps, long long int* flops){

    fesetround(FE_DOWNWARD);
    somatorio(&pi[0], tolerance, n, aproximate_absolut_error, ulps, flops);
    fesetround(FE_UPWARD);
    *aproximate_absolut_error = 0;
    *n = 0;
    somatorio(&pi[1], tolerance, n, aproximate_absolut_error, ulps, flops);
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
    int n = 0; 
    double aproximate_absolut_error = 0;
    double exact_absolut_error = 0;
    double pi[2] = {0, 0};
    int64_t *ptr_pi_down = NULL;
    int64_t *ptr_pi_up = NULL;
    int64_t *ptr_abe = NULL;
    int64_t *ptr_ebe = NULL;
    
    int ulps = 0;
    long long int flops = 0; 

    if(argc > 1){
        //printf("INPUT: %s!\n", argv[1]);
    }else{
        printf("No input detected!\n");
        return 1;
    }

    tolerance = strtod(argv[1], NULL);
    // printf("Tolerance: %.15e\n", tolerance);
    pi_calc(pi, tolerance, &n, &aproximate_absolut_error, &ulps, &flops);
    exact_absolut_error = fabs(M_PI - pi[1]);


    ptr_pi_down = (int64_t *) &pi[0];
    ptr_pi_up = (int64_t *) &pi[1];
    ptr_abe = (int64_t *) &aproximate_absolut_error;
    ptr_ebe = (int64_t *) &exact_absolut_error;
    ulps = *ptr_pi_up - *ptr_pi_down;

    printf("%d\n", n);
    printf("%.15e %lx\n", aproximate_absolut_error, *ptr_abe);

    printf("%.15e %lx\n", exact_absolut_error, *ptr_ebe);
    
    printf("%.15e %lx\n", pi[0], *ptr_pi_down); // pi_down

    printf("%.15e %lx\n", pi[1], *ptr_pi_up); // pi_up

    printf("%d\n", ulps);
    printf("%lld\n", flops);

    return 0;
}
